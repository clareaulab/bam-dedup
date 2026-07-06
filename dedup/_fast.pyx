# cython: boundscheck=False, wraparound=False, cdivision=True, language_level=3
"""
Cython acceleration for the profiled hot spots of bam-dedup.

picardlike (PCR/optical duplicate marking):
  1. sum_of_base_qualities  -- tight loop over the base-quality array
  2. parse_location         -- Illumina read-name tile/x/y parsing
  3. optical_flags_graph    -- union-find optical-duplicate clustering (O(n^2))

fgbiolike (UMI molecular-consensus deduplication):
  4. consensus_call_molecule -- per-base Kahan-summed consensus likelihood loop

Each function is behavior-identical to the pure-Python reference in the
corresponding module; only the mechanics are lowered to C.  Import is optional --
the callers fall back to pure Python if this module is not built.
"""

from libc.stdlib cimport malloc, free
from libc.math cimport log, log1p, exp, expm1, floor, fabs
from cpython.mem cimport PyMem_Malloc, PyMem_Free


# ---------------------------------------------------------------------------
# 1. sum of base qualities >= 15  (htsjdk getSumOfBaseQualities)
# ---------------------------------------------------------------------------
cpdef int sum_of_base_qualities(const unsigned char[:] quals):
    cdef Py_ssize_t i, n = quals.shape[0]
    cdef int s = 0
    cdef unsigned char q
    for i in range(n):
        q = quals[i]
        if q >= 15:
            s += q
    return s


# ---------------------------------------------------------------------------
# 2. read-name location parsing (ReadNameParser default regex behavior)
#    Require exactly 5 or 7 ':'-separated fields; take last three as tile/x/y,
#    each parsed with rapidParseInt semantics (leading [-]digits, stop at
#    first non-digit).  Returns (tile, x, y) or None.
# ---------------------------------------------------------------------------
cdef inline long _rapid_parse(const char* s, Py_ssize_t start, Py_ssize_t end, bint* ok):
    cdef Py_ssize_t i = start
    cdef bint neg = False
    cdef long val = 0
    cdef bint has = False
    cdef char c
    ok[0] = True
    if i < end and s[i] == b'-':
        neg = True
        i += 1
    while i < end:
        c = s[i]
        if c >= b'0' and c <= b'9':
            val = val * 10 + (c - 48)
            has = True
            i += 1
        else:
            break
    if not has:
        ok[0] = False
        return 0
    return -val if neg else val


def parse_location(str read_name):
    cdef bytes b = read_name.encode('ascii', 'replace')
    cdef const char* s = b
    cdef Py_ssize_t n = len(b)
    cdef Py_ssize_t i
    cdef int nfields = 1
    for i in range(n):
        if s[i] == b':':
            nfields += 1
    if nfields != 5 and nfields != 7:
        return None

    # locate the boundaries of the last three fields
    cdef Py_ssize_t last_colon = -1, second_colon = -1, third_colon = -1
    cdef int seen = 0
    for i in range(n - 1, -1, -1):
        if s[i] == b':':
            seen += 1
            if seen == 1:
                last_colon = i
            elif seen == 2:
                second_colon = i
            elif seen == 3:
                third_colon = i
                break
    if last_colon < 0 or second_colon < 0 or third_colon < 0:
        return None

    cdef bint ok = True
    cdef long tile = _rapid_parse(s, third_colon + 1, second_colon, &ok)
    if not ok:
        return None
    cdef long x = _rapid_parse(s, second_colon + 1, last_colon, &ok)
    if not ok:
        return None
    cdef long y = _rapid_parse(s, last_colon + 1, n, &ok)
    if not ok:
        return None
    return (tile, x, y)


# ---------------------------------------------------------------------------
# 3. optical-duplicate union-find graph clustering.
#    Inputs are parallel arrays over the duplicate-set members.  tile == -1
#    means "no location" (skipped).  Returns a Python list[bool] flagging
#    optical duplicates, identical to OpticalDuplicateFinder's graph path.
# ---------------------------------------------------------------------------
cdef int _find(int* parent, int a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]
        a = parent[a]
    return a


def optical_flags_graph(list xs, list ys, list tiles, list rgs,
                        int keeper_index, int dist):
    cdef int n = len(xs)
    cdef bint* flags = <bint*> malloc(n * sizeof(bint))
    cdef int* parent = <int*> malloc(n * sizeof(int))
    cdef int* ax = <int*> malloc(n * sizeof(int))
    cdef int* ay = <int*> malloc(n * sizeof(int))
    cdef int* atile = <int*> malloc(n * sizeof(int))
    cdef int* arg = <int*> malloc(n * sizeof(int))
    if not flags or not parent or not ax or not ay or not atile or not arg:
        raise MemoryError()

    cdef int i, j, ra, rb
    for i in range(n):
        flags[i] = False
        parent[i] = i
        ax[i] = xs[i]
        ay[i] = ys[i]
        atile[i] = tiles[i]
        arg[i] = rgs[i]

    # Union within same (readGroup, tile) group and within pixel distance.
    # Grouping is by connectivity only, so a direct O(n^2) scan restricted to
    # matching (rg, tile) yields the identical clustering as the tile-bucketed
    # version, without needing a hash map.
    for i in range(n):
        if atile[i] == -1:
            continue
        for j in range(i + 1, n):
            if atile[j] == -1:
                continue
            if arg[i] == arg[j] and atile[i] == atile[j]:
                if _iabs(ax[i] - ax[j]) <= dist and _iabs(ay[i] - ay[j]) <= dist:
                    ra = _find(parent, i)
                    rb = _find(parent, j)
                    if ra != rb:
                        parent[ra] = rb

    # Cluster representative assignment (min-x, then min-y), keeper protected.
    # cluster_rep[root] = current representative index; -1 = unset.
    cdef int* cluster_rep = <int*> malloc(n * sizeof(int))
    for i in range(n):
        cluster_rep[i] = -1

    cdef int keeper_cluster = -1
    if keeper_index >= 0:
        keeper_cluster = _find(parent, keeper_index)
        cluster_rep[keeper_cluster] = keeper_index

    cdef int c, rep
    cdef bint in_keeper_cluster
    for i in range(n):
        c = _find(parent, i)
        if cluster_rep[c] != -1 and i != keeper_index:
            rep = cluster_rep[c]
            in_keeper_cluster = (keeper_index >= 0 and c == keeper_cluster)
            if (not in_keeper_cluster and
                    (ax[i] < ax[rep] or (ax[i] == ax[rep] and ay[i] < ay[rep]))):
                flags[rep] = True
                cluster_rep[c] = i
            else:
                flags[i] = True
        else:
            cluster_rep[c] = i

    result = [False] * n
    for i in range(n):
        result[i] = bool(flags[i])

    free(flags); free(parent); free(ax); free(ay); free(atile); free(arg)
    free(cluster_rep)
    return result


cdef inline int _iabs(int v):
    return v if v >= 0 else -v


# ===========================================================================
# 4. Molecular consensus calling (fgbio VanillaUmiConsensusCaller, fragment path)
#    Mirrors ConsensusCaller.ConsensusBaseBuilder including Kahan (compensated)
#    summation of per-base log-likelihoods, so it is numerically identical to the
#    pure-Python reference (and to fgbio, modulo libm's last ULP).
# ===========================================================================
cdef double _LN10 = log(10.0)
cdef double _LN2 = log(2.0)
cdef double _LOG_FOUR_THIRDS = log(4.0) - log(3.0)
cdef double _NEG_INF = -1e308 * 10
cdef double _EPSILON = 2.0 ** -52
cdef double _MAX_VALUE_AS_LOG = _LN10 * 93.0 / -10.0
cdef int _PHRED_MIN = 2
cdef int _PHRED_MAX = 93
cdef double _PRECISION = 0.001
cdef unsigned char _A = 65, _C = 67, _G = 71, _T = 84, _N = 78


cdef inline double _log1pexp(double v) nogil:
    if v <= -37:
        return exp(v)
    elif v <= 18:
        return log1p(exp(v))
    elif v <= 33.3:
        return v + exp(-v)
    else:
        return v


cdef inline double _log1mexp(double v) nogil:
    if v <= _LN2:
        return log(-expm1(-v))
    else:
        return log1p(-exp(-v))


cdef inline double _lor(double a, double b) nogil:
    cdef double t
    if a == _NEG_INF:
        return b
    if b == _NEG_INF:
        return a
    if b < a:
        t = a; a = b; b = t
    return a + _log1pexp(b - a)


cdef inline double _a_or_not_b(double a, double b) nogil:
    if b == _NEG_INF:
        return a
    if a == b:
        return _NEG_INF
    return a + _log1mexp(a - b)


cdef inline double _lnot(double a) nogil:
    if 0.0 < a:
        return _NEG_INF
    return _a_or_not_b(0.0, a)


cdef inline double _prob_error_two_trials(double a, double b) nogil:
    cdef double t, term1, term2
    if a < b:
        t = a; a = b; b = t
    if a - b >= 6:
        return a
    term1 = _lor(a, b)
    term2 = _LOG_FOUR_THIRDS + a + b
    return _a_or_not_b(term1, term2)


cdef inline double _lor4(double x0, double x1, double x2, double x3) nogil:
    cdef double vals[4]
    vals[0] = x0; vals[1] = x1; vals[2] = x2; vals[3] = x3
    cdef int min_i = 0
    cdef double min_v = vals[0]
    cdef int i
    for i in range(1, 4):
        if vals[i] < min_v:
            min_v = vals[i]; min_i = i
    cdef double s = min_v
    for i in range(4):
        if i != min_i:
            s = _lor(s, vals[i])
    return s


cdef inline int _phred_from_logprob(double lp) nogil:
    if lp < _MAX_VALUE_AS_LOG:
        return _PHRED_MAX
    return <int>floor(-10.0 * (lp / _LN10) + _PRECISION)


cdef inline int _cap_phred(int q) nogil:
    if q < _PHRED_MIN:
        return _PHRED_MIN
    if q > _PHRED_MAX:
        return _PHRED_MAX
    return q


def consensus_call_molecule(list p_err_third, list p_truth, double ln_pre,
                            list capped, int length, int min_reads, int min_consensus_q,
                            bases_out, quals_out, depths_out, errors_out):
    """Fill the output arrays for one molecule's consensus read.

    p_err_third, p_truth : length-127 lists of log-probabilities (per input qual)
    ln_pre               : ln(pre-UMI error rate)
    capped               : list of source reads (objects with .bases bytes, .quals, .length)
    *_out                : preallocated arrays (bytearray, array('b'/'B'), array('h'), array('h'))
    """
    cdef int n = len(capped)
    cdef int i, pos, r

    cdef double pet[127]
    cdef double ptt[127]
    for i in range(127):
        pet[i] = p_err_third[i]
        ptt[i] = p_truth[i]

    cdef int* lengths = <int*>PyMem_Malloc(n * sizeof(int))
    cdef int* offsets = <int*>PyMem_Malloc(n * sizeof(int))
    cdef int total = 0
    for r in range(n):
        lengths[r] = capped[r].length
        offsets[r] = total
        total += lengths[r]
    cdef unsigned char* basebuf = <unsigned char*>PyMem_Malloc(total * sizeof(unsigned char))
    cdef unsigned char* qualbuf = <unsigned char*>PyMem_Malloc(total * sizeof(unsigned char))

    cdef bytes bb
    cdef object qq
    cdef int off, L, k
    for r in range(n):
        bb = capped[r].bases
        qq = capped[r].quals
        off = offsets[r]
        L = lengths[r]
        for k in range(L):
            basebuf[off + k] = bb[k]
            qualbuf[off + k] = qq[k]

    cdef unsigned char[:] bo = bases_out
    cdef signed char[:] qo = quals_out
    cdef short[:] do = depths_out
    cdef short[:] eo = errors_out

    cdef double ll0, ll1, ll2, ll3, c0, c1, c2, c3
    cdef double y, t, term, pt, pe
    cdef int obs0, obs1, obs2, obs3, depth
    cdef unsigned char base
    cdef int q
    cdef double ll_sum, max_v, max_posterior, p_cons_err, p
    cdef int max_i, assigned, raw_qual, errors
    cdef unsigned char raw_base

    for pos in range(length):
        ll0 = 0.0; ll1 = 0.0; ll2 = 0.0; ll3 = 0.0
        c0 = 0.0; c1 = 0.0; c2 = 0.0; c3 = 0.0
        obs0 = 0; obs1 = 0; obs2 = 0; obs3 = 0

        for r in range(n):
            if lengths[r] > pos:
                base = basebuf[offsets[r] + pos]
                if base == _N:
                    continue
                q = qualbuf[offsets[r] + pos]
                pe = pet[q]
                pt = ptt[q]
                term = pt if base == _A else pe
                y = term - c0; t = ll0 + y; c0 = (t - ll0) - y; ll0 = t
                term = pt if base == _C else pe
                y = term - c1; t = ll1 + y; c1 = (t - ll1) - y; ll1 = t
                term = pt if base == _G else pe
                y = term - c2; t = ll2 + y; c2 = (t - ll2) - y; ll2 = t
                term = pt if base == _T else pe
                y = term - c3; t = ll3 + y; c3 = (t - ll3) - y; ll3 = t
                if base == _A: obs0 += 1
                elif base == _C: obs1 += 1
                elif base == _G: obs2 += 1
                elif base == _T: obs3 += 1

        depth = obs0 + obs1 + obs2 + obs3

        ll_sum = _lor4(ll0, ll1, ll2, ll3)
        max_v = -1.7976931348623157e308
        max_i = -1
        assigned = 0
        if (assigned == 0) or (ll0 > max_v):
            max_v = ll0; max_i = 0; assigned = 1
        if ll1 > max_v:
            max_v = ll1; max_i = 1
        elif fabs(ll1 - max_v) <= _EPSILON:
            max_i = -1
        if ll2 > max_v:
            max_v = ll2; max_i = 2
        elif fabs(ll2 - max_v) <= _EPSILON:
            max_i = -1
        if ll3 > max_v:
            max_v = ll3; max_i = 3
        elif fabs(ll3 - max_v) <= _EPSILON:
            max_i = -1

        if max_i == -1:
            raw_base = _N
            raw_qual = _PHRED_MIN
        else:
            max_posterior = max_v - ll_sum
            p_cons_err = _lnot(max_posterior)
            p = _prob_error_two_trials(ln_pre, p_cons_err)
            raw_qual = _cap_phred(_phred_from_logprob(p))
            if max_i == 0: raw_base = _A
            elif max_i == 1: raw_base = _C
            elif max_i == 2: raw_base = _G
            else: raw_base = _T

        if raw_base == _N:
            errors = depth
        elif raw_base == _A:
            errors = depth - obs0
        elif raw_base == _C:
            errors = depth - obs1
        elif raw_base == _G:
            errors = depth - obs2
        else:
            errors = depth - obs3

        if depth < min_reads:
            bo[pos] = _N; qo[pos] = 0
        elif raw_qual < min_consensus_q:
            bo[pos] = _N; qo[pos] = 2
        else:
            bo[pos] = raw_base; qo[pos] = <signed char>raw_qual

        do[pos] = <short>(depth if depth <= 32767 else 32767)
        eo[pos] = <short>(errors if errors <= 32767 else 32767)

    PyMem_Free(lengths)
    PyMem_Free(offsets)
    PyMem_Free(basebuf)
    PyMem_Free(qualbuf)
