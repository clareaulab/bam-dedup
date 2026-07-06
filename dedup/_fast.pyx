# cython: boundscheck=False, wraparound=False, cdivision=True, language_level=3
"""
Cython acceleration for the three hot spots identified by profiling pymarkdup:

  1. sum_of_base_qualities  -- tight loop over the base-quality array
  2. parse_location         -- Illumina read-name tile/x/y parsing
  3. optical_flags_graph    -- union-find optical-duplicate clustering (O(n^2))

Each function is behavior-identical to the pure-Python reference in pymarkdup.py;
only the mechanics are lowered to C.  Import is optional -- pymarkdup falls back
to pure Python if this module is not built.
"""

from libc.stdlib cimport malloc, free


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
