"""
bam-dedup: a fast, JVM-free reimplementation of Picard MarkDuplicates.
"""
from setuptools import find_packages, setup, Extension

# Build the Cython extension from the .pyx when Cython is available (the normal
# path, since pyproject.toml requires it at build time); otherwise fall back to
# a pre-generated .c shipped in the sdist.
try:
    from Cython.Build import cythonize
    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False

ext = ".pyx" if HAVE_CYTHON else ".c"
extensions = [Extension("dedup._fast", ["dedup/_fast" + ext])]
if HAVE_CYTHON:
    extensions = cythonize(extensions, compiler_directives={"language_level": "3"})

dependencies = ["pysam"]

setup(
    name="bam-dedup",
    version="0.2.0",
    url="https://github.com/caleblareau/bam-dedup",
    license="MIT",
    author="Caleb Lareau",
    author_email="caleb.lareau@gmail.com",
    description="Fast, JVM-free BAM deduplication: PCR/optical (Picard-like) "
                "and molecular/UMI consensus (fgbio-like).",
    long_description=__doc__,
    packages=find_packages(exclude=["tests"]),
    ext_modules=extensions,
    include_package_data=True,
    zip_safe=False,
    platforms="any",
    python_requires=">=3.8",
    setup_requires=["Cython"],
    install_requires=dependencies,
    extras_require={"test": ["pytest"]},
    entry_points={
        "console_scripts": [
            "bam-dedup = dedup.cli:main",
            "bam-consensus = dedup.fgbiolike:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Cython",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
