from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("dedup_umi.pyx"),
    zip_safe=False,
)
