import os
import setuptools

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname), "rt", encoding="utf-8").read()

__version__ = "1.0rc1"

install_requires = [
    "jax>=0.2.13",
    "numba>=0.53.1",
    "regex"
    ]

ext_modules = []

