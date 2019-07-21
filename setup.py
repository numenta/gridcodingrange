import sys
from setuptools import setup, Extension


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


sources = [
    'src/gridcodingrange_module.cpp',
    'src/GridUniqueness.cpp',
]


compile_args = ["-std=c++14"]
link_args = []

debug_mode = False
if debug_mode:
    compile_args += ["-O0", "-D NTA_ASSERTIONS_ON"]
else:
    compile_args += ["-g0"]


if sys.platform == "darwin":
    compile_args += ["-std=c++14", "-mmacosx-version-min=10.10"]
    link_args += ["-stdlib=libc++", "-mmacosx-version-min=10.10"]


module = Extension(
    '_gridcodingrange',
    sources=sources,
    extra_compile_args=compile_args,
    extra_link_args=link_args,
    include_dirs=['./src/external',
                  get_pybind_include(),
                  get_pybind_include(user=True)]
)

setup(name='gridcodingrange',
      version='1.0',
      description='This is a demo package',
      packages=['gridcodingrange'],
      setup_requires=["pybind11"],
      ext_modules=[module])
