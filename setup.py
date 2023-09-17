from setuptools import setup, Extension
import pybind11

ext_modules = [Extension(
                        'closefriends',
                        sources=['src/closefriends.cpp'],
                        include_dirs=[pybind11.get_include()],
                        language='c++',
                        extra_compile_args=["-O3"]
                       )]

setup(
    name='closefriends',
    version='0.1',
    ext_modules=ext_modules
)
