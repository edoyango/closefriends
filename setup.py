from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
from pybind11 import get_include

ext_modules = [Pybind11Extension(
                        'closefriends',
                        sources=['src/closefriends.cpp'],
                        include_dirs=[get_include()],
                        language='c++',
                        extra_compile_args=["-O3"],
                        cpp_standard=11
                       )]

dev_requirements = [
    "scipy",
    "pytest"
]

setup(
    name='closefriends',
    version='0.0.1',
    description="A package to find point pairs in a point cloud within a fixed distance using a Cell lists data-structure.",
    author="Edward Yang",
    author_email="edward_yang_125@hotmail.com",
    ext_modules=ext_modules,
    install_requires=[
        "numpy",
        'importlib-metadata; python_version == "3.10"'
    ],
    extras_requires={
        "dev": dev_requirements
    },
    cmdclass={"build_ext": build_ext}
)
