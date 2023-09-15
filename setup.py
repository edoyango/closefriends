from numpy.distutils.core import setup, Extension

ext = Extension(name='closefriends', 
                sources=['closefriends/src/closefriends.F90', 'closefriends/src/sort.cpp'],
                libraries=['gfortran'],
                extra_compile_args=['-std=c++11', '-O3'])

setup(
    name='closefriends',
    ext_modules=[ext]
)