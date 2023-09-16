from numpy.distutils.core import setup, Extension

ext = Extension('closefriends.closefriends',
                       sources=['closefriends/src/closefriends.F90', 'closefriends/src/sort.cpp'],
                       libraries=['gfortran'],
                       extra_compile_args=['-O3'])

setup(
    name='closefriends',
    packages=['closefriends'],
    ext_modules=[ext]
)