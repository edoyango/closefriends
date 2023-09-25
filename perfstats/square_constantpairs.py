import scipy as sp
import seaborn as sb
import numpy as np
import closefriends as cf
from timeit import timeit
import pandas as pd

dims_iterator = range(1, 7); ndims = len(dims_iterator)
npoints_iterator = [2**i for i in range(7, 20)]; nnpoints = len(npoints_iterator)

times_sp_square = np.zeros((ndims, nnpoints))
times_cf_square = np.zeros((ndims, nnpoints))

def setup_stmt_generator(seed, dim, npoints):
    return \
f"""
import numpy as np, scipy as sp, closefriends as cf

def generateRandomPoints(dim, npoints):
    rng = np.random.default_rng({seed})
    x = rng.random(dim * npoints).reshape((npoints, dim), order="C")
    cutoff = 2.0 * npoints ** (-1.0 / dim)
    return x, cutoff

x, cutoff = generateRandomPoints({dim}, {npoints})
maxnpair = 330*{npoints}
"""

seed = 0
for i, dim in enumerate(dims_iterator):
    inner_times = []
    for j, npoints in enumerate(npoints_iterator):
        seed += 1
        setup_stmt = setup_stmt_generator(seed, dim, npoints)
        run_stmt = \
'''
tree = sp.spatial.cKDTree(x)
pairs = tree.query_pairs(cutoff, output_type='ndarray')
'''
        times_sp_square[i, j] = timeit(run_stmt, setup_stmt, number = 5)

seed = 0
for i, dim in enumerate(dims_iterator):
    inner_times = []
    for j, npoints in enumerate(npoints_iterator):
        seed += 1
        setup_stmt = setup_stmt_generator(seed, dim, npoints)
        run_stmt = "pairs = cf.query_pairs(x, cutoff, maxnpair, output_type='ndarray')"
        times_cf_square[i, j] = timeit(run_stmt, setup_stmt, number = 5)

df_sp_square = pd.DataFrame(times_sp_square, index=dims_iterator, columns=npoints_iterator)
df_cf_square = pd.DataFrame(times_cf_square, index=dims_iterator, columns=npoints_iterator)
df_sp_square.to_csv("times_sp_square.csv")
df_cf_square.to_csv("times_cf_square.csv")