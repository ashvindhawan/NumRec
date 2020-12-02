from utils import *
import numc as nc 
import dumbpy as dp
dp_mat1, nc_mat1 = rand_dp_nc_matrix(3, 3, seed=0)
expect = dp_mat1 ** 2
actual = nc_mat1 ** 2
print("Should be", expect, "it was", actual)