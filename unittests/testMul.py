from utils import *
import numc as nc 
import dumbpy as dp

dp_mat1, nc_mat1 = rand_dp_nc_matrix(5, 5, seed=0)
dp_mat2, nc_mat2 = rand_dp_nc_matrix(5, 5, seed=1)

expect = dp_mat1 * dp_mat2
actual = nc_mat1 * nc_mat2
print("Should be\n", expect, "\nit was\n", actual)