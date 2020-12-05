from utils import *
import numc as nc 
import dumbpy as dp
dp_mat1, nc_mat1 = rand_dp_nc_matrix(64, 64, seed=0)
expect = dp_mat1 ** 1
actual = nc_mat1 ** 1
print("Should be:", expect[0])
print("it was   :", actual[0])