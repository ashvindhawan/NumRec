from utils import *
import numc as nc 
import dumbpy as dp
dp_mat1, nc_mat1 = rand_dp_nc_matrix(2, 2, seed=0)
dp_mat2, nc_mat2 = rand_dp_nc_matrix(2, 2, seed=1)
nc_mat1.set(1,1,1.0)
print(nc_mat1)