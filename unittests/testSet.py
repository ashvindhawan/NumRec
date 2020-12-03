from utils import *
import numc as nc 
import dumbpy as dp
dp_mat1, nc_mat1 = rand_dp_nc_matrix(3, 4, seed=0)
dp_mat2, nc_mat2 = rand_dp_nc_matrix(2, 2, seed=1)
print(nc_mat1)
nc_mat1.set(2,2,7.3)
print(nc_mat1)