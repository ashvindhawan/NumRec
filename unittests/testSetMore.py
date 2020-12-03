from utils import *
import numc as nc 
import dumbpy as dp
dp_mat1, nc_mat1 = rand_dp_nc_matrix(2, 2, seed=0)
dp_mat2, nc_mat2 = rand_dp_nc_matrix(2, 2, seed=1)
nc_mat1.set(1,1,1)
dp_mat1.set(1,1,1)
print("Ours",  nc_mat1, "Correct", dp_mat1)
print("Equal?", nc_mat1==dp_mat1)
#Type Error: num args != 3, i and j not integer, or val not a float or int
#Index Error if i, j, or both are out of range
dp_mat2, nc_mat2 = rand_dp_nc_matrix(10, 10, seed=12)
nc_mat2.set(3,4,5)
dp_mat2.set(3,4,5)
nc_mat2.set(5,8,1.0)
dp_mat2.set(5,8,1.0)
nc_mat2.set(9,9,1)
dp_mat2.set(9,9,1)
print("Ours",  nc_mat2, "Correct", dp_mat2)
print("Equal?", nc_mat2==dp_mat2)

nc_mat2.set(1,1,1)