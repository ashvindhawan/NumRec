from utils import *
import numc as nc 
import dumbpy as dp

# dp_mat1, nc_mat1 = rand_dp_nc_matrix(64, 30, seed=0)
# dp_mat2, nc_mat2 = rand_dp_nc_matrix(30, 12, seed=1)

dp_mat1, nc_mat1 = rand_dp_nc_matrix(2, 4, seed=0)
dp_mat2, nc_mat2 = rand_dp_nc_matrix(4, 2, seed=1)

# expect = dp_mat1 * dp_mat2
# actual = nc_mat1 * nc_mat2
# print("Expected \n", expect[0])
# print("Actual   \n", actual[0])

dp_lst = [dp_mat1, dp_mat2]
nc_lst = [nc_mat1, nc_mat2]

is_correct, speed_up = compute(dp_lst, nc_lst, "mul")

print("Matches: ", is_correct)