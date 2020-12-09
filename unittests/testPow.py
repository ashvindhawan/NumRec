from utils import *
import numc as nc 
import dumbpy as dp
dp_mat, nc_mat = rand_dp_nc_matrix(4, 4, seed=0)
expect = dp_mat ** 10
actual = nc_mat ** 10
print("Should be:", expect[0])
print("it was   :", actual[0])

# is_correct, speed_up = compute([dp_mat, 6], [nc_mat, 6], "pow")


# print("Matches: ", is_correct)