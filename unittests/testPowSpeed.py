from utils import *
from unittest import TestCase

class TestPow(TestCase):
    # def test_small_pow(self):
    #     # TODO: YOUR CODE HERE
    #     dp_mat, nc_mat = rand_dp_nc_matrix(64, 64, seed=0)
    #     is_correct, speed_up = compute([dp_mat, 1], [nc_mat, 1], "pow")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

    def test_medium_pow(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(512, 512, seed=0)
        is_correct, speed_up = compute([dp_mat, 500], [nc_mat, 500], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    # def test_large_pow(self):
    #     # TODO: YOUR CODE HERE
    #     dp_mat, nc_mat = rand_dp_nc_matrix(1024, 1024, seed=0)
    #     is_correct, speed_up = compute([dp_mat, 7], [nc_mat, 7], "pow")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)