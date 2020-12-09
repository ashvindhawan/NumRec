from utils import *
from unittest import TestCase

class TestPow(TestCase):

    # def test_0064x0250_pow(self):
    #     # TODO: YOUR CODE HERE
    #     dp_mat, nc_mat = rand_dp_nc_matrix(64, 64, seed=0)
    #     is_correct, speed_up = compute([dp_mat, 250], [nc_mat, 250], "pow")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

    # def test_0256x0250_pow(self):
    #     dp_mat, nc_mat = rand_dp_nc_matrix(256, 256, seed=0)
    #     is_correct, speed_up = compute([dp_mat, 250], [nc_mat, 250], "pow")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

    # def test_0512x0250_pow(self):
    #     # TODO: YOUR CODE HERE
    #     dp_mat, nc_mat = rand_dp_nc_matrix(512, 512, seed=0)
    #     is_correct, speed_up = compute([dp_mat, 250], [nc_mat, 250], "pow")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

    # def test_0064x0500_pow(self):
    #     # TODO: YOUR CODE HERE
    #     dp_mat, nc_mat = rand_dp_nc_matrix(64, 64, seed=0)
    #     is_correct, speed_up = compute([dp_mat, 500], [nc_mat, 500], "pow")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

    def test_0256x0500_pow(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(256, 256, seed=0)
        is_correct, speed_up = compute([dp_mat, 500], [nc_mat, 500], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    # def test_0512x0500_pow(self):
    #     # TODO: YOUR CODE HERE
    #     dp_mat, nc_mat = rand_dp_nc_matrix(512, 512, seed=0)
    #     is_correct, speed_up = compute([dp_mat, 500], [nc_mat, 500], "pow")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

    # def test_1024x0100_pow(self):
    #     # TODO: YOUR CODE HERE
    #     dp_mat, nc_mat = rand_dp_nc_matrix(1024, 1024, seed=0)
    #     is_correct, speed_up = compute([dp_mat, 100], [nc_mat, 100], "pow")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)