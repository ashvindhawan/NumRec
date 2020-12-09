from utils import *
import unittest
from unittest import TestCase

class TestNeg(TestCase):
    def test_1024_neg(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1024, 1024, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_2048_neg(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(2048, 2048, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_4096_neg(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(4096, 4096, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)