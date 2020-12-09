from utils import *
import unittest
from unittest import TestCase

unittest.TestLoader.sortTestMethodsUsing = None

N = 2

class TestMul(TestCase):
    def test_0064_mul(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(256, 256, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(256, 256, seed=1)

        dp_lst = [dp_mat1, dp_mat2]
        nc_lst = [nc_mat1, nc_mat2]

        is_correct, speed_up = compute(dp_lst, nc_lst, "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_0256_mul(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(256, 256, seed=i)

        is_correct, speed_up = compute(dp_lst, nc_lst, "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_0512_mul(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(512, 512, seed=i)

        is_correct, speed_up = compute(dp_lst, nc_lst, "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_1024_mul(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(1024, 1024, seed=i)
        is_correct, speed_up = compute(dp_lst, nc_lst, "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_2048_mul(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(2048, 2048, seed=i)
        is_correct, speed_up = compute(dp_lst, nc_lst, "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)