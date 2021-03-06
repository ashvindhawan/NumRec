from utils import *
import unittest
from unittest import TestCase

unittest.TestLoader.sortTestMethodsUsing = None

N = 2

class Testsub(TestCase):
    def test_00064_sub(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(64, 64, seed=i)

        is_correct, speed_up = compute(dp_lst, nc_lst, "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_00256_sub(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(256, 256, seed=i)

        is_correct, speed_up = compute(dp_lst, nc_lst, "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_00512_sub(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(512, 512, seed=i)

        is_correct, speed_up = compute(dp_lst, nc_lst, "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_01024_sub(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(1024, 1024, seed=i)
        is_correct, speed_up = compute(dp_lst, nc_lst, "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_02048_sub(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(2048, 2048, seed=i)
        is_correct, speed_up = compute(dp_lst, nc_lst, "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_04096_sub(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(4096, 4096, seed=i)
        is_correct, speed_up = compute(dp_lst, nc_lst, "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_08192_sub(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(8192, 8192, seed=i)
        is_correct, speed_up = compute(dp_lst, nc_lst, "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_16384_sub(self):
        # TODO: YOUR CODE HERE
        dp_lst, nc_lst = [[] for _ in range(N)], [[] for _ in range(N)]
        for i in range(N):
            dp_lst[i], nc_lst[i] = rand_dp_nc_matrix(16384, 16384, seed=i)
        is_correct, speed_up = compute(dp_lst, nc_lst, "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        