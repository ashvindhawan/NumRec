from utils import *
import numc as nc 
import dumbpy as dp
import unittest 

dp_mat2, nc_mat2 = rand_dp_nc_matrix(4, 2, seed=12)
  
class SetTestCase2(unittest.TestCase): 

    def test_1(self): 
        with self.assertRaises(TypeError): 
            dp_mat2.set(2147483648, 0, 1)
        with self.assertRaises(TypeError): 
            nc_mat2.set(2147483648, 0, 1)
        with self.assertRaises(TypeError): 
            nc_mat2.set(-2147483649, 0, 1)
        with self.assertRaises(TypeError): 
            nc_mat2.set(-2147483649, 0, 1)
        return

    def test_2(self):
        with self.assertRaises(IndexError): 
            dp_mat2.set(2147483647, 0, 1)
        with self.assertRaises(IndexError): 
            dp_mat2.set(2147483647, 0, 1)
        with self.assertRaises(IndexError): 
            dp_mat2.set(-2147483648, 0, 1)
        with self.assertRaises(IndexError): 
            dp_mat2.set(-2147483648, 0, 1)
    
if __name__ == '__main__':  
    unittest.main()