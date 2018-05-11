# -*- coding: utf-8 -*-
"""
Created on Thu May  3 10:42:52 2018

@author: beccalove
"""

import unittest

import ingenos
import numpy as np

class ComputeConcordanceStratTestCase(unittest.TestCase):
        
    def test_standard(self):
        ##tests standard use
        test_alleles = np.array([0,1,2,0,1,2,0,1,2])
        test_is_called = np.array([False,True,True,True,False,
                                   True,True,True,False])
        test_karyos = np.array([0,0,0,1,1,1,2,2,2])
        
        test_output = (0,0,0,0,(2/3),(2/3),(2/3),(2/3))
        
        self.assertAlmostEqual(ingenos.compute_concordance_strat(test_alleles, 
                         test_is_called, test_karyos), test_output)
        
    def test_no_called(self):
        ##tests samples not being called
        test_alleles = np.array([0, 0, 0, 0])
        
        test_is_called = np.array([False, False, False, False])
        
        test_karyos = np.array([0, 1, 0, 0])
        
        test_output = (np.nan, np.nan, np.nan, np.nan, 0, 0, np.nan, 0)
        
        np.testing.assert_almost_equal(
                ingenos.compute_concordance_strat(test_alleles,
                            test_is_called, test_karyos), test_output)
        
    def test_incorrect_shape(self):
        ##tests objects of incorrect shapes
        
        test_alleles = np.array([0, 1, 2])
        
        test_is_called = np.array([False, False])
        
        test_karyos = np.array([0, 1, 2])
        
        self.assertRaises(AssertionError, ingenos.compute_concordance_strat,
                          test_alleles, test_is_called, test_karyos)
        
        
    def test_karyos_too_big(self):
        ##tests karyotypes > 2
        
        test_alleles = np.array([2, 2, 2])
        
        test_is_called = np.array([True, True, True])
        
        test_karyos = np.array([3, 2, 2])
    
        self.assertRaises(AssertionError, ingenos.compute_concordance_strat,
                          test_alleles, test_is_called, test_karyos)
        
        
    def test_karyos_too_small(self):
        ##test negative karyotypes
        
        test_alleles = np.array([2, 2, 2])
        
        test_is_called = np.array([True, True, True])
        
        test_karyos = np.array([2, -5, 2])
        
        self.assertRaises(AssertionError, ingenos.compute_concordance_strat,
                          test_alleles, test_is_called, test_karyos)
        
    def test_alleles_too_small(self):
        ##test negative alt allele counts
        
        test_alleles = np.array([-1, 2, 0])
        
        test_is_called = np.array([True, True, True])
        
        test_karyos = np.array([1, 2, 0])
        
        self.assertRaises(AssertionError, ingenos.compute_concordance_strat,
                          test_alleles, test_is_called, test_karyos)
        
        
    def test_alleles_too_large(self):
        ##test alt allele counts > 2
        
        test_alleles = np.array([0, 30, 0])
        
        test_is_called = np.array([True, True, True])
        
        test_karyos = np.array([1, 2, 0])
        
        self.assertRaises(AssertionError, ingenos.compute_concordance_strat,
                          test_alleles, test_is_called, test_karyos)        

if __name__ == '__main__':
        unittest.main()
        