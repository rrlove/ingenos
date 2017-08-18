import unittest
import ingenos
import numpy

class PruneByLDTestCase(unittest.TestCase):
    
    def setUp(self):
        
        self.allele_counts = numpy.random.randint(0,2,(1000,25))
        
    def test_window_size_is_positive(self):
        
        self.assertRaises(ValueError, ingenos.prune_by_LD, self.allele_counts, -50, 100, 0.2)
        
    def test_step_size_is_positive(self):
        
        self.assertRaises(ValueError, ingenos.prune_by_LD, self.allele_counts, 50, -5, 0.2)
        
    def test_r2_is_positive(self):
        
        self.assertRaises(ValueError, ingenos.prune_by_LD, self.allele_counts, 50, 5, -0.2)
        
    def test_output_array_is_no_larger_than_input(self):
        
        output = ingenos.prune_by_LD(self.allele_counts, 50, 5, 0.2)
        self.assertLessEqual(len(output), len(self.allele_counts))

if __name__ == '__main__':
	unittest.main()
