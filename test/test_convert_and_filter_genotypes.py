import unittest
import ingenos
import numpy
import allel

class FilterGenotypesTestCase(unittest.TestCase):

    def setUp(self):
        self.genotypes =\
        allel.GenotypeArray(numpy.random.randint(0,2,(500,25,2)), dtype='i1')
        self.misshapen_sites_boolean = numpy.random.choice([True,False], 10)
        self.sites_boolean =  numpy.random.choice([True,False], 500)
        self.misshapen_samples_boolean = numpy.random.choice([True,False], 75)
    
    def test_input_sites_bool_does_not_match_input_genotype_array(self):
        
        self.assertRaises(ValueError, ingenos.filter_and_convert_genotypes,
                         self.genotypes, self.misshapen_sites_boolean)
                
    def test_input_samples_bool_does_not_match_input_genotype_array(self):
        
        self.assertRaises(ValueError, ingenos.filter_and_convert_genotypes,
                         self.genotypes, self.sites_boolean, 
                         self.misshapen_samples_boolean)


if __name__ == '__main__':
        unittest.main()
