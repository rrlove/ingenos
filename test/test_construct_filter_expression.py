import unittest
import ingenos

class ConstructFilterExpressionTestCase(unittest.TestCase):
    
    def setUp(self):
        self.name="2Rq"
        self.mock_inversion_dict = {}
        self.mock_inversion_dict[self.name] = ingenos.Inversion(b'2R', 15,20,750,775)
        
    def test_complex_case(self):
         self.assertEqual(ingenos.construct_filter_expression(self.name,self.mock_inversion_dict,buffer=0),
                         '( ( POS > 15 & POS < 20 ) | ( POS > 750 & POS < 775 ) )') 
    
    def test_simple_case(self):
        self.assertEqual(ingenos.construct_filter_expression(self.name,self.mock_inversion_dict,buffer=0,whole_inversion=True),
                        '( POS > 15 & POS < 775 )')
        
    def test_negative_coordinates_fail(self):
        self.assertRaises(ValueError, ingenos.construct_filter_expression, self.name, self.mock_inversion_dict)

if __name__ == '__main__':
	unittest.main()
