import unittest
from src.vortex_data_parse import hello, get_os, hours_mins_2_mins, dec_2_deg

class TestVortexDataParse(unittest.TestCase):



    def test_1(self):
        self.assertTrue(True)



    def test_2(self):
        self.hello = hello()
        self.assertEqual(self.hello, 'Hello from vortex_data_parse!')



        # Test get_os()
    def test_3(self):
        self.get_os = get_os()
        self.assertEqual(self.get_os, 'linux')



        # Test hours_mins_2_mins()
    def test_4(self):
        self.hrs_2_mins = hours_mins_2_mins('0204')
        self.assertEqual(self.hrs_2_mins, 124)
        self.hrs_2_mins = hours_mins_2_mins('1030')
        self.assertEqual(self.hrs_2_mins, 630)



        # Test dec_2_deg()
    def test_5(self):
        self.dec_2_deg = dec_2_deg('25.12')
        self.assertEqual(self.dec_2_deg, '25.2')
        self.dec_2_deg = dec_2_deg('094.41')
        self.assertEqual(self.dec_2_deg, '-94.683')
        self.dec_2_deg = dec_2_deg('18.23')
        self.assertEqual(self.dec_2_deg, '18.383')
        self.dec_2_deg = dec_2_deg('064.24')
        self.assertEqual(self.dec_2_deg, '-64.4')
        self.dec_2_deg = dec_2_deg('31.48')
        self.assertEqual(self.dec_2_deg, '31.8')
        self.dec_2_deg = dec_2_deg('071.39')
        self.assertEqual(self.dec_2_deg, '-71.65')



if __name__ == '__main__':
    unittest.main()
