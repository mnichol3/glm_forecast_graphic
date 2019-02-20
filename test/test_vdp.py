import unittest
from src.vortex_data_parse import (hello, get_os, hours_mins_2_mins,
                                    dec_2_deg, calc_min_list)

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


        # Test calc_min_list
    def test_6(self):
        # Basic case
        self.min_list = calc_min_list('09121200', '09121205')
        self.assertEqual(self.min_list, ['09121200', '09121201', '09121202',
                                        '09121203', '09121204', '09121205'])

        # Hour change case
        self.min_list = calc_min_list('09121158', '09121201')
        self.assertEqual(self.min_list, ['09121158', '09121159', '09121200',
                                        '09121201'])

        # Day change case
        self.min_list = calc_min_list('09122358', '09130002')
        self.assertEqual(self.min_list, ['09122358', '09122359', '09130000',
                                        '09130001', '09130002'])

        # Month change case (Sep -> Oct)
        self.min_list = calc_min_list('09302358', '10010002')
        self.assertEqual(self.min_list, ['09302358', '09302359', '10010000',
                                        '10010001', '10010002'])

        # Another month change case (Aug -> Sep)
        self.min_list = calc_min_list('08312358', '09010002')
        self.assertEqual(self.min_list, ['08312358', '08312359', '09010000',
                                        '09010001', '09010002'])




if __name__ == '__main__':
    unittest.main()
