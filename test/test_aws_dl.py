import unittest
from src.aws_dl import get_os, padding_zero, date_time_chunk, calc_julian_day

class TestAWS_DL(unittest.TestCase):


    def test_1(self):
        self.assertTrue(True)



        # Test get_os()
    def test_2(self):
        self.get_os = get_os()
        self.assertEqual(self.get_os, 'linux')



        # Test padding_zero()
    def test_3(self):
        self.padding_zero = padding_zero(10, 10)
        self.assertEqual(self.padding_zero, '10')

        self.padding_zero = padding_zero(9, 10)
        self.assertEqual(self.padding_zero, '09')

        self.padding_zero = padding_zero(1, 10)
        self.assertEqual(self.padding_zero, '01')

        self.padding_zero = padding_zero(10, 100)
        self.assertEqual(self.padding_zero, '010')

        self.padding_zero = padding_zero(1, 100)
        self.assertEqual(self.padding_zero, '001')

        self.padding_zero = padding_zero(00, 10)
        self.assertEqual(self.padding_zero, '00')

        self.padding_zero = padding_zero(20, 10)
        self.assertEqual(self.padding_zero, '20')

        self.padding_zero = padding_zero(1, 1000)
        self.assertEqual(self.padding_zero, '0001')

        self.padding_zero = padding_zero('9', 10)
        self.assertEqual(self.padding_zero, '09')

        self.padding_zero = padding_zero('i', 10)
        self.assertEqual(self.padding_zero, '0i')



        # Test date_time_chunk()
    def test_4(self):
        self.chunk = date_time_chunk('2018091218', '2018091305')
        self.assertEqual(self.chunk, ['201809121800', '201809121900',
            '201809122000', '201809122100', '201809122200', '201809122300',
            '201809130000', '201809130100', '201809130200', '201809130300',
            '201809130400', '201809130500'])

        self.chunk = date_time_chunk('2018093023', '2018100102')
        self.assertEqual(self.chunk, ['201809302300', '201810010000',
            '201810010100', '201810010200'])

        self.chunk = date_time_chunk('2018123123', '2019010102')
        self.assertEqual(self.chunk, ['201812312300', '201901010000',
            '201901010100', '201901010200'])



            # Test calc_julian_day()
    def test_5(self):
        self.jd = calc_julian_day('20180912')
        self.assertEqual(self.jd, '255')

        self.jd = calc_julian_day('20180925')
        self.assertEqual(self.jd, '268')

        self.jd = calc_julian_day('20181105')
        self.assertEqual(self.jd, '309')

        self.jd = calc_julian_day('20180101')
        self.assertEqual(self.jd, '001')

        self.jd = calc_julian_day('20181231')
        self.assertEqual(self.jd, '365')

        self.jd = calc_julian_day('20180310')
        self.assertEqual(self.jd, '069')

        # Some fun with leap year cases
        self.jd = calc_julian_day('20160201')
        self.assertEqual(self.jd, '032')

        self.jd = calc_julian_day('20160916')
        self.assertEqual(self.jd, '260')

        self.jd = calc_julian_day('20161231')
        self.assertEqual(self.jd, '366')

        self.jd = calc_julian_day('20160309')
        self.assertEqual(self.jd, '069')




if __name__ == '__main__':
    unittest.main()
