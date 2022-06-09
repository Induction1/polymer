import unittest

from tools import down_sample


class ProcessFileTestCase(unittest.TestCase):
    def test_process_file(self):
        down_sample(in_file_name='data14.csv', out_file_name='out14.csv')




