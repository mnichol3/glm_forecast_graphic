import unittest
glm_fc_graphic = __import__("glm-forecast-graphic")
from glm_fc_graphic import vortex_data_parse

class LearningCase(unittest.TestCase):
    def test_starting_out(self):
        self.assertEqual(1, 1)


def main():
    unittest.main()

if __name__ == "__main__":
    main()
