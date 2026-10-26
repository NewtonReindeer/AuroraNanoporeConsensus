import unittest
from cmv_preconsensus_qc import classify_indel, homopolymer_run

class TestPreconsensusQC(unittest.TestCase):
    def test_frame_classes(self):
        self.assertEqual(classify_indel(1, 5), ('frame_disrupting', 'homopolymer_context'))
        self.assertEqual(classify_indel(3, 1), ('in_frame', 'other_context'))
        self.assertEqual(classify_indel(6, 4), ('in_frame', 'homopolymer_context'))
    def test_homopolymer(self):
        self.assertEqual(homopolymer_run('ACCCCT', 2), 4)
        self.assertEqual(homopolymer_run('ACCCCT', 0), 1)
    def test_invalid(self):
        with self.assertRaises(ValueError):
            classify_indel(0, 0)

if __name__ == '__main__':
    unittest.main()
