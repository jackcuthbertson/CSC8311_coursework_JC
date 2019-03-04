from build import correct_letters
from build import seq_to_string
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import unittest


class TestCoursework(unittest.TestCase):
    def setUp(self):
        self.seq = "AGCT"

    def test_correct_letters(self):
        self.assertTrue(correct_letters(self.seq))
        self.assertFalse(correct_letters("THFB"))

    def test_seq_to_string(self):
        self.assertTrue(isinstance(SeqRecord(Seq(self.seq)), SeqRecord))
        self.assertTrue(isinstance(seq_to_string(SeqRecord(Seq(self.seq))), str))


if __name__ == '__main__':
    unittest.main()
