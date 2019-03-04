from build import correct_letters
from build import seq_to_string
from build import replace_codon
from build import remove_restriction_sites
from build import add_chosen_restriction_sites
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import unittest


class TestCoursework(unittest.TestCase):
    def setUp(self):
        self.seq = "AGCT"

    def test_correct_letters(self):
        self.assertTrue(correct_letters(self.seq))
        self.assertFalse(correct_letters("THFB"))

    def test_seq_to_string(self):
        self.assertTrue(isinstance(seq_to_string(SeqIO.index("test.seq", "fasta")), str))

    def test_replace_codon(self):
        self.assertEqual(replace_codon(0, "TTC"), "TTT")
        self.assertEqual(replace_codon(1, "ATTC"), "ATTT")

    def test_add_chosen_restriction_sites(self):
        restriction_sites = []
        add_chosen_restriction_sites("1")
        self.assertEqual(restriction_sites, ["GAATTC"])

    def test_remove_restriction_sites(self):
        add_chosen_restriction_sites("1")
        self.assertNotEqual(remove_restriction_sites("GAATTC"), "GAATTC")


if __name__ == '__main__':
    unittest.main()
