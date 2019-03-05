from build import correct_letters
from build import seq_to_string
from build import replace_codon
from build import remove_restriction_sites
from build import add_chosen_restriction_sites
from build import create_seqrecord
from build import restriction_sites
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import unittest


class TestCoursework(unittest.TestCase):
    def setUp(self):
        self.seq = "AGCT"
        self.file = SeqIO.index("test.seq", "fasta")

    def test_correct_letters(self):
        self.assertTrue(correct_letters(self.seq))
        self.assertFalse(correct_letters("THFB"))

    def test_seq_to_string(self):
        self.assertFalse(isinstance(self.file, str))
        self.assertTrue(isinstance(seq_to_string(self.file), str))
        self.file.close()

    def test_replace_codon(self):
        self.assertEqual(replace_codon(0, "TTC"), "TTT")
        self.assertEqual(replace_codon(1, "ATTC"), "ATTT")
        self.assertNotEqual(replace_codon(0, "GAT"), "GAT")

    def test_add_chosen_restriction_sites(self):
        add_chosen_restriction_sites("1")
        self.assertEqual(restriction_sites, ["GAATTC"])

    def test_remove_restriction_sites(self):
        add_chosen_restriction_sites("1")
        self.assertNotEqual(remove_restriction_sites("GAATTC"), "GAATTC")

    def test_create_seqrecord(self):
        self.assertTrue(isinstance(self.seq, str))
        self.assertFalse(isinstance(create_seqrecord(self.seq, "none", "none", "none"), str))
        self.assertTrue(isinstance(create_seqrecord(self.seq, "none", "none", "none"), SeqRecord))


if __name__ == '__main__':
    unittest.main()
