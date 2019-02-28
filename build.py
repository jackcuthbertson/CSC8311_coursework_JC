# imports
from Bio import SeqIO
from Bio.Seq import Seq


# functions
def replace_codon(input_restriction_site):
    for i in range(0, len(input_restriction_site)):
        checked_codon = input_restriction_site[i:i+3]
        replaced_codon = input_restriction_site[i:i+3]
        if checked_codon == "TTT":
            replaced_codon = "TTC"
        if checked_codon == "TTC":
            replaced_codon = "TTT"
        if input_restriction_site[i:i+3] != replaced_codon:
            input_restriction_site = input_restriction_site.replace(checked_codon, replaced_codon)
    return input_restriction_site


def remove_restriction_sites(input_DNA):
    for strand, nuc in [(+1, input_DNA)]:
        for restriction_site in restriction_sites:
            for i in range(len(nuc)):
                r = len(restriction_site)
                if nuc[i:i+r] == restriction_site:
                    nuc = nuc[:i] + nuc[i:i+r].replace(nuc[i:i+r], replace_codon(nuc[i:i+r])) + nuc[i+r:]
    return nuc


# main
restriction_sites = ["GAATTC"]

print(remove_restriction_sites("AAGGCCTTAAGAATTCAA"))
