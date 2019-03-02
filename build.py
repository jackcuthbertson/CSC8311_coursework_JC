# imports
from Bio import SeqIO
import random


# functions
# This function will determine if a sequence only has AGCT characters in it. if it does it returns true if it doesnt it
# returns false.
def correct_letters(seq):
    x = True
    for letter in seq:
        if letter in ["A", "T", "C", "G", "a", "c", "t", "g"]:
            pass
        else:
            print("The sequence entered ({}) does not contain correct characters. Please try again.".format(seq))
            x = False
            break
    return x


# This function asks the user what file contains their desired sequence and what file type. The function is tested with
# fasta and genbank files
def enter_file_name_type():
    filename = input("What is the name of the file you wish to enter? (include file extension)\n")
    filetype = input("Which file format do you wish to use? 1) fasta 2) genbank\n"
                     "Enter the number or the first initial of your desired format.\n")
    if filetype in ["1", "f", "F", "fasta", "Fasta", "FASTA"]:
        filetype = "fasta"
    elif filetype in ["2", "g", "G", "genbank", "Genbank", "GENBANK"]:
        filetype = "genbank"
    try:
        input_sequence = SeqIO.index(filename, filetype)
    except ValueError:
        print("Sorry the file format you entered was incorrect please try again.")
        input_sequence = 0
    except FileNotFoundError:
        print("Sorry the file you were looking for was not found.")
        input_sequence = 0
    if input_sequence == 0:
        return enter_file_name_type()
    else:
        output = seq_to_string(input_sequence)
        return output


# Here an input fasta or genbank file will be converted to a string.
def seq_to_string(input_sequence):
    sequence = 0
    for keys, values in input_sequence.items():
        sequence = values.seq
    str(sequence)
    return sequence


# this function will add restriction sites to be removed to a list based on user input
def add_chosen_restriction_sites(chosen_site):
    if chosen_site in ["1", "e", "E"]:
        restriction_sites.append("GAATTC")
    elif chosen_site in ["2", "x", "X"]:
        restriction_sites.append("TCTAGA")
    elif chosen_site in ["3", "s", "S"]:
        restriction_sites.append("ACTAGT")
    elif chosen_site in ["4", "p", "P"]:
        restriction_sites.append("CTGCAG")
    else:
        y = 1
        if len(chosen_site) > 10 or len(chosen_site) < 2:
            print("Your input restriction site ({}) is either too long or too short.\n"
                  "Your restriction site must be between 3 and 9 bases long\n"
                  "Please try again.".format(chosen_site))
            y = 0
        if correct_letters(chosen_site) is False:
            y = 0
        if y == 1:
            chosen_site = chosen_site.upper()
            restriction_sites.append(chosen_site)


# When a restriction site is found this function will be called and pass through the site finding a codon it can change
# while attempting to not change the amino acids expressed in the codon sequence
def replace_codon(frame, input_restriction_site):
    for i in range(frame, len(input_restriction_site), 3):
        checked_codon = input_restriction_site[i:i+3]
        replaced_codon = input_restriction_site[i:i+3]
        for aa in codon_table.values():
            for codon in aa:
                if checked_codon == codon:
                    aa.remove(checked_codon)
                    replaced_codon = random.choice(aa)
                    aa.append(checked_codon)
                    print(checked_codon, replaced_codon)
        if input_restriction_site[i:i+3] != replaced_codon:
            input_restriction_site = input_restriction_site.replace(checked_codon, replaced_codon)
            return input_restriction_site


# This function will pass through a DNA sequence until it finds a restriction site. it will then call the replace_codon
# function which will remove it. It then replaces the old site with the new one and continues searching.
def remove_restriction_sites(input_DNA):
    restriction_sites_changed = 1
    for strand, nuc in [(+1, input_DNA)]:
        for restriction_site in restriction_sites:
            x = 0
            for i in range(len(nuc)):
                r = len(restriction_site)
                stop_codons = ["TAA", "TAG", "TGA"]
                if "ATG" in nuc[i:i+r] and x == 0:
                    x = 1
                    start_codon = nuc.index("ATG", i, i+r)
                for stop_codon in stop_codons:
                    if stop_codon in nuc[i:i+r] and x == 1:
                        if (nuc.index(stop_codon, i, i+r) - start_codon) % 3 == 0:
                            x = 0
                if nuc[i:i+r] == restriction_site:
                    if x == 1:
                        print(str(restriction_sites_changed) + " : " + str(nuc[i:i+r]) + "->" +
                              nuc[i:i+r].replace(nuc[i:i+r], replace_codon((i-start_codon) % 3, nuc[i:i+r])) +
                              " at position {} to {}".format(i, i+r))
                        nuc = nuc[:i] + nuc[i:i+r].replace(nuc[i:i+r], replace_codon((i-start_codon) % 3, nuc[i:i+r])) \
                            + nuc[i+r:]
                        restriction_sites_changed += 1
                    else:
                        print(str(restriction_sites_changed) + " : " + str(nuc[i:i+r]) + "->" +
                              nuc[i:i+r].replace(nuc[i:i+r], replace_codon(0, nuc[i:i+r])) +
                              " at position {} to {}".format(i, i+r))
                        nuc = nuc[:i] + nuc[i:i+r].replace(nuc[i:i+r], replace_codon(0, nuc[i:i+r])) + nuc[i+r:]
                        restriction_sites_changed += 1
    return nuc


# main
# This dictionary contains all of the codons for amino acids which have more than one codon. These are fed into the
#
codon_table = {
    "Phe": ["TTC", "TTT"],
    "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "Ile": ["ATT", "ATC", "ATA"],
    "Val": ["GTT", "GTC", "GTA", "GTG"],
    "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "Pro": ["CCT", "CCC", "CCA", "CCG"],
    "Thr": ["ACT", "ACC", "ACA", "ACG"],
    "Ala": ["GCT", "GCC", "GCA", "GCG"],
    "Tyr": ["TAT", "TAC"],
    "His": ["CAT", "CAC"],
    "Gln": ["CAA", "CAG"],
    "Asn": ["AAT", "AAC"],
    "Lys": ["AAA", "AAG"],
    "Asp": ["GAT", "GAC"],
    "Glu": ["GAA", "GAG"],
    "Cys": ["TGT", "TGC"],
    "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Gly": ["GGT", "GGC", "GGA", "GGG"]
}
# This dictionary contains codons for amino acids with just one codon. These codons cannot be changed in  protein coding
# sequences so the replace_codon function will have to choose another.
unchangeable_codon_table = {
    "Met": ["ATG"],
    "Trp": ["TGG"]
}


sequence1 = enter_file_name_type()

# this is an empty list that will contain restriction sites that need to be changed.
restriction_sites = []
changed_restriction_sites = {}
# Here the user is presented with options that will determine which restriction sites are placed into the list
flag = True
first_choice = input("What restriction sites would you like to remove from your DNA sequence? \n"
                     "You can choose from 1) EcoR1, 2) Xba1, 3) Spe1 and 4) Pst1.\n"
                     "If you would like to choose one of these enzymes enter the corresponding "
                     "number or first initial\n"
                     "If you would prefer a different restriction enzyme please enter the recognition sequence.\n")

add_chosen_restriction_sites(first_choice)

while flag is True:
    user_input = input("Would you like to add another restriction site? (Y/N)\n")
    if user_input in ["Y", "Yes", "y", "yes"]:
        add_chosen_restriction_sites(input("Which restriction site would you like to add?\n"
                                           "1) EcoR1, 2) Xba1, 3) Spe1 and 4) Pst1\n"))
    elif user_input in ["N", "No", "n", "no"]:
        flag = False
    else:
        print("Sorry your input did not work please try again.")
        pass


print("Your chosen restriction sites are " + str(restriction_sites))

print(remove_restriction_sites(str(sequence1)))
