# imports
from Bio import SeqIO
from Bio.Seq import Seq


# functions
# this function will add restriction sites to be removed to a list based on user input
def add_chosen_restriction_sites(chosen_site):
    if chosen_site in ["1", "e", "E"]:
        restriction_sites.append("GAATTC")
    elif chosen_site in ["2", "x", "X"]:
        restriction_sites.append("TCTAGA")
    elif chosen_site in ["3", "s", "S"]:
        restriction_sites.append("ACTAGT")
    elif chosen_site is ["4", "p", "P"]:
        restriction_sites.append("CTGCAG")
    else:
        print("Your input does not work please try again")


# When a restriction site is found this function will  be called and pass through the site finding a codon it can change
# while attempting to not change the amino acids expressed in the codon sequence
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


# This function will pass through a DNA sequence until it finds a restriction site. it will then call the replace_codon
# function which will remove it. It then replaces the old site with the new one and continues searching.
def remove_restriction_sites(input_DNA):
    for strand, nuc in [(+1, input_DNA)]:
        for restriction_site in restriction_sites:
            for i in range(len(nuc)):
                r = len(restriction_site)
                if nuc[i:i+r] == restriction_site:
                    nuc = nuc[:i] + nuc[i:i+r].replace(nuc[i:i+r], replace_codon(nuc[i:i+r])) + nuc[i+r:]
    return nuc


# main
restriction_sites = []
flag = True
first_choice = input("What restriction sites would you like to remove from your DNA sequence? \n"
                     "You can choose from 1) EcoR1, 2) Xba1, 3) Spe1 and 4) Pst1.\n"
                     "If you would like to choose one of these enzymes enter the corresponding "
                     "number or first initial\n")
print(first_choice)
add_chosen_restriction_sites(first_choice)
while flag is True:
    user_input = input("Would you like to add another restriction site? (Y/N)\n")
    if user_input in ["Y", "Yes", "y", "yes"]:
        add_chosen_restriction_sites(input("Which restriction site would you like to add?\n"
                                           "1) EcoR1, 2) Xba1, 3) Spe1 and 4) Pst1\n"))
    elif user_input in ["N", "No", "n", "no"]:
        flag = False
    else:
        pass

print(restriction_sites)
#print(remove_restriction_sites("AAGGCCTTAAGAATTCAA"))
