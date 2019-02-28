# imports
from Bio import SeqIO
from Bio.Seq import Seq


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
# this is an empty list that will contain restriction sites that need to be changed.
restriction_sites = []
# Here the user is presented with options that will determine which restriction sites are placed into the list
flag = True
first_choice = input("What restriction sites would you like to remove from your DNA sequence? \n"
                     "You can choose from 1) EcoR1, 2) Xba1, 3) Spe1 and 4) Pst1.\n"
                     "If you would like to choose one of these enzymes enter the corresponding "
                     "number or first initial\n"
                     "If you would prefer a different restriction enzyme please enter the recognition sequence.\n")
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
        print("Sorry your input did not work please try again.")
        pass

print(restriction_sites)
#print(remove_restriction_sites("AAGGCCTTAAGAATTCAA"))
