# imports
from Bio import SeqIO
import gzip
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo


# functions
def find_proteins():
    desired_protein = input("please input the protein you are looking for? ")
    # User will input what protein they are looking for.
    limit_results = int(input("how many proteins do you want? "))
    # User will input how many proteins they are looking for.
    find_synonym = input("Do you wish to find proteins with secondary relation to your search? (y/n) ")
    # this input command will determine if the function will find only the proteins by their primary annotation or
    # secondary annotation as well
    if find_synonym == "y" or "yes":
        flag = True
    elif find_synonym == "n" or "no":
        flag = False
    found_proteins = []  # creates an empty list to store contained proteins
    count = 0  # variable to count how many entries have been found
    x = False
    with gzip.open("uniprot_sprot.xml.gz") as f:
        for entry in SeqIO.parse(f, "uniprot-xml"):
            if x is True:
                break
            found = False
            if "gene_name_primary" in entry.annotations.keys():
                value = entry.annotations["gene_name_primary"].lower()
                found = desired_protein in value
                print("id:{}, Primary: {}".format(entry.id, value)) if found else ""
            if not found and "gene_name_synonym" in entry.annotations.keys() and flag is True:
                values = [x.lower() for x in entry.annotations["gene_name_synonym"]]
                for value in values:
                    if desired_protein in value:
                        found = True
                        break
                print("Id:{}, Secondary:{}".format(entry.id, str(values))) if found else ""
            if found:
                found_proteins.append(entry)
                count += 1
                if count == limit_results:
                    x = True
    with open("found_proteins.fasta", "w") as f:
        SeqIO.write(found_proteins, f, "fasta")


def alignment():
    cline = ClustalwCommandline("clustalw", infile="found_proteins.fasta")
    print(cline)
    stdout, stderr = cline()

    align = AlignIO.read("found_proteins.aln", "clustal")
    print(align)

    tree = Phylo.read("found_proteins.dnd", "newick")
    Phylo.draw_ascii(tree)


# main
find_proteins()
alignment()
