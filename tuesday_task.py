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
    found_proteins = []
    count = 0
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
            if not found and "gene_name_synonym" in entry.annotations.keys():
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
