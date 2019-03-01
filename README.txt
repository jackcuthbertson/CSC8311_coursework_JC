Restriction Site Remover

This program reads a DNA sequence that is fed in as a fasta file and will remove unwanted restriction sites. 
This program can be used to produce BioBrick compatible sequences, but is also helpful with other DNA assembly methods that require restriction sites.
Protein coding sequences are taken into account when using this program. Inside coding sequences the program will make sure that the bases that are changed will not change the amino acid sequence.
Careful! Outside of protein coding sequences this program will not care which bases are changed. This means that promotors or ribosome binding sites etc. will possibly have function removed.

How to use.

Enter the sequence you wish to remove the restriction sites from as a fasta file when prompted.
Run build.py with python3. The program has been tested in the terminal and using pycharm professional edition.
The program will ask you which restriction sites you wish to use.
It has 4 restriction sites built in automatically, EcoR1, Xba1, Spe1, Pst1.
If you wish to look for another restriction site please type in what the recognition site is.
The program will check to make sure the entered sequence is valid or not and will prompt you again if it is not.
