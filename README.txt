Restriction Site Remover

This program reads a DNA sequence that is fed in as a fasta or genbank file and will remove unwanted restriction sites. 
This program can be used to produce BioBrick compatible sequences, but is also helpful with other DNA assembly methods that require restriction sites.
Protein coding sequences are taken into account when using this program. Inside coding sequences the program will make sure that the bases that are changed will not change the amino acid sequence.
Careful! Outside of protein coding sequences this program will not care which bases are changed. This means that promotors or ribosome binding sites etc. will possibly have function removed.
Due to limitations in codon degeneracy, Tryptophan and methionine residues will not be changed by this program so if the restriction site ATGTGG for example, this program will not be able to change that
since it would be changing the protein coding sequence.
The program also has functionality to make your sequence completely BioBrick compatible if you wish by adding a BioBrick prefix and suffix to your sequence

How to use.

Before starting make sure that Biopython and random python modules are installed.
When prompted, follow the instructions to enter the name of the file and which file type you are using into the console.
There are two test files in the repository, one fasta and one genbank file. These can be used to test if the script is working correctly.
The test fasta file is test.seq and the genbank file is test_sequence.gb.
Run build.py with python3. The program has been tested in the terminal and using pycharm professional edition.
The program will ask you which restriction sites you wish to use.
There are 4 restriction sites built in, EcoR1, Xba1, Spe1, Pst1.
If you wish to look for another restriction site please type in what the recognition site is.
The program will check to make sure the entered sequence is valid or not and will prompt you again if it is not.
Your sequence will be checked for restriction recognition sites and if they are found they will be changed.
If you wish to add a BioBrick compatible prefix and suffix it will add those to the start and end of your sequence.
You will be asked to enter the details of your new SeqRecord object (id, name, description).
A prompt will come up to enter the filename for your new SeqRecord object.
The SeqRecord will then be saved as a fasta file.
For a detailed run through of the program check the image test_image.png found in the repository.

Unittesting.

To test the functions in the main build file run coursework_test.py.
This will run through the functions written in the build.py file to test if they work.
It is recommended that you run this file before running the build file to make sure all the functions work and all the correct modules are installed

