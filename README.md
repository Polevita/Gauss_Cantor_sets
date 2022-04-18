# Gauss_Cantor_sets
Computing Hausdorff dimension of Gauss-Cantor sets with forbidden words of length < 12
Requires the Arb library https://arblib.org/ 
The file setup12.data specifies the set and has the following format 
letters of the alphabet written in first line separated by space 
the precision of the computation required (number of digits)
list of forbidden words; one word per line, letters in the word separated by space 

bash script setup-build.sh takes the setup12.data file and converts it into the file which is necessary for subsequent combinatorics setup
- number of letters in the alphabet 
- alphabet
- the word length for the dynamical system
- number of forbidden words
- list of forbidden words each word is preceded by the number of letters it contains
- the accuracy (number of digits) desired

run: setup-build.h < setup12.data > setup12.day

the file dimsetup.h takes the number of the setup??.day file and performs combinatorics setup. The output is written in several files which contain information about the system: words12.dat, multmatrix12.dat, comultmatrix12.dat,  smallmatrices12.dat, minimatrix12.dat, and init.dat                             
run: ./dimsetup.out 12

Finally, the file computedim.c contains the program in C which computes the Hausdorff dimension of the Markov IFS specified by the files produces by dimsetup. It takes the files ID number (12 in the example above) from the file init.dat. 

run: ./computedim.out
