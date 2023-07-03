[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Exact Matrix Factorization Updates for Nonlinear Programming

The folder contains the source co\de for the ROU executable; its functionality and usage are defined in the README.md file within the main directory. The folder also contains the header files that define the algorithms associated with the ROU executable. There are a total of 5 header files; their brief descriptions are below:

1. GF.h: Defines various basic functions (copy, multiply, print) specific to mpz and mpq vector and matrix data structures.

2. Exact_LU_Alg.h: Defines the exact matrix factorization objects from which the REF LU factorization and the exactional rational arithmetic factorizations are built; it also includes basic initialization and output validation subroutines.

3. IPGE.h: Defines integer-preserving algorithms that form the core of the REF factorization framework. The main functionality of the defined subroutines includes (a) initializing and constructing the REF LU factorization of a square matrix, (b) using the REF LU factorization to solve systems of linear equations via REF forward and backward substitution, (c) performing a REF LU rank-one update, (d) performing a column replacement update via the ROU algorithms, and (e) performing a column replacement update via the push-and-swap approach defined in Escobedo and Moreno-Centeno (2017).

4. QGE.h: Defines rational arithmetic algorithms that form the core of the exact rational arithmetic LU factorization subroutines tested in the paper. The main functionality of the defined subroutines includes (a) initializing and constructing the rational LU factorization of a square matrix via the Doolittle and Crout LU algorithms, (b) using the rational factorization to solve systems of linear equations via rational forward and backward substitution, (c) performing a column replacement update via the Bartels and Golub (1969) algorithm.

5. cmdOpt.h: Defines the object for handling command line inputs associated with the ROU executable. 