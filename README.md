[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Exact Matrix Factorization Updates for Nonlinear Programming

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data associated with the paper (%%%) by Adolfo R. Escobedo.

The snapshot is based on 
[this SHA](https://github.com/adolfoescobedo/2021.0331) in the development repository.  


**Important: This code is being periodically updated to add more features at 
https://github.com/adolfoescobedo/REF-ROU. Please go there if you would like to
get a more recent version or would like support**

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/%%%

https://doi.org/%%%

Below is the BibTex for citing this snapshot of the repository.

```
@article{escobedo2023exact,
  author=	{Escobedo, Adolfo R},
  publisher=	{INFORMS Journal on Computing},
  title=	{Exact Matrix Factorization Updates for Nonlinear Programming},
  year=		{2023},  
  doi=		{https://doi.org/%%%},
  url=		{https://github.com/adolfoescobedo/2021.0331},
}
```

# REF-ROU
# The REF-ROU repository consists of integer-preserving algorithms for efficiently updating a roundoff-error-free LU or Cholesky factorization. The main executable (ROU.exe) can be used to replicate the three experiments in the paper and for additional computational testing.

## Dependencies

- gcc-g++ (11.0.0 or compatible version)
- make
- gmp (required for working in arbitrary precision arithmetic)


# Installation

Before installing the REF-ROU algorithms, make sure to have the following dependencies installed:

1. gcc-g++:
   - Install compatible version of gcc-g++ compiler.

2. make:
   - Install make using your package manager or from the official website: https://www.gnu.org/software/make/#download

3. gmp:
   - The REF-ROU algorithms depend on the gnu gmp library for working in arbitrary precision arithmetic
   - Install gmp by following the instructions provided in the GNU Multiple Precision Arithmetic Library documentation: [gmp download and installation](https://gmplib.org/#DOWNLOAD)


## Getting Started

The following instructions will help set up the REF-ROU algorithms on your system.


### Installation

1. Clone the repository 
```
git clone https://github.com/adolfoescobedo/2021.0331.git

Install dependencies, as needed

```

2. Configure and build the project:

```
make 

```


## Usage

The main executable associated with the installation of REF-ROU is ROU.exe. There are two possible usages of the ROU update algorithms, described as follows. 

#General ROU update. 
This usage choice of the algorithm (a) builds the REF LU factorization of a dense matrix A, (b) performs a rank-one update of REF-LU(A), and (3) builds the REF LU factorization of the updated matrix, Ah, from scratch (for comparison).

To implement this usage, type <./src/ROU> with the command line argument <-a 0>. The list of optional command line arguments (and the relevant default values) defined for this usage are:

	-c,--check_sol  [0,1]           >default=0<
        -f,--file_out   [string]        >runs_ROU.out<
        -i,--upd_idx    [int] [int]     >default=-1,-1< (relevant only for Up_v4)
                        (if upd_idx_1 = -1, the update vector v is not forced to be linearly dependent)
                        (if upd_idx_1 = 0, v is forced to be linearly dependent using random indices)
                        (if upd_idx_1 > 0, v is forced to be linearly dependent using first up_idx1 columns and first up_idx_2 rows of original matrix)

        -n,--num_vars   [int]           >default=5<
        -p,--print      [0,1]           >default=0<
        -r,--range      [int] [int]     >default=-9,9<
        -s,--seeds      [int] [int]     >default=11,1400<

Note that this usage is relevant to Experiments 1 and 2 of the paper, which can be replicated by running ./scripts/Exp1.sh and ./scripts/Exp2.sh, respectively, using a Linux-based batch processing system.

#Special column replacement update. 
This usage choice of the algorithm (a) builds the REF LU factorization of a dense matrix A, (b) performs a column replacement update of REF-LU(A) using the REF rank-one update algorithm, and (3) performs a column replacement update of REF-LU(A) using the alternative push-and-swap approach of Escobedo and Moreno-Centeno (2017) (for comparison).

To implement this usage, type <./src/ROU> with the command line argument <-a 1>. The list of optional command line arguments (and the relevant default values) defined for this usage are:

	-c,--check_sol  [0,1]           >default=0<
        -f,--file_out   [string]        >runs_CR_ROU.out.out<
        -n,--num_vars   [int]           >default=5<
        -p,--print      [0,1]           >default=0<
        -r,--range      [int] [int]     >default=-9,9<
        -s,--seeds      [int] [int]     >default=11,1400<

Note that this usage is relevant to Experiment 3 of the paper, which can be replicated by running ./scripts/Exp3.sh using a Linux-based batch processing system.


## Ongoing Development

This code is being periodically updated to add more features at the author's  
[Github site](https://github.com/adolfoescobedo/REF-ROU). Please email the author if you have related questions**
