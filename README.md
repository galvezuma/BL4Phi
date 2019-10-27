# BLVector
BLVector is a BLAST Like program originally designed for Xeon-Phi, but that can be executed any current Intel x86 microprocessor with the basic [AVX-512](https://en.wikipedia.org/wiki/AVX-512#CPUs_with_AVX-512) instruction set. In addition to Xeon-Phi, this includes those manufacturated with architecture Skylake-SP, Skylake-X, Cannon Lake and Cascade Lake.
BLVector produces a list in a similar way than BLAST. However, instead of giving a score based on a criteria of its own, BLVector provides a score based on the result of a Smith-Waterman algorithm with affine gaps which, in turn, is the result of applying a score matrix like BLOSUM62 or PAM30. Although BLPhi could be used to align DNA sequences, its potential is obtained when it is applied to peptides. In fact, it uses the Farrar's implementation.
The optional parameters allowed by BLPhi are the next:

*   `-n nearby`: This is an integer that specifies hoy many times a 4-mer must be found in a window of 16 amino acids in the subject to be considered as a potential hit. This integer must be between 1 and 15. Default value: `4`.

*   `-t threshold`: This is an integer that represents the threshold score to be used in the Smith-Waterman algorithm. Those alignment whose final score exceeds this threshold will be reported as hits. Default value: `80`.

*   `-m matrix`: This is the name of the score matrix to be used by the Smith-Waterman algorithm. Default value: `BLOSUM62`.

*   `-g open extend`: These are two integers that represent the costs of opening a new gap or extending an already existing during the pairwise alignment (Smith-Waterman with affine gaps). Default values: `10` and `1`.

*   `-c cluster`: This integer must be between 15 and 21. To improve sensitivity, those amino acids with lowest probability to appear in a peptide can be collapsed into the same letter. Original matrices have 21 different letter and this parameters indicates the final length of the matrix to be applied. Default value: `21` (no cluster).

*   `-p threads`: This integer indicates how many threads must be created in the current BLVector execution. Default value: `228`.

*   `-f`: Fast (non-exhaustive) modifier. When this modifier is included BLVector does not quadruplicate subject sequences but only deduplicates them. This improves execution time but introduces a random component that makes results less accurate, i.e., important hits could be missed. Usually, this is used along to -b modifier. Default value: not set.

*   `-b`: Best (hit only). This retrieves only the best hit after execution. This does not reduce time exceution but only filters the final results, if any, to that with the best score. Default value: not set.

The only mandatory parameters of BLPhi are:

*   `query_filename`: Name of the file with the query Fasta sequence to search for. It may be a multifasta file.

*   `subject_filename`: Name of the Fasta file name with the sequences to search in. It is usually a multifasta file.
