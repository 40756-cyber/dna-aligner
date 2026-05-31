# DNA Sequence Aligner

A Python implementation of the Needleman-Wunsch algorithm
for pairwise DNA sequence alignment.

## What it does
- Reads two DNA sequences from FASTA files
- Builds a scoring matrix using match, mismatch and gap penalties
- Traces back through the matrix to find the optimal alignment
- Displays the alignment with match/mismatch markers

## Sequences compared
- Human insulin gene (NM_000207)
- Mouse insulin gene (NM_008386)

## Key finding
The coding regions of human and mouse insulin are nearly identical,
showing how evolution preserves critical genes across species.
The untranslated end regions show far more variation.

## Tools used
- Python 3

## How to run
1. Download two gene sequences as FASTA files from NCBI
2. Save them as `gene1.fasta` and `gene2.fasta`
3. Run `python dna_aligner.py`

## Algorithm
Needleman-Wunsch (1970) — global pairwise sequence alignment.
Scoring: match +2, mismatch -1, gap -1
