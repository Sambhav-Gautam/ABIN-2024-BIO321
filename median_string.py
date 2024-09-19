'''
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
The associated textbook is Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.
The course is run on Coursera and the assignments and textbook are hosted on Stepic.

Problem Title: Median String Problem
Assignment #: 03
Problem ID: B
URL: https://beta.stepic.org/Bioinformatics-Algorithms-2/From-Motif-Finding-to-Finding-a-Median-String-158/#step-7
'''

from itertools import product

def DNA_to_RNA(dna):
    '''Translates DNA to RNA'''
    return dna.replace('T', 'U')

def RNA_to_DNA(rna):
    '''Translates RNA to DNA'''
    return rna.replace('U', 'T')

def ReverseComplementDNA(nucleic_acid):
    '''Returns the reverse complement of a given DNA strand.'''
    nucleotide = 'ATCG'
    complement = 'TAGC'
    transtab = str.maketrans(nucleotide, complement)  # Updated for Python 3
    
    return nucleic_acid.translate(transtab)[::-1].lstrip()

def ReverseComplementRNA(nucleic_acid):
    '''Returns the reverse complement of a given RNA strand.'''
    nucleotide = 'AUCG'
    complement = 'UAGC'
    transtab = str.maketrans(nucleotide, complement)  # Updated for Python 3
    
    return nucleic_acid.translate(transtab)[::-1].lstrip()

def HammingDistance(seq1, seq2):
    'Return the Hamming distance between equal-length sequences.'
    if len(seq1) != len(seq2):
        raise ValueError('Undefined for sequences of unequal length.')
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))

def motif_score(pattern, motif):
    '''Returns the score of d(pattern, motif).'''
    return min([HammingDistance(motif[i:i + len(pattern)], pattern) for i in range(len(motif) - len(pattern) + 1)])

with open('stepic_3b.txt') as input_data:
    k = int(input_data.readline())
    dna_list = [line.strip() for line in input_data.readlines()]

# Initialize the best pattern score as one greater than the maximum possible score.
best_score = k * len(dna_list) + 1
best_patterns = []

# Check the scores of all k-mers.
for pattern in product('ACGT', repeat=k):
    current_score = sum([motif_score(''.join(pattern), dna) for dna in dna_list])
    if current_score < best_score:
        best_score = current_score
        best_patterns = [''.join(pattern)]  # Start a new list with the current best pattern
    elif current_score == best_score:
        best_patterns.append(''.join(pattern))  # Add the pattern to the list

# Print and save the answer.
best_patterns.sort()
print("Best Patterns:", best_patterns)
print("Best Score:", best_score)

# Write to output
with open('Assignment_03B.txt', 'w') as output_data:
    for pattern in best_patterns:
        output_data.write(pattern + '\n')
