from itertools import product

# Sample dataset (DNA strings) - you can replace these with your input strings
D = [
    "gcggaagagggcactagcccatgtgagagggcaaggacca",
    "atctttctcttaaaaataacataattcagggccaggatgt",
    "gtcacgagctttatcctacagatgatgaatgcaaatcagc",
    "taaaagataatatcgaccctagcgtggcgggcaaggtgct",
    "gtagattcgggtaccgttcataaaagtacgggaatttcgg",
]

for i in range(len(D)):
    D[i] = D[i].upper()
l = 6  # Length of the l-mer strings to search for

# Function to compute Hamming distance between two strings
def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# Function to compute the TotalDistance of a given l-mer v to all sequences in D
def total_distance(v, D):
    total_dist = 0
    for seq in D:
        min_dist = float('inf')
        # Iterate over all possible starting positions in each sequence
        for i in range(len(seq) - len(v) + 1):
            subseq = seq[i:i + len(v)]
            dist = hamming_distance(v, subseq)
            if dist < min_dist:
                min_dist = dist
        total_dist += min_dist
    return total_dist

# Alphabet for DNA sequences
alphabet = "ACGT"

# Generate all possible l-mer strings from the alphabet
all_lmers = [''.join(p) for p in product(alphabet, repeat=l)]
print(len(all_lmers))
# Initialize variables to store the best l-mers and the minimum total distance
best_lmers = []
min_total_distance = float('inf')

# Iterate over each l-mer and compute its TotalDistance
for lmer in all_lmers:
    dist = total_distance(lmer, D)
    if dist < min_total_distance:
        min_total_distance = dist
        best_lmers = [lmer]  # Start a new list with the new best l-mer
    elif dist == min_total_distance:
        best_lmers.append(lmer)  # Add to the list of best l-mers

# Output the result
print(f"Best l-mers: {best_lmers}")
print(f"Minimum Total Distance: {min_total_distance}")
