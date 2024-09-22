import time
import random
import string

def brute_force_motif_search(dna, t, n, l):
    bestScore = 0
    bestMotif = None

    s = [0] * t  # Initialize starting indices to 0 for all sequences
    total_iterations = (n - l + 1) ** t  # Calculate total iterations based on the length of the sequences and the motif length

    for _ in range(total_iterations):
        current_motifs = []
        profile_matrix = {'A': [0]*l, 'C': [0]*l, 'G': [0]*l, 'T': [0]*l}  # Initialize profile matrix with empty lists
        for i in range(t):
            current_sequence = dna[i]
            start = s[i]  # Start index
            current_motifs.append(current_sequence[start:start + l])

        # Profile matrix
        for j in range(l):
            counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            for motif in current_motifs:
                counts[motif[j].upper()] += 1  # Increment the letter found
            for nucleotide in profile_matrix:  # Append the counts found
                profile_matrix[nucleotide][j] = counts[nucleotide]

        # Consensus
        consensus = []
        score = 0
        for i in range(l):
            max_in_column, nucleotide = max((profile_matrix[nucleotide][i], nucleotide) for nucleotide in profile_matrix)
            score += max_in_column
            consensus.append(nucleotide)

        # Construct consensus string with correct case
        consensus_str = ''.join([consensus[i].upper() if current_motifs[0][i].isupper() else consensus[i].lower() for i in range(l)])

        # Calculate best score
        if score > bestScore:
            bestScore = score
            bestMotif = consensus_str

        j = t - 1  # j is the last dna sequence in the list
        while j >= 0:  # A loop that iterates over each index in the s array
            if s[j] < n - l:  # Check that the starting index is less than the maximum possible starting index
                s[j] += 1  # Increment start index
                break
            else:  # Reset start index to zero
                s[j] = 0
                j -= 1

    print("Best Motif:", bestMotif)
    print("Best Score:", bestScore)


def hamming_dist(word, seq):
    minDist = float('inf')

    for i in range(len(seq) - len(word) + 1):  # Iterate over seq
        pattern = seq[i:i + len(word)]  # Get the pattern from the seq
        mismatch = 0

        for j in range(len(word)):
            if word[j] != pattern[j]:  # Compare
                mismatch += 1  # Increment mismatch

        minDist = min(minDist, mismatch)
    return minDist


def MedianStringSearch(sequences, num_sequences, sequence_length, motif_length):
    best_motif = None
    best_distance = float('inf')
    motifs_by_distance = {}  # Dictionary to store motifs by their distance

    for i in range(sequence_length - motif_length + 1):  # Iterate on all possible sequences
        motif = sequences[0][i:i + motif_length]  # Extract motif

        total_distance = sum(hamming_dist(motif, seq) for seq in sequences)  # Calculate distance from current motif

        # Store motifs by their distance
        if total_distance not in motifs_by_distance:
            motifs_by_distance[total_distance] = []
        motifs_by_distance[total_distance].append(motif)

        # Update best motif found
        if total_distance < best_distance:
            best_distance = total_distance
            best_motif = motif

    print("Best Median Motif:", best_motif)
    print("Best Hamming Distance:", best_distance)
    
    # Print all motifs with the same Hamming distance
    for distance, motifs in motifs_by_distance.items():
        motifs.sort()
        print(f"Hamming Distance: {distance}, Motifs: {motifs}")


print("Motif Finding Problem Solver\n")

file_path = 'rawDNA.txt'
file = open(file_path, 'r')
file_lines = file.readlines()

# Remove \n
file_lines = [line.strip() for line in file_lines]

# Reading t, n, l from the first line
sequences, nucleotides, pattern_length = file_lines[0].split()
sequences = int(sequences)
nucleotides = int(nucleotides)
pattern_length = int(pattern_length)

# Remove the first line
DNA = file_lines[1:]

print("\nBrute Force Algorithm: ")
brute_force_motif_search(DNA, sequences, nucleotides, pattern_length)
# print("\nMedian String Algorithm: ")
# MedianStringSearch(DNA, sequences, nucleotides, pattern_length)
