from itertools import product

def score(motifs, DNA):
    """
    Calculate the score of a motif.
    A higher score indicates that the motifs are more conserved.
    """
    l = len(motifs[0])  # Length of each motif
    t = len(motifs)     # Number of sequences
    score = 0
    
    # For each position in the motifs
    for i in range(l):
        counts = {}
        # Count each nucleotide at the ith position of the motif
        for j in range(t):
            nucleotide = motifs[j][i]
            if nucleotide in counts:
                counts[nucleotide] += 1
            else:
                counts[nucleotide] = 1
        # Find the most frequent nucleotide
        max_count = max(counts.values())
        # Add the number of sequences where the most frequent nucleotide is found
        score += max_count
    return score

def brute_force_motif_search(DNA, t, n, l):
    """
    Brute force search for the best motif in the DNA sequences.
    DNA: list of t DNA sequences
    t: number of DNA sequences
    n: length of each DNA sequence
    l: length of the motif to search for
    """
    bestScore = 0
    bestMotif = []
    
    # Generate all possible starting positions for motifs in each sequence
    for starts in product(range(n-l+1), repeat=t):
        # Extract motifs from each sequence starting at the positions in starts
        motifs = [DNA[i][starts[i]:starts[i]+l] for i in range(t)]
        
        # Calculate the score of the current motif set
        currentScore = score(motifs, DNA)
        
        # Update bestScore and bestMotif if current motif is better
        if currentScore > bestScore:
            bestScore = currentScore
            bestMotif = motifs
            
    return bestMotif

# Example usage
DNA = [
    "AGCTAGCTAG",
    "CGTAGCTAGC",
    "TAGCTAGCTA",
    "GCTAGCTAGC"
]

t = len(DNA)  # Number of sequences
n = len(DNA[0])  # Length of each sequence
l = 3  # Length of the motif

best_motif = brute_force_motif_search(DNA, t, n, l)
print("Best Motif:", best_motif)