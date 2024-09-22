import itertools
from collections import Counter

def Score(S, DNA):
    motifs = [DNA[i][S[i]-1:S[i]-1+l] for i in range(len(S))]
    
    profile = {base: [0] * l for base in 'ACGT'}
    
    for motif in motifs:
        for position in range(l):
            profile[motif[position]][position] += 1

    score = 0
    for position in range(l):
        max_count = max(profile[base][position] for base in 'ACGT')
        score += max_count

    # Construct the consensus string
    consensus = ''
    for position in range(l):
        max_base = max('ACGT', key=lambda base: profile[base][position])
        consensus += max_base

    return score, consensus

def BruteForceMotifSearch(DNA, t, n, l):
    best_scores = {}
    
    for s in itertools.product(range(1, n - l + 2), repeat=t):
        currentScore, consensus = Score(s, DNA)
        if currentScore not in best_scores or best_scores[currentScore][0] < s:
            best_scores[currentScore] = (s, consensus)

    return best_scores

# Example usage
if __name__ == "__main__":
    # Sample DNA sequences
    DNA = [
        "gcggaagagggcactagcccatgtgagagggcaaggacca",
        "atctttctcttaaaaataacataattcagggccaggatgt",
        "gtcacgagctttatcctacagatgatgaatgcaaatcagc",
        "taaaagataatatcgaccctagcgtggcgggcaaggtgct",
        "gtagattcgggtaccgttcataaaagtacgggaatttcgg"
    ]
    for i in range(len(DNA)):
        DNA[i] = DNA[i].upper()
    t = len(DNA)  # number of sequences
    n = len(DNA[0])  # length of each sequence
    l = 6  # length of the motif

    best_scores = BruteForceMotifSearch(DNA, t, n, l)
    
    for score, (positions, consensus) in sorted(best_scores.items(), reverse=True):
        print(f"Score: {score}, Starting Positions: {positions}, Consensus: {consensus}")
