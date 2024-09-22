import random
import matplotlib.pyplot as plt

# Convert symbol to corresponding number
def symbolToNumber(symbol): 
    if symbol == "A":
        return 0
    if symbol == "C":
        return 1
    if symbol == "G":
        return 2
    if symbol == "T":
        return 3

# Convert number back to corresponding symbol
def numberToSymbol(x): 
    if x == 0:
        return "A"
    if x == 1:
        return "C"
    if x == 2:
        return "G"
    if x == 3:
        return "T"

# Randomly select a k-mer based on profile probabilities
def profileRandom(k, profile, text): 
    probs = [] 
    for i in range(0, len(text) - k + 1):  # Iterate over each k-mer in text
        prob = 1.0
        pattern = text[i:i+k]  # Extract k-mer
        for j in range(k):  # For each symbol in the k-mer
            l = symbolToNumber(pattern[j])  # Get corresponding index for the symbol
            prob *= profile[l][j]  # Multiply probabilities from profile
        probs.append(prob)  # Store the probability for the current k-mer
    r = myRandom(probs)  # Select a random k-mer based on probabilities
    return r

# Form a profile from motifs
def profileForm(motifs): 
    k = len(motifs[0])  # k is the length of the k-mer
    t = len(motifs)  # t is the number of motifs
    profile = [[1 for _ in range(k)] for _ in range(4)]  # Initialize profile with pseudocounts (1)
    
    # Count each symbol in each position
    for motif in motifs: 
        for i in range(len(motif)): 
            j = symbolToNumber(motif[i])  # Get index for the symbol
            profile[j][i] += 1  # Increment count in the profile
    
    # Normalize counts to probabilities
    for i in range(4): 
        for j in range(k): 
            profile[i][j] /= (t + 4)  # Pseudocount division
    
    return profile

# Get consensus string from the profile
def consensus(profile): 
    cons = ""
    for i in range(len(profile[0])): 
        max_prob = 0
        loc = 0
        for j in range(4):  # Find the symbol with the highest probability at position i
            if profile[j][i] > max_prob:
                loc = j
                max_prob = profile[j][i]
        cons += numberToSymbol(loc)  # Add the symbol to the consensus string
    return cons

# Calculate the score of the motifs
def score(motifs): 
    profile = profileForm(motifs)  # Build profile from motifs
    cons = consensus(profile)  # Get the consensus motif
    score = 0 
    for motif in motifs:  # Compare each motif to the consensus
        for i in range(len(motif)): 
            if cons[i] != motif[i]:  # Increase score for mismatches
                score += 1
    return score

# Select an index based on probabilities
def myRandom(dist): 
    s = sum(dist)  # Calculate total sum of probabilities
    i = random.random()  # Generate random number between 0 and 1
    partial = 0.0
    for idx, p in enumerate(dist):  # For each element in dist
        partial += p
        if partial / s >= i:  # If cumulative probability exceeds i
            return idx

# Gibbs sampler for motif finding
def gibbsSampler(dna, k, t, n): 
    bestMotifs = []
    motifs = []
    
    # Initialize random motifs from each string in dna
    for seq in dna:  
        i = random.randint(0, len(seq) - k)  # Select random starting index
        motifs.append(seq[i:i+k])  # Add random k-mer
    
    bestMotifs = motifs[:]  # Copy motifs as bestMotifs
    
    for _ in range(n):  # Iterate n times
        j = random.randint(0, t-1)  # Randomly select a motif to update
        profile = profileForm(motifs[:j] + motifs[j+1:])  # Form profile excluding the jth motif
        r = profileRandom(k, profile, dna[j])  # Get new motif from the profile
        motifs[j] = dna[j][r:r+k]  # Replace the jth motif
        
        if score(motifs) < score(bestMotifs):  # Update best motifs if current motifs are better
            bestMotifs = motifs[:]
    
    return bestMotifs

# Hardcoded input data
k = 6  # Length of k-mers
t = 26  # Number of DNA sequences
n = 500  # Number of iterations for Gibbs sampling

dna = [
    "gcggaagagggcactagcccatgtgagagggcaaggacca",
    "atctttctcttaaaaataacataattcagggccaggatgt",
    "gtcacgagctttatcctacagatgatgaatgcaaatcagc",
    "taaaagataatatcgaccctagcgtggcgggcaaggtgct",
    "gtagattcgggtaccgttcataaaagtacgggaatttcgg",
    "tatacttttaggtcgttatgttaggcgagggcaaaagtca",
    "ctctgccgattcggcgagtgatcgaagagggcaatgcctc",
    "aggatggggaaaatatgagaccaggggagggccacactgc",
    "acacgtctagggctgtgaaatctctgccgggctaacagac",
    "gtgtcgatgttgagaacgtaggcgccgaggccaacgctga",
    "atgcaccgccattagtccggttccaagagggcaactttgt",
    "ctgcgggcggcccagtgcgcaacgcacagggcaaggttta",
    "tgtgttgggcggttctgaccacatgcgagggcaacctccc",
    "gtcgcctaccctggcaattgtaaaacgacggcaatgttcg",
    "cgtattaatgataaagaggggggtaggaggtcaactcttc",
    "aatgcttataacataggagtagagtagtgggtaaactacg",
    "tctgaaccttctttatgcgaagacgcgagggcaatcggga",
    "tgcatgtctgacaacttgtccaggaggaggtcaacgactc",
    "cgtgtcatagaattccatccgccacgcggggtaatttgga",
    "tcccgtcaaagtgccaacttgtgccggggggctagcagct",
    "acagcccgggaatatagacgcgtttggagtgcaaacatac",
    "acgggaagatacgagttcgatttcaagagttcaaaacgtg",
    "cccgataggactaataaggacgaaacgagggcgatcaatg",
    "ttagtacaaacccgctcacccgaaaggagggcaaatacct",
    "agcaaggttcagatatacagccaggggagacctataactc",
    "gtccacgtgcgtatgtactaattgtggagagcaaatcatt"
]

for i in range(len(dna)):
    dna[i] = dna[i].upper()

# Set a random seed for reproducibility
random.seed(42)

# Run Gibbs sampler and get initial best motifs
best = gibbsSampler(dna, k, t, n) 
best_score = score(best)
all_scores = []

# Run Gibbs sampler multiple times (500+ iterations) to find better motifs
for _ in range(10):  
    sample = gibbsSampler(dna, k, t, n)
    sample_score = score(sample)
    all_scores.append(sample_score)
    # print(_ ,sample , " " , sample_score)
    if sample_score < best_score:  # Update best motifs if a better one is found
        best_score = sample_score
        best = sample[:]

# Visualization function
def visualize_scores(dnas, k, t, n, num_runs):
    all_scores = []
    
    for _ in range(num_runs):
        sample = gibbsSampler(dnas, k, t, n)
        score_value = score(sample)
        all_scores.append(score_value)

    # Plot the average score over iterations
    avg_scores = [sum(all_scores) / num_runs]
    
    plt.figure(figsize=(10, 6))
    plt.plot(range(num_runs), all_scores, marker='o', label='Individual Scores')
    plt.axhline(y=avg_scores[0], color='r', linestyle='--', label='Average Score')
    plt.title("Consensus Motif Scores Over Iterations")
    plt.xlabel("Iteration")
    plt.ylabel("Consensus Score")
    plt.legend()
    plt.grid()
    plt.show()

# Call the visualization function
# visualize_scores(dna, k, t, n, 500)

# Output the best motifs and score
plt.figure(figsize=(10, 6))
plt.hist(all_scores, bins=20, alpha=0.7, color='blue', edgecolor='black')
plt.title("Distribution of Consensus Motif Scores")
plt.xlabel("Consensus Score")
plt.ylabel("Frequency")
plt.grid()
plt.show()

print("Best Motifs:")
for motif in best:
    print(motif)
print("Best Score:", best_score)
