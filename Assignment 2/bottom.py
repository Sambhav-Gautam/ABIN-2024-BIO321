def edit_distance_bottom_up(seq1, seq2, ins_cost, del_cost, gap_penalty):
    """
    Bottom-up approach to compute minimum edit distance with gap penalties.

    Args:
    seq1 (str): First DNA sequence.
    seq2 (str): Second DNA sequence.
    ins_cost (int): The cost of an isolated insertion.
    del_cost (int): The cost of an isolated deletion.
    gap_penalty (int): The penalty for consecutive insertions or deletions.

    Returns:
    int: The minimum cost of transforming seq1 into seq2.
    """
    m, n = len(seq1), len(seq2)
    
    # DP table to store the minimum edit distances
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Base cases: one of the sequences is exhausted
    for i in range(1, m + 1):
        dp[i][0] = del_cost + (i - 1) * gap_penalty  # Deleting all characters from seq1
    for j in range(1, n + 1):
        dp[0][j] = ins_cost + (j - 1) * gap_penalty  # Inserting all characters into seq1
    
    # Fill the DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:  # No cost if the characters match
                dp[i][j] = dp[i - 1][j - 1]
            else:
                # Calculate costs for insertion, deletion, and substitution
                if j == 1 or dp[i][j - 1] == dp[i][j - 2] + gap_penalty:
                    insert_cost = dp[i][j - 1] + ins_cost
                else:
                    insert_cost = dp[i][j - 1] + gap_penalty
                
                if i == 1 or dp[i - 1][j] == dp[i - 2][j] + gap_penalty:
                    delete_cost = dp[i - 1][j] + del_cost
                else:
                    delete_cost = dp[i - 1][j] + gap_penalty
                
                substitute_cost = dp[i - 1][j - 1] + 1

                # Take the minimum of the three possible operations
                dp[i][j] = min(insert_cost, delete_cost, substitute_cost)
    
    # Printing the DP matrix
    print("DP Matrix:")
    for row in dp:
        print(row)

    return dp[m][n]

# Example usage
seq1 = "AGTACG"
seq2 = "GTTAC"

ins_cost = 2  # Cost of a single isolated insertion
del_cost = 2  # Cost of a single isolated deletion
gap_penalty = 1  # Reduced cost for consecutive insertions or deletions

min_cost = edit_distance_bottom_up(seq1, seq2, ins_cost, del_cost, gap_penalty)
print(f"Minimum edit distance: {min_cost}")
