def edit_distance_top_down(seq1, seq2, ins_cost, del_cost, gap_penalty):
    """
    Top-down recursive approach with memoization to compute minimum edit distance with gap penalties.
    Prints the computed DP matrix in a formatted grid.
    
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

    # Memoization table to store the already computed values
    memo = {}

    # Create a DP matrix (initialized with None to indicate unvisited cells)
    dp_matrix = [[None for _ in range(n + 1)] for _ in range(m + 1)]

    def helper(i, j, prev_op):
        """
        Recursive helper function with memoization to compute minimum edit distance.

        Args:
        i (int): Current index in seq1.
        j (int): Current index in seq2.
        prev_op (str): Previous operation ("insert", "delete", or None) to check for gap penalty.

        Returns:
        int: The minimum edit distance for the given indices.
        """
        # Base cases: one of the sequences is exhausted
        if i == len(seq1):
            cost = (gap_penalty) * (len(seq2) - j )
            dp_matrix[i][j] = cost  # Store the value in the matrix
            return cost
        if j == len(seq2):
            cost = (gap_penalty) * (len(seq1) - i )
            dp_matrix[i][j] = cost  # Store the value in the matrix
            return cost

        # Check if the result is already computed
        if (i, j, prev_op) in memo:
            return memo[(i, j, prev_op)]

        # If characters match, no operation needed
        if seq1[i] == seq2[j]:
            result = helper(i + 1, j + 1, None)
        else:
            # Consider three operations: insert, delete, and substitute
            # Insert operation (add character from seq2 to seq1)
            if prev_op == "insert":
                insert_cost = gap_penalty + helper(i, j + 1, "insert")
            else:
                insert_cost = ins_cost + helper(i, j + 1, "insert")

            # Delete operation (remove character from seq1)
            if prev_op == "delete":
                delete_cost = gap_penalty + helper(i + 1, j, "delete")
            else:
                delete_cost = del_cost + helper(i + 1, j, "delete")

            # Substitute operation (replace character in seq1 with character from seq2)
            substitute_cost = 1 + helper(i + 1, j + 1, None)

            # Choose the minimum of the three operations
            result = min(insert_cost, delete_cost, substitute_cost)

        # Memoize and store the result in the matrix
        memo[(i, j, prev_op)] = result
        dp_matrix[i][j] = result
        return result

    # Start recursion with no previous operation and indices at the start of both sequences
    result = helper(0, 0, None)

    # Printing the DP matrix with formatting
    print("\nMemoized States (DP Matrix):")
    print("   " + "   ".join(f"{ch2}" for ch2 in (" ",) + tuple(seq2)))  # Header with seq2 characters
    for i, row in enumerate(dp_matrix):
        seq1_char = seq1[i-1] if i > 0 else " "  # Add seq1 characters along rows
        formatted_row = [f"{x if x is not None else ' '}".rjust(3) for x in row]  # Right-align values in the row
        print(f"{seq1_char}  " + " ".join(formatted_row))

    return result

# Example usage
seq1 = "GATTACA"
seq2 = "GCATGCU"

ins_cost = 1  # Cost of a single isolated insertion
del_cost = 1  # Cost of a single isolated deletion
gap_penalty = 0.5  # Reduced cost for consecutive insertions or deletions

min_cost = edit_distance_top_down(seq1, seq2, ins_cost, del_cost, gap_penalty)
print(f"\nMinimum edit distance: {min_cost}")
