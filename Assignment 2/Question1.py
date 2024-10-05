def edit_distance_top_down(seq1, seq2, ins_cost, del_cost, gap_penalty):
    """
    Top-down recursive approach with memoization to compute minimum edit distance with gap penalties.

    Args:
    seq1 (str): First DNA sequence.
    seq2 (str): Second DNA sequence.
    ins_cost (int): The cost of an isolated insertion.
    del_cost (int): The cost of an isolated deletion.
    gap_penalty (int): The penalty for consecutive insertions or deletions.

    Returns:
    int: The minimum cost of transforming seq1 into seq2.
    """
    # Memoization table to store the already computed values
    memo = {}

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
            return (ins_cost) + (gap_penalty) * (len(seq2) - j - 1)
        if j == len(seq2):
            return (del_cost) + (gap_penalty) * (len(seq1) - i - 1)

        # Check if the result is already computed
        if (i, j, prev_op) in memo:
            return memo[(i, j, prev_op)]

        # If characters match, no operation needed
        if seq1[i] == seq2[j]:
            memo[(i, j, prev_op)] = helper(i + 1, j + 1, None)
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
            memo[(i, j, prev_op)] = min(insert_cost, delete_cost, substitute_cost)

        return memo[(i, j, prev_op)]

    # Start recursion with no previous operation and indices at the start of both sequences
    return helper(0, 0, None)

# Example usage with bigger sequences:
seq1 = "ABC"
seq2 = "ABCDE"
ins_cost = 2  # Cost of a single isolated insertion
del_cost = 2  # Cost of a single isolated deletion
gap_penalty = 1  # Reduced cost for consecutive insertions or deletions

min_cost = edit_distance_top_down(seq1, seq2, ins_cost, del_cost, gap_penalty)
print(f"Minimum edit distance : {min_cost}")
