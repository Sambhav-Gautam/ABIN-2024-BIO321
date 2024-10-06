import math

def hybrid_global_local_alignment(A, B, match_score, mismatch_penalty, gap_open, gap_extend, local_threshold):
    m, n = len(A), len(B)
    
    # Initialize matrices for scores and directions
    M = [[-math.inf] * (n + 1) for _ in range(m + 1)]
    I = [[-math.inf] * (n + 1) for _ in range(m + 1)]
    D = [[-math.inf] * (n + 1) for _ in range(m + 1)]
    traceback = [[None] * (n + 1) for _ in range(m + 1)]
    
    # Base case: starting point
    M[0][0] = 0
    I[0][0] = D[0][0] = -math.inf
    
    # Initialize first row for insertions
    for j in range(1, n + 1):
        I[0][j] = gap_open + (j - 1) * gap_extend
        traceback[0][j] = 'I'
    
    # Initialize first column for deletions
    for i in range(1, m + 1):
        D[i][0] = gap_open + (i - 1) * gap_extend
        traceback[i][0] = 'D'

    # Fill matrices with local-global hybrid logic
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Match or mismatch score
            score = match_score if A[i-1] == B[j-1] else mismatch_penalty
            
            # Update M[i][j]
            M[i][j] = max(
                (M[i-1][j-1] + score, 'M'),
                (I[i-1][j-1] + score, 'I'),
                (D[i-1][j-1] + score, 'D'),
                key=lambda x: x[0]
            )[0]
            
            # Apply local alignment reset if below threshold
            if M[i][j] < local_threshold:
                M[i][j] = 0
            
            # Update I[i][j] (insertion)
            I[i][j] = max(
                M[i][j-1] + gap_open,
                I[i][j-1] + gap_extend
            )
            
            # Local alignment reset
            if I[i][j] < local_threshold:
                I[i][j] = 0

            # Update D[i][j] (deletion)
            D[i][j] = max(
                M[i-1][j] + gap_open,
                D[i-1][j] + gap_extend
            )
            
            # Local alignment reset
            if D[i][j] < local_threshold:
                D[i][j] = 0
    
            # Store the best direction in traceback
            if M[i][j] >= I[i][j] and M[i][j] >= D[i][j]:
                traceback[i][j] = 'M'
            elif I[i][j] >= M[i][j] and I[i][j] >= D[i][j]:
                traceback[i][j] = 'I'
            else:
                traceback[i][j] = 'D'

    # Get the final global-local alignment score
    alignment_score = max(M[m][n], I[m][n], D[m][n])
    
    # Traceback to get the alignment
    i, j = m, n
    aligned_A, aligned_B = [], []
    
    while i > 0 or j > 0:
        if traceback[i][j] == 'M':
            aligned_A.append(A[i-1])
            aligned_B.append(B[j-1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 'I':
            aligned_A.append('-')
            aligned_B.append(B[j-1])
            j -= 1
        elif traceback[i][j] == 'D':
            aligned_A.append(A[i-1])
            aligned_B.append('-')
            i -= 1
    
    aligned_A = ''.join(reversed(aligned_A))
    aligned_B = ''.join(reversed(aligned_B))
    
    return alignment_score, aligned_A, aligned_B

# Example Usage
# A = "ACGTACGT"
# B = "AGTACG"
A = "DONE"
B = "REDO"
match_score = 2
mismatch_penalty = -1
gap_open = -2
gap_extend = -0.5
local_threshold = -3

score, aligned_A, aligned_B = hybrid_global_local_alignment(A, B, match_score, mismatch_penalty, gap_open, gap_extend, local_threshold)

print(f"Alignment Score: {score}")
print(f"Aligned A: {aligned_A}")
print(f"Aligned B: {aligned_B}")