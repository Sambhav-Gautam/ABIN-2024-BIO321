from collections import defaultdict

def bwt_transform(s):
    """Perform Burrows-Wheeler Transform on string s."""
    s = s + '$'
    rotations = sorted(s[i:] + s[:i] for i in range(len(s)))
    bwt_result = ''.join(rotation[-1] for rotation in rotations)
    return bwt_result, rotations.index(s)

def reverse_bwt(bwt_string, original_index):
    """Reverse the Burrows-Wheeler Transform to obtain original string."""
    table = [""] * len(bwt_string)
    for _ in bwt_string:
        table = sorted(bwt_string[i] + table[i] for i in range(len(bwt_string)))
    return table[original_index]

def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    """Perform the Smith-Waterman local alignment."""
    # Initialize scoring and traceback matrices
    m, n = len(seq1), len(seq2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    max_score, max_pos = 0, (0, 0)

    # Fill score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete = score_matrix[i - 1][j] + gap
            insert = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(0, match_score, delete, insert)
            
            # Update max score position for traceback
            if score_matrix[i][j] >= max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Trace back to find the optimal local alignment
    align1, align2 = "", ""
    i, j = max_pos
    while score_matrix[i][j] != 0:
        current_score = score_matrix[i][j]
        if current_score == score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif current_score == score_matrix[i - 1][j] + gap:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2

def compress_and_align(seq1, seq2):
    # Compress the first sequence using BWT
    compressed_seq, original_index = bwt_transform(seq1)

    # Perform alignment with the Smith-Waterman algorithm
    align_compressed, align_uncompressed = smith_waterman(compressed_seq, seq2)

    # Decompress the aligned sequence
    decompressed_align = reverse_bwt(align_compressed, original_index)

    return decompressed_align, align_uncompressed

# Example usage
seq1 = "ACGTACGTG"
seq2 = "CGTACGT"

aligned_seq1, aligned_seq2 = compress_and_align(seq1, seq2)
print("Aligned Sequence 1 (Decompressed):", aligned_seq1)
print("Aligned Sequence 2:", aligned_seq2)
