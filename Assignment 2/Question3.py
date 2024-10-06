def compute_suffix_array(input_text):
    len_text = len(input_text)
    # Array of tuples to store rotations and their indexes
    suff = [(i, input_text[i:]) for i in range(len_text)]
    # Sorts rotations lexicographically
    suff.sort(key=lambda x: x[1])
    # Stores the indexes of sorted rotations
    suffix_arr = [i for i, _ in suff]
    return suffix_arr

def find_last_char(input_text, suffix_arr):
    n = len(input_text)
    # BWT result
    bwt_arr = ""
    for i in range(n):
        # Find last character of each rotation
        j = (suffix_arr[i] - 1) % n
        bwt_arr += input_text[j]
    return bwt_arr

def build_last_to_first(bwt_arr):
    # Maps each position in the BWT to the original suffix array order
    sorted_bwt = sorted([(char, i) for i, char in enumerate(bwt_arr)])
    return [i[1] for i in sorted_bwt]

def pattern_match(pattern, bwt_arr, suffix_arr, last_to_first):
    # Finds exact matches of 'pattern' in the BWT-transformed text
    top = 0
    bottom = len(bwt_arr) - 1
    for char in reversed(pattern):
        top_index = next((i for i in range(top, bottom + 1) if bwt_arr[i] == char), None)
        bottom_index = next((i for i in range(bottom, top - 1, -1) if bwt_arr[i] == char), None)

        # Check if character not found in the specified range
        if top_index is None or bottom_index is None:
            return []

        # Update top and bottom to the corresponding last-to-first positions
        top = last_to_first[top_index]
        bottom = last_to_first[bottom_index]

    return [suffix_arr[i] for i in range(top, bottom + 1)]

# Driver code to test the functions
input_text = "AG"  # DNA sequence with end marker $
input_text+= "$"
suffix_arr = compute_suffix_array(input_text)
bwt_arr = find_last_char(input_text, suffix_arr)
last_to_first = build_last_to_first(bwt_arr)

patterns = ["AGT", "AAG"]
print("Input sequence:", input_text)
print("Burrows-Wheeler Transform:", bwt_arr)
print("Suffix Array:", suffix_arr)
print("Last-to-First Mapping:", last_to_first)

# Finding positions of each pattern in the original sequence
for pattern in patterns:
    match_positions = pattern_match(pattern, bwt_arr, suffix_arr, last_to_first)
    print(f"Positions of pattern '{pattern}':", match_positions)
