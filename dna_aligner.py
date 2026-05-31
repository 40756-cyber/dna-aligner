def read_fasta(filename):
    with open(filename) as f:
        lines = f.readlines()
    return "".join(line.strip() for line in lines[1:]).upper()

seq1 = read_fasta("gene1.fasta")
seq2 = read_fasta("gene2.fasta")
print("Sequence 1:", seq1)
print("Sequence 2:", seq2)
print("Length 1:", len(seq1))
print("Length 2:", len(seq2))
# Scoring rules
match_score = 2      # reward for two bases matching
mismatch = -1        # penalty for two bases being different
gap = -1             # penalty for inserting a gap

# Build the grid (rows = seq2, columns = seq1)
rows = len(seq2) + 1
cols = len(seq1) + 1

# Fill grid with zeros
grid = [[0] * cols for _ in range(rows)]

# Fill first row and column with gap penalties
for i in range(rows):
    grid[i][0] = i * gap
for j in range(cols):
    grid[0][j] = j * gap

print("Grid initialized!")
print("Grid size:", rows, "x", cols)
# Fill in the rest of the grid
for i in range(1, rows):
    for j in range(1, cols):
        # Check if bases match or mismatch
        if seq1[j-1] == seq2[i-1]:
            diagonal = grid[i-1][j-1] + match_score
        else:
            diagonal = grid[i-1][j-1] + mismatch

        # Three choices — pick the best score
        up   = grid[i-1][j] + gap
        left = grid[i][j-1] + gap

        grid[i][j] = max(diagonal, up, left)

# Print the finished grid
print("\nScoring grid:")
for row in grid:
    print([str(x).rjust(3) for x in row])
# Traceback to find the best alignment
# Traceback — walk backwards through the grid
align1 = ""
align2 = ""
i = rows - 1
j = cols - 1

while i > 0 or j > 0:
    if i > 0 and j > 0 and grid[i][j] == grid[i-1][j-1] + (match_score if seq1[j-1] == seq2[i-1] else mismatch):
        align1 = seq1[j-1] + align1
        align2 = seq2[i-1] + align2
        i -= 1
        j -= 1
    elif i > 0 and grid[i][j] == grid[i-1][j] + gap:
        align1 = "-" + align1
        align2 = seq2[i-1] + align2
        i -= 1
    else:
        align1 = seq1[j-1] + align1
        align2 = "-" + align2
        j -= 1

# Print the alignment
print("\nAlignment result:")
print("Seq1:", align1)
match_line = "".join("|" if a == b else "X" for a, b in zip(align1, align2))
print("     ", match_line)
print("Seq2:", align2)
print("\nAlignment score:", grid[rows-1][cols-1])
# Show alignment with line breaks every 60 characters
print("\nFull alignment:")
chunk = 60
for i in range(0, len(align1), chunk):
    a1_chunk = align1[i:i+chunk]
    a2_chunk = align2[i:i+chunk]
    match_chunk = "".join("|" if a == b else "X" for a, b in zip(a1_chunk, a2_chunk))
    print(f"\nPos {i+1}-{i+len(a1_chunk)}")
    print("Human:", a1_chunk)
    print("      ", match_chunk)
    print("Mouse:", a2_chunk)