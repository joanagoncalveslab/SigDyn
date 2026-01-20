
def aggregate_neighbors(input_files, out_file, K):
    # Dictionary to store the top-K scores for each cell using a min-heap
    scores = defaultdict(list)

    # Process each input file to collect similarity scores
    for file in input_files:
        print(f"Processing {file}...")
        with open(file, encoding="utf-8") as f:
            header = next(f).strip().split("\t")
            assert header[0] == "Score", f"Unexpected header: {header[0]}"

            for line in f:
                cols = line.strip().split("\t")
                cell1 = cols[0]

                for i in range(1, len(cols)):
                    cell2 = header[i]
                    if cell1 == cell2:
                        continue  # Skip if comparing cell to itself

                    # Extract the score part
                    score = float(cols[i].split(':')[0])

                    # Use a min-heap to maintain only the top K neighbors
                    for cell, neighbor in [(cell1, cell2), (cell2, cell1)]:
                        heapq.heappush(scores[cell], (score, neighbor))
                        if len(scores[cell]) > K:
                            heapq.heappop(scores[cell])

    # Write the aggregated neighbors file
    with open(out_file, 'w', encoding="utf-8") as handle:
        # Create header for the output file
        neighbors_header = [f"N{n + 1:02d}" for n in range(K)]
        header = ["Cell"] + neighbors_header
        print("\t".join(header), file=handle)

        for cell in sorted(scores):
            # Sort the top K scores in descending order
            top_neighbors = sorted(scores[cell], key=lambda x: -x[0])
            neighbor_cells = [neighbor for _, neighbor in top_neighbors]

            # Pad with empty strings if fewer than K neighbors
            if len(neighbor_cells) < K:
                neighbor_cells.extend([""] * (K - len(neighbor_cells)))

            print("\t".join([cell] + neighbor_cells), file=handle)

if __name__ == '__main__':
    import glob
    from collections import defaultdict
    import heapq  # For maintaining a fixed-size top-K list

    input_files = snakemake.input
    out_file = snakemake.output[0]
    K = snakemake.params.K

    aggregate_neighbors(input_files, out_file, K)

