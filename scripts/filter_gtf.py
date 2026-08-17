import sys

# Load bad transcript IDs
bad_tx_ids = set()
with open(sys.argv[1], 'r') as f:
    for line in f:
        bad_tx_ids.add(line.strip())

print(f"[INFO] Loaded {len(bad_tx_ids)} bad transcript IDs to filter", file=sys.stderr)

# Read GTF and track which genes to keep
# A gene is kept only if ALL its transcripts are clean
gene_to_transcripts = {}  # gene_id -> set of transcript_ids
bad_genes = set()

# First pass: identify which genes have bad transcripts
with open(sys.argv[2], 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        # Parse attributes
        attrs = {}
        for attr in fields[8].split(';'):
            attr = attr.strip()
            if not attr:
                continue
            if ' ' in attr:
                key, val = attr.split(' ', 1)
                attrs[key] = val.strip('"')

        gene_id = attrs.get('gene_id')
        tx_id = attrs.get('transcript_id')

        if gene_id and tx_id:
            if gene_id not in gene_to_transcripts:
                gene_to_transcripts[gene_id] = set()
            gene_to_transcripts[gene_id].add(tx_id)

            # If this transcript is bad, mark the whole gene as bad
            if tx_id in bad_tx_ids:
                bad_genes.add(gene_id)

print(f"[INFO] Found {len(bad_genes)} genes with bad transcripts", file=sys.stderr)

# Second pass: write clean GTF (skip bad genes entirely)
filtered_count = 0
kept_count = 0

with open(sys.argv[2], 'r') as f_in, open(sys.argv[3], 'w') as f_out:
    for line in f_in:
        if line.startswith('#'):
            f_out.write(line)
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            f_out.write(line)
            continue

        # Parse attributes
        attrs = {}
        for attr in fields[8].split(';'):
            attr = attr.strip()
            if not attr:
                continue
            if ' ' in attr:
                key, val = attr.split(' ', 1)
                attrs[key] = val.strip('"')

        gene_id = attrs.get('gene_id')

        # Skip entire gene if it has any bad transcripts
        if gene_id in bad_genes:
            filtered_count += 1
            continue

        f_out.write(line)
        kept_count += 1

print(f"[INFO] Filtered {filtered_count} lines", file=sys.stderr)
print(f"[INFO] Kept {kept_count} lines", file=sys.stderr)