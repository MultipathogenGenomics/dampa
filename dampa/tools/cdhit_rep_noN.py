from Bio import SeqIO
from collections import defaultdict
import sys
def longest_or_fewest_ns_representatives(fasta_file, clstr_file, output_file):
    # Load sequences
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    shorttolong = {x[:19]:x for x in seqs.keys()}
    # Parse .clstr file
    clusters = defaultdict(list)
    cluster_id = -1
    with open(clstr_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id += 1
            else:
                seq_id = line.split('>')[1].split('...')[0]
                clusters[cluster_id].append(seq_id)

    # Choose representative per cluster
    representatives = []
    for cluster, members in clusters.items():
        seq_objs = [seqs[shorttolong[sid]] for sid in members if shorttolong[sid] in seqs]
        if len(seq_objs) == 0:
            print("no sequences in cluster", cluster, "skipping")
            continue

        # Try to get N-free sequences
        nfree = [s for s in seq_objs if 'N' not in s.seq.upper()]
        if nfree:
            best = max(nfree, key=lambda s: len(s.seq))
        else:
            # Fallback: fewest Ns, tie-breaker is longest
            best = min(seq_objs, key=lambda s: (s.seq.upper().count('N'), -len(s.seq)))

        representatives.append(best)

    # Write results
    with open(output_file, "w") as out_f:
        SeqIO.write(representatives, out_f, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python cdhit_rep_noN.py <fasta_file> <clstr_file> <output_file>")
        sys.exit(1)
    longest_or_fewest_ns_representatives(sys.argv[1], sys.argv[2], sys.argv[3])
