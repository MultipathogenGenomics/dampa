from collections import defaultdict
import subprocess
from Bio import SeqIO
import re
import os

def runblast(targets,blastout,ident):
    terminal_command = (f"makeblastdb -in {targets} -dbtype nucl -logfile /dev/null")
    finished_process = subprocess.run(terminal_command, shell=True)#, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if finished_process.returncode != 0:
        print(f'\nERROR: makeblastdb terminated with errors while indexing target seqs (Error code: {finished_process.returncode}).\n')
        print(finished_process.stderr)
        print(finished_process.stdout)
        exit(1)
    terminal_command = (f"blastn -db {targets} -query {targets} -task megablast -perc_identity {ident} -dust no -soft_masking false -outfmt '6 qseqid sseqid pident qstart qend qlen length' > {blastout}")
    finished_process = subprocess.run(terminal_command, shell=True)#, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if finished_process.returncode != 0:
        print(f'\nERROR: blastn terminated with errors while aligning probes to target seqs (Error code: {finished_process.returncode}).\n')
        exit(1)
    for ext in [".nhr", ".nin", ".nsq", ".nal", ".ndb", ".nts",".not",".njs",".nto",".ntf"]:
        db_file = targets + ext
        if os.path.exists(db_file):
            os.remove(db_file)

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for s,e in intervals[1:]:
        ls,le = merged[-1]
        if s <= le+1:
            merged[-1] = (ls, max(le,e))
        else:
            merged.append((s,e))
    return merged

def compute_coverage(seq_id, kept,lengths,alignments):
    ivals = [(qs, qe) for s, qs, qe in alignments.get(seq_id, []) if s in kept]
    merged = merge_intervals(ivals)
    covered = sum(e - s + 1 for s,e in merged)
    return covered / lengths[seq_id] if lengths[seq_id] > 0 else 0

def has_clean_stretch(seq, min_len=90):
    pattern = rf"[ACGTacgt]{{{min_len},}}"
    return bool(re.search(pattern, seq))

def union_coverage(targets,outtargets,min_ident=99.0,min_covered_frac=99.0,non_n_stretch=90):

    base = outtargets.replace(".fasta","").replace(".fa","")
    targets_nfilt = base+"_nfilt.fasta"
    tsv = base+"all_vs_all.tsv"

    intargets = {s.id:s for s in targets}

    usabletargets = [s for s in intargets.values() if has_clean_stretch(str(s.seq), non_n_stretch)]
    nfilt = len(intargets) - len(usabletargets)
    usable_dict = {s.id: s for s in usabletargets}
    if not usabletargets:
        return 0,0,0,nfilt,set(),[]
    SeqIO.write(usabletargets, targets_nfilt, "fasta")

    runblast(targets_nfilt, tsv, min_ident)

    alignments = defaultdict(list)
    with open(tsv) as f:
        for line in f:
            q, s, pident, qstart, qend, qlen_i, alen = line.strip().split('\t')
            if q == s:
                continue
            pident = float(pident)
            if pident < min_ident:
                continue
            qs, qe = int(qstart), int(qend)
            if qs > qe:
                qs, qe = qe, qs
            alignments[q].append((s, qs, qe))

    seq_order = sorted(usable_dict.keys(), key=lambda k: len(usable_dict[k].seq))

    lengths = {s.id: len(s.seq) for s in usable_dict.values()}

    kept = set(seq_order)
    removed = set()

    changed = True
    while changed:
        changed = False
        for seq_id in sorted(kept, key=lambda x: lengths[x]):  # smallest first
            frac_covered = compute_coverage(seq_id, kept - {seq_id}, lengths, alignments)
            len_multicovered = lengths[seq_id] * frac_covered
            minmulticovered = lengths[seq_id] * min_covered_frac

            if len_multicovered > minmulticovered:
                kept.remove(seq_id)
                removed.add(seq_id)
                changed = True
                break


    keep = []
    for k in kept:
        keepseq = usable_dict[k]
        keep.append(keepseq)
    # SeqIO.write(keep, outtargets, "fasta")
    if os.path.exists(targets_nfilt):
        os.remove(targets_nfilt)
    if os.path.exists(tsv):
        os.remove(tsv)
    return len(lengths), len(kept), len(removed), nfilt,kept,keep