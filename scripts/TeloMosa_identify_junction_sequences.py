#!/usr/bin/env python3
import gzip
import sys

G_MOTIF = "TTAGGG"
C_MOTIF = "TAACCC"

def rc(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def open_file(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def count_canonical_hits(seq: str) -> int:
    hits = 0
    for i in range(len(seq) - 5):
        k = seq[i:i+6]
        if k == G_MOTIF or k == C_MOTIF:
            hits += 1
    return hits

def choose_orientation(seq: str) -> str:
    s1 = seq
    s2 = rc(seq)
    return s1 if count_canonical_hits(s1) >= count_canonical_hits(s2) else s2

def iter_hits(seq: str):
    """Yield (pos, type) for canonical 6-mer hits only."""
    for i in range(len(seq) - 5):
        k = seq[i:i+6]
        if k == G_MOTIF:
            yield i, "G"
        elif k == C_MOTIF:
            yield i, "C"

def build_runs(hits, min_hits=3, max_gap=6):
    
    if not hits:
        return []

    runs = []
    cur_t = hits[0][1]
    cur_start = hits[0][0]
    cur_last = hits[0][0]
    cur_n = 1

    for pos, t in hits[1:]:
        if t == cur_t and (pos - cur_last) <= max_gap:
            cur_last = pos
            cur_n += 1
        else:
            if cur_n >= min_hits:
                runs.append((cur_t, cur_start, cur_last + 6, cur_n))
            cur_t = t
            cur_start = pos
            cur_last = pos
            cur_n = 1

    if cur_n >= min_hits:
        runs.append((cur_t, cur_start, cur_last + 6, cur_n))

    return runs

def main():
    if len(sys.argv) < 2:
        print(
            "Usage: junction_point_recognition_canonical.py <reads.fastq|reads.fastq.gz> [min_hits] [max_gap] [keep_empty(0/1)] [max_junction_len]",
            file=sys.stderr
        )
        sys.exit(1)

    fq = sys.argv[1]
    min_hits = int(sys.argv[2]) if len(sys.argv) >= 3 else 3
    max_gap = int(sys.argv[3]) if len(sys.argv) >= 4 else 6
    keep_empty = int(sys.argv[4]) if len(sys.argv) >= 5 else 1
    max_junction_len = int(sys.argv[5]) if len(sys.argv) >= 6 else 12

    total_reads = 0
    total_junctions = 0

    print("ReadID\tDirection\tJunctionSeq\tJunctionLen")

    with open_file(fq) as f:
        read_id = None
        for i, line in enumerate(f):
            mod = i % 4

            if mod == 0:
                read_id = line.strip().split()[0][1:]

            elif mod == 1:
                total_reads += 1
                seq0 = line.strip().upper()
                seq = seq0

                hits = list(iter_hits(seq))
                if len(hits) < min_hits:
                    continue

                runs = build_runs(hits, min_hits=min_hits, max_gap=max_gap)
                if len(runs) < 2:
                    continue

                for (t1, s1, e1, n1), (t2, s2, e2, n2) in zip(runs, runs[1:]):
                    if t1 == t2:
                        continue

                    jseq = seq[e1:s2]

                    if (not jseq) and (not keep_empty):
                        continue

                    if len(jseq) > max_junction_len:
                        continue

                    direction = f"{t1}to{t2}"
                    js = jseq if jseq != "" else "<EMPTY>"
                    jlen = len(jseq)

                    print(f"{read_id}\t{direction}\t{js}\t{jlen}")
                    total_junctions += 1

    print(
        f"[QC] reads={total_reads} junctions_reported={total_junctions} "
        f"(canonical TTAGGG/TAACCC only, max_junction_len={max_junction_len})",
        file=sys.stderr
    )

if __name__ == "__main__":
    main()

    