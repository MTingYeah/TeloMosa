#!/usr/bin/env python3

import gzip
import argparse


MOTIF = "TTAGGG"
WINDOW_REPEAT = 2
MIN_CONTINUOUS_REPEAT = 4
MAX_MISMATCH = 2
STEP = 6

WINDOW_SIZE = len(MOTIF) * WINDOW_REPEAT
MIN_TELOMERIC_BP = len(MOTIF) * MIN_CONTINUOUS_REPEAT


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode)
    return open(path, mode=mode)


def revcomp(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def rotations(s):
    return [s[i:] + s[:i] for i in range(len(s))]


def hamming(a, b):
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def build_templates():
    templates = []

    fwd_long = MOTIF * WINDOW_REPEAT
    rev_long = revcomp(MOTIF) * WINDOW_REPEAT

    templates.extend(rotations(fwd_long))
    templates.extend(rotations(rev_long))

    return set(templates)


TEMPLATES = build_templates()


def is_telomeric_window(window):
    for template in TEMPLATES:
        if hamming(window, template) <= MAX_MISMATCH:
            return True
    return False


def estimate_telomeric_repeat_length(seq):
    seq = seq.upper()

    if len(seq) < WINDOW_SIZE:
        return 0, "NA", "NA"

    hit_positions = []

    for i in range(0, len(seq) - WINDOW_SIZE + 1, STEP):
        window = seq[i:i + WINDOW_SIZE]
        if is_telomeric_window(window):
            hit_positions.append(i)

    if not hit_positions:
        return 0, "NA", "NA"

    blocks = []
    start = hit_positions[0]
    prev = hit_positions[0]

    for pos in hit_positions[1:]:
        if pos <= prev + STEP:
            prev = pos
        else:
            blocks.append((start, prev + WINDOW_SIZE))
            start = pos
            prev = pos

    blocks.append((start, prev + WINDOW_SIZE))

    valid_blocks = [
        (start, end) for start, end in blocks
        if (end - start) >= MIN_TELOMERIC_BP
    ]

    if not valid_blocks:
        return 0, "NA", "NA"

    longest_block = max(valid_blocks, key=lambda x: x[1] - x[0])
    tel_start, tel_end = longest_block
    tel_bp = tel_end - tel_start

    return tel_bp, tel_start, tel_end


def read_fastq_chunks(handle):
    while True:
        name = handle.readline()
        if not name:
            return

        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()

        if not qual:
            return

        yield name, seq, plus, qual


def clean_read_id(name_line):
    return name_line.strip().split()[0].replace("@", "")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Estimate telomeric repeat length from ONT FASTQ reads using fixed criteria: "
            "12-bp telomeric windows and minimum 4 continuous TTAGGG repeats."
        )
    )

    parser.add_argument(
        "input",
        help="Input FASTQ file (.fastq, .fq, .fastq.gz, or .fq.gz)"
    )

    parser.add_argument(
        "-o", "--out-prefix",
        required=True,
        help="Output prefix"
    )

    args = parser.parse_args()

    out_tsv = args.out_prefix + ".telomeric_repeat_length.tsv"

    total = 0
    reported = 0

    with open_maybe_gz(args.input, "rt") as fin, open(out_tsv, "w") as fout:
        fout.write(
            "read_id\tread_len\ttelomeric_repeat_bp\ttelomeric_repeat_kb\t"
            "telomeric_fraction\ttel_start\ttel_end\t"
            "motif\twindow_bp\tmin_continuous_repeat\tmax_mismatch\tstep\n"
        )

        for chunk in read_fastq_chunks(fin):
            total += 1

            name, seq_line, plus, qual = chunk
            read_id = clean_read_id(name)
            seq = seq_line.strip()
            read_len = len(seq)

            tel_bp, tel_start, tel_end = estimate_telomeric_repeat_length(seq)

            if tel_bp >= MIN_TELOMERIC_BP:
                reported += 1
                tel_kb = tel_bp / 1000
                tel_fraction = tel_bp / read_len if read_len > 0 else 0

                fout.write(
                    f"{read_id}\t{read_len}\t{tel_bp}\t{tel_kb:.6f}\t"
                    f"{tel_fraction:.6f}\t{tel_start}\t{tel_end}\t"
                    f"{MOTIF}\t{WINDOW_SIZE}\t{MIN_CONTINUOUS_REPEAT}\t"
                    f"{MAX_MISMATCH}\t{STEP}\n"
                )

    print(f"Total reads processed: {total}")
    print(f"Reads with telomeric repeat block >= {MIN_TELOMERIC_BP} bp: {reported}")
    print(f"Window size: {WINDOW_SIZE} bp")
    print(f"Minimum continuous telomeric repeat length: {MIN_TELOMERIC_BP} bp")
    print(f"Output: {out_tsv}")


if __name__ == "__main__":
    main()
