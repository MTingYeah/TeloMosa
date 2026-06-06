#!/usr/bin/env python3

import argparse
import gzip
import re
import logging
import numpy as np
from Bio import SeqIO

logging.basicConfig(level=logging.ERROR)


expectedTeloCompositionQ = [
    ["GGG", 3 / 6],
]

expectedTeloCompositionP = [
    ["CCC", expectedTeloCompositionQ[0][1]],
]

teloNPTeloCompositionGStrand = [
    ["(?!GGG)[ATC]{2}.GGG|(?!AAA)[AT]{3}AAA|TTAGG.", 6 / 6, 6]
]

teloNPTeloCompositionCStrand = [
    ["CTTCTT|CCTGG|CCC(?!CCC)[ATCG]{3}", 6 / 6, 6]
]

areaDiffsThreshold = 0.2

errorReturns = {
    "init": -1,
    "fusedRead": -10,
    "strandType": -20,
    "seqNotFound": -1000,
}

FUSED_READ_CODE = errorReturns["fusedRead"]


def is_regex_pattern(input_string):
    return bool(re.search(r"[^a-zA-Z]", input_string))


def getGraphArea(offsets, targetColumn, windowSize):
    data = np.array(offsets)
    transposed_data = data.T
    areaList = []
    row = transposed_data[targetColumn, :]

    for i in range(0, len(row) - windowSize, 1):
        area = row[i:i + windowSize].sum()
        areaList.append(area / windowSize)

    return areaList


def validate_seq_teloWindow(seq, teloWindow):
    if not isinstance(teloWindow, int) or teloWindow < 6:
        raise ValueError("teloWindow should be an int greater than or equal to 6")

    if len(seq) < teloWindow:
        raise Warning("sequence length is less than teloWindow")


def validate_parameters(
    seq,
    isGStrand,
    composition,
    teloWindow=100,
    windowStep=6,
    plateauDetectionThreshold=-15,
    changeThreshold=-5,
    targetPatternIndex=-1,
    nucleotideGraphAreaWindowSize=500,
    showGraphs=False,
):
    validate_seq_teloWindow(seq, teloWindow)

    if not isinstance(isGStrand, bool) and not isinstance(isGStrand, np.bool_):
        raise ValueError("isGStrand should be a boolean, or numpy boolean")

    if not isinstance(windowStep, int) or windowStep < 1:
        raise ValueError("windowStep should be an int greater than or equal to 1")

    if not isinstance(targetPatternIndex, int) or targetPatternIndex > len(composition):
        raise ValueError(
            "targetPatternIndex should be an int and within the range of the composition list"
        )

    if not isinstance(nucleotideGraphAreaWindowSize, int) or nucleotideGraphAreaWindowSize < 1:
        raise ValueError(
            "nucleotideGraphAreaWindowSize should be an int greater than or equal to 1"
        )

    if not isinstance(showGraphs, bool):
        raise ValueError("showGraphs should be a boolean")


def getIsGStrandFromSeq(
    seq,
    GStrandPatternIn,
    CStrandPatternIn,
    searchStrandRepeats=4,
    minTeloCountDiff=1,
    fusedReadTeloRepeatThreshold=20,
    maxTelomereGap=150,
):
    CStrandPattern = CStrandPatternIn[0]
    GStrandPattern = GStrandPatternIn[0]

    boundaryReg = ".{0,6}"

    CStrandPattern = ("(" + CStrandPattern + ")" + boundaryReg) * searchStrandRepeats
    GStrandPattern = ("(" + GStrandPattern + ")" + boundaryReg) * searchStrandRepeats

    CStrandPattern = CStrandPattern[:-len(boundaryReg)]
    GStrandPattern = GStrandPattern[:-len(boundaryReg)]

    seq_upper = str(seq).upper()

    allCStrands = [match for match in re.finditer(CStrandPattern, seq_upper)]

    cStrandCount = 0
    cStrandMatchLengths = 0
    lastMatch = 0

    for match in allCStrands:
        if lastMatch == 0:
            lastMatch = match.start()

        if match.start() - lastMatch <= maxTelomereGap:
            cStrandCount += 1
            cStrandMatchLengths += match.end() - match.start()

        lastMatch = match.start()

    allGStrands = [match for match in re.finditer(GStrandPattern, seq_upper)]

    gStrandCount = 0
    gStrandMatchLengths = 0
    lastMatch = 0

    for match in allGStrands[::-1]:
        if lastMatch == 0:
            lastMatch = match.start()

        if lastMatch - match.start() <= maxTelomereGap:
            gStrandCount += 1
            gStrandMatchLengths += match.end() - match.start()

        lastMatch = match.start()

    if max(cStrandCount, gStrandCount) <= minTeloCountDiff and abs(cStrandCount - gStrandCount) <= minTeloCountDiff:
        if cStrandMatchLengths > gStrandMatchLengths:
            return False
        elif gStrandMatchLengths > cStrandMatchLengths:
            return True

        return errorReturns["strandType"]

    if (
        min(cStrandCount, gStrandCount) > 100
        or abs(cStrandCount - gStrandCount) <= 0.6 * min(cStrandCount, gStrandCount)
        or min(cStrandCount, gStrandCount) >= fusedReadTeloRepeatThreshold
    ):
        return errorReturns["fusedRead"]

    if cStrandCount > gStrandCount:
        return False
    else:
        return True


def getTeloBoundary(
    seq,
    isGStrand=None,
    compositionGStrand=None,
    compositionCStrand=None,
    teloWindow=100,
    windowStep=6,
    changeThreshold=-20,
    plateauDetectionThreshold=-50,
    targetPatternIndex=-1,
    nucleotideGraphAreaWindowSize=500,
    returnLastDiscontinuity=False,
    secondarySearch=False,
):
    if compositionGStrand is None:
        compositionGStrand = expectedTeloCompositionQ

    if compositionCStrand is None:
        compositionCStrand = expectedTeloCompositionP

    seq = str(seq)

    boundaryPoint = -1
    ntOffsets = []
    graphAreaWindowSize = int(nucleotideGraphAreaWindowSize / windowStep)

    try:
        validate_seq_teloWindow(seq, teloWindow)
    except Warning:
        return errorReturns["init"]

    if isGStrand is None:
        isGStrand = getIsGStrandFromSeq(
            seq,
            compositionGStrand[targetPatternIndex],
            compositionCStrand[targetPatternIndex],
        )

        if isGStrand < 0 or (
            not isinstance(isGStrand, bool)
            and not isinstance(isGStrand, np.bool_)
        ):
            if isGStrand == errorReturns["fusedRead"]:
                return errorReturns["fusedRead"]

            return errorReturns["strandType"]

    composition = compositionGStrand if isGStrand else compositionCStrand

    try:
        validate_parameters(
            seq,
            isGStrand,
            composition,
            teloWindow,
            windowStep,
            plateauDetectionThreshold,
            changeThreshold,
            targetPatternIndex,
            nucleotideGraphAreaWindowSize,
            False,
        )
    except Warning:
        return errorReturns["init"]

    for i in range(0, len(seq) - teloWindow, windowStep):
        if isGStrand:
            teloSeq = seq[(len(seq) - i) - teloWindow:len(seq) - i]
        else:
            teloSeq = seq[i:i + teloWindow]

        teloLen = len(teloSeq)
        teloSeqUpper = str(teloSeq.upper())

        currentOffsets = []

        for ntPatternEntry in composition:
            ntPattern = ntPatternEntry[0]
            patternComposition = ntPatternEntry[1]
            patternCount = len(re.findall(ntPattern, teloSeqUpper))

            if is_regex_pattern(ntPattern):
                if len(ntPatternEntry) != 3:
                    raise ValueError(
                        "Regex pattern requires target length as third list item."
                    )

                regexTargetLength = ntPatternEntry[2]
                rawOffsetValue = (
                    (patternCount * regexTargetLength) / teloLen
                    - patternComposition
                )
            else:
                rawOffsetValue = (
                    (patternCount * len(ntPattern)) / teloLen
                    - patternComposition
                )

            if patternComposition != 0:
                percentOffsetValue = (rawOffsetValue / patternComposition) * 100
            else:
                percentOffsetValue = rawOffsetValue * 100

            currentOffsets.append(percentOffsetValue)

        ntOffsets.append(currentOffsets)

    areaList = getGraphArea(ntOffsets, targetPatternIndex, graphAreaWindowSize)
    areaDiffs = np.diff(areaList)
    indexAtThreshold = -1

    if returnLastDiscontinuity:
        indexAtThreshold = next(
            (
                y for y in range(len(areaList) - 2, 0, -1)
                if areaList[y] > changeThreshold
                and (0 > (areaList[y + 1] - areaList[y]))
            ),
            indexAtThreshold,
        )

        if indexAtThreshold != -1:
            indexAtThreshold = next(
                (
                    y for y in range(indexAtThreshold, len(areaList) - 2)
                    if areaList[y] < plateauDetectionThreshold
                    and (0 > (areaList[y + 1] - areaList[y]))
                ),
                indexAtThreshold,
            )
        else:
            indexAtThreshold = next(
                (
                    y for y in range(len(areaList) - 2)
                    if areaList[y] < plateauDetectionThreshold
                    and (0 > (areaList[y + 1] - areaList[y]))
                ),
                indexAtThreshold,
            )
    else:
        indexAtThreshold = next(
            (
                y for y in range(len(areaList) - 2)
                if areaList[y] < plateauDetectionThreshold
                and (0 > (areaList[y + 1] - areaList[y]))
            ),
            indexAtThreshold,
        )

    if indexAtThreshold == -1:
        return errorReturns["init"]

    for x in range(indexAtThreshold, len(areaDiffs) - 1):
        if (
            abs(areaDiffs[x]) < areaDiffsThreshold
            or areaDiffs[x] < areaDiffsThreshold
            and areaDiffs[x + 1] > areaDiffsThreshold
        ):
            boundaryPoint = x * windowStep
            break

    if boundaryPoint == -1:
        boundaryPoint = len(areaDiffs) * windowStep

    if secondarySearch:
        telomereOffset = 500
        subTelomereOffset = 1000
        telomereOffsetRE = 30
        subTelomereOffsetRE = 30

        ntPattern = composition[targetPatternIndex][0]
        patternComposition = composition[targetPatternIndex][1]

        if isGStrand:
            lowerIndex = len(seq) - (boundaryPoint + subTelomereOffset)
            upperIndex = len(seq) - (boundaryPoint - telomereOffset)

            lowerIndex = max(lowerIndex, 0)
            upperIndex = min(upperIndex, len(seq))

            scanSeq = seq[lowerIndex:upperIndex]

            if len(scanSeq) >= 100:
                secBoundary = getTeloBoundary(
                    scanSeq,
                    isGStrand,
                    compositionGStrand=composition,
                    compositionCStrand=composition,
                    teloWindow=90,
                    windowStep=6,
                    changeThreshold=changeThreshold,
                    plateauDetectionThreshold=-60,
                    targetPatternIndex=-1,
                    nucleotideGraphAreaWindowSize=100,
                    returnLastDiscontinuity=returnLastDiscontinuity,
                    secondarySearch=False,
                )

                if secBoundary >= 0:
                    tempBoundary = boundaryPoint + secBoundary - telomereOffset

                    if upperIndex == len(seq):
                        tempBoundary = secBoundary

                    lower = len(seq) - (tempBoundary + subTelomereOffsetRE)
                    upper = len(seq) - (tempBoundary - telomereOffsetRE)

                    scanSeq = seq[max(lower, 0):min(upper, len(seq))]
                    scan_pattern = "(" + ntPattern + ")" + "(" + ntPattern + ")"

                    match = re.search(str(scan_pattern), str(scanSeq))

                    if match:
                        secBoundary = len(scanSeq) - match.span()[0]
                        boundaryPoint = tempBoundary + secBoundary - telomereOffsetRE
                    else:
                        boundaryPoint = tempBoundary

        else:
            lowerIndex = boundaryPoint - telomereOffset
            upperIndex = boundaryPoint + subTelomereOffset

            lowerIndex = max(lowerIndex, 0)
            upperIndex = min(upperIndex, len(seq))

            scanSeq = seq[lowerIndex:upperIndex]

            if len(scanSeq) >= 100:
                secBoundary = getTeloBoundary(
                    scanSeq,
                    isGStrand,
                    compositionGStrand=composition,
                    compositionCStrand=composition,
                    teloWindow=90,
                    windowStep=6,
                    changeThreshold=changeThreshold,
                    plateauDetectionThreshold=-60,
                    targetPatternIndex=-1,
                    nucleotideGraphAreaWindowSize=100,
                    returnLastDiscontinuity=returnLastDiscontinuity,
                    secondarySearch=False,
                )

                if secBoundary >= 0:
                    tempBoundary = boundaryPoint + secBoundary - telomereOffset

                    if lowerIndex == 0:
                        tempBoundary = secBoundary

                    scanSeq = seq[
                        max(tempBoundary - telomereOffsetRE, 0):
                        min(tempBoundary + subTelomereOffsetRE, len(seq))
                    ]

                    if patternComposition == 1:
                        scan_pattern = "(" + ntPattern + ")" + "(" + ntPattern + ")"
                        matches = [match for match in re.finditer(scan_pattern, str(scanSeq))]

                        if matches:
                            teloEnd = matches[-1].end()
                            boundaryPoint = tempBoundary + teloEnd - telomereOffsetRE
                        else:
                            boundaryPoint = tempBoundary
                    else:
                        boundaryPoint = tempBoundary

    return boundaryPoint



def getTeloNPBoundary(
    seq,
    isGStrand=None,
    compositionCStrandIn=teloNPTeloCompositionCStrand,
    compositionGStrandIn=teloNPTeloCompositionGStrand,
    teloWindow=100,
    windowStep=6,
    changeThreshold=-20,
    plateauDetectionThreshold=-60,
    targetPatternIndex=-1,
    nucleotideGraphAreaWindowSize=750,
    returnLastDiscontinuity=True,
    secondarySearch=True,
):
    return getTeloBoundary(
        seq,
        isGStrand=isGStrand,
        compositionCStrand=compositionCStrandIn,
        compositionGStrand=compositionGStrandIn,
        teloWindow=teloWindow,
        windowStep=windowStep,
        changeThreshold=changeThreshold,
        plateauDetectionThreshold=plateauDetectionThreshold,
        targetPatternIndex=targetPatternIndex,
        nucleotideGraphAreaWindowSize=nucleotideGraphAreaWindowSize,
        returnLastDiscontinuity=returnLastDiscontinuity,
        secondarySearch=secondarySearch,
    )



def open_input_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def open_output_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")

def has_consecutive_telomeric_repeats(
    seq,
    min_repeats=6
):
    seq_upper = str(seq).upper()

    g_ok = re.search(
        rf"(TTAGGG){{{min_repeats},}}",
        seq_upper
    )

    c_ok = re.search(
        rf"(CCCTAA){{{min_repeats},}}",
        seq_upper
    )

    return g_ok is not None, c_ok is not None

def extract_fused_reads(input_fastq, output_fastq, output_ids=None):
    total = 0
    fused = 0
    init_error = 0
    strand_error = 0
    seq_not_found = 0
    positive = 0

    ids_handle = open(output_ids, "w") if output_ids else None

    with open_input_fastq(input_fastq) as handle, open_output_fastq(output_fastq) as out:
        for record in SeqIO.parse(handle, "fastq"):
            total += 1
            seq = str(record.seq)

            if seq is None or len(seq) == 0:
                result = errorReturns["seqNotFound"]
            else:
                result = getTeloNPBoundary(seq)

            if result == errorReturns["fusedRead"]:
                g_ok, c_ok = has_consecutive_telomeric_repeats(
                    seq,
                    min_repeats=6
                )
                if g_ok and c_ok:
                    fused += 1
                    SeqIO.write(record, out, "fastq")

                    if ids_handle is not None:
                        ids_handle.write(record.id + "\n")
            elif result == errorReturns["init"]:
                init_error += 1
            elif result == errorReturns["strandType"]:
                strand_error += 1
            elif result == errorReturns["seqNotFound"]:
                seq_not_found += 1
            elif result >= 0:
                positive += 1

    if ids_handle is not None:
        ids_handle.close()

    print("Done.")
    print(f"Input reads: {total}")
    print(f"Fused reads : {fused}")
    if output_ids:
        print(f"Output read IDs: {output_ids}")


def main():
    parser = argparse.ArgumentParser(
        description="Standalone script to extract teloNP fused reads (-10) from FASTQ/FASTQ.gz."
    )

    parser.add_argument(
        "input_fastq",
        help="Input FASTQ or FASTQ.gz",
    )

    parser.add_argument(
        "output_fastq",
        help="Output FASTQ or FASTQ.gz containing only teloNP fused reads (-10)",
    )

    parser.add_argument(
        "--output_ids",
        default=None,
        help="Optional output txt file for fused read IDs",
    )

    args = parser.parse_args()

    extract_fused_reads(
        input_fastq=args.input_fastq,
        output_fastq=args.output_fastq,
        output_ids=args.output_ids,
    )


if __name__ == "__main__":
    main()