#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Parse masking coordinates from TRF output.
# ==============================================================


import argparse
import csv
import re
import sys
from enum import Enum


class State(Enum):
    START = 1
    SEQUENCE_START = 2
    READING_REPEATS = 3


def extractFeatureGtf(text, feature):
    regex = feature + ' "([^"]+)"'
    return re.search(regex, text).groups()[0]


def extractFeatureGff(text, feature):
    regex = feature + '=([^;]+);'
    return re.search(regex, text).groups()[0]


def printStatistics(prefix, periodCounts, periodSequenceLengths, baseCounts):
    with open(prefix + ".periodCounts", "w") as f:
        f.write("period\tcount\n")
        for key, value in sorted(periodCounts.items(),
                                 key=lambda item: int(item[0])):
            f.write(key + "\t" + str(periodCounts[key]) + "\n")

    with open(prefix + ".periodSeqLengths", "w") as f:
        f.write("period\tSeqLen\n")
        for key, value in sorted(periodSequenceLengths.items(),
                                 key=lambda item: int(item[0])):
            f.write(key + "\t" + str(periodSequenceLengths[key]) + "\n")

    if len(baseCounts) == 0:
        return

    with open(prefix + ".GC", "w") as f:
        f.write("Repeat period\tGC (%)\n")

        for key, value in sorted(baseCounts.items(),
                                 key=lambda item: int(item[0])):
            gc = 0
            if baseCounts[key]["G"] + baseCounts[key]["C"] != 0:
                gc = ((baseCounts[key]["G"] + baseCounts[key]["C"]) /
                      (baseCounts[key]["G"] + baseCounts[key]["C"] +
                       baseCounts[key]["A"] + baseCounts[key]["T"]))
            f.write(key + "\t" + str(gc) + "\n")


def addToDict(dictionary, key, item):
    if key not in dictionary:
        dictionary[key] = item
    else:
        dictionary[key] += item


def initCounts():
    bases = "ACGT"
    counts = {}
    for base in bases:
        counts[base] = 0
    return counts


def addToGc(baseCounts, period, seq):

    if period not in baseCounts:
        baseCounts[period] = initCounts()

    for base in seq:
        if base == "N":
            continue
        baseCounts[period][base] += 1


def parse(trfOutput, minCopies, gc):
    state = State.START
    periodCounts = {}
    periodSequenceLengths = {}
    baseCounts = {}
    sequence = ""
    line = 1
    for row in csv.reader(open(trfOutput), delimiter=' '):
        line += 1
        if len(row) == 0:
            continue
        first = row[0]

        if state == State.START:
            if first == "Sequence:":
                state = State.SEQUENCE_START
                sequence = row[1]
        elif state == State.SEQUENCE_START:
            if first == "Parameters:":
                state = State.READING_REPEATS
            else:
                sys.exit("Error when parsing line " + str(line) + ": expected"
                         " a list of Parameters after sequence definition.")
        elif state == State.READING_REPEATS:
            if first == "Sequence:":
                state = State.SEQUENCE_START
                sequence = row[1]
            elif first.isnumeric():
                start = row[0]
                end = row[1]
                copies = float(row[3])
                period = row[2]
                if copies >= minCopies:
                    print("\t".join([sequence, "trf", "repeat", start, end,
                                     ".", ".", ".", "."]))
                    addToDict(periodCounts, period, 1)
                    addToDict(periodSequenceLengths, period,
                              int(end) - int(start) + 1)
                    if gc:
                        addToGc(baseCounts, period, row[14])
            else:
                sys.exit("Error when parsing line " + str(line) + ": unexpected"
                         "entry")
    return periodCounts, periodSequenceLengths, baseCounts


def main():
    args = parseCmd()
    periodCounts, periodSequenceLengths, baseCounts = parse(args.trfOutput,
                                                            args.minCopies,
                                                            args.gc)
    if args.statisticsPrefix:
        printStatistics(args.statisticsPrefix, periodCounts,
                        periodSequenceLengths, baseCounts)


def parseCmd():

    parser = argparse.ArgumentParser(description='Parse masking coordinates\
                                     from TRF output.')

    parser.add_argument('trfOutput', type=str,
                        help='Trf output .dat file')

    parser.add_argument('--minCopies', type=int, default=1,
                        help='Minimum required copies of a repeat to be\
                        printed')

    parser.add_argument('--statisticsPrefix', type=str,
                        help='If specified, save statistics about repeats\
                        into files with this prefix.')

    parser.add_argument('--gc',  default=False, action='store_true',
                        help='Compute GC content per period')


    return parser.parse_args()


if __name__ == '__main__':
    main()
