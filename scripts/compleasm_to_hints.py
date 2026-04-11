#!/usr/bin/env python3

"""
Description: This script runs compleasm on a genome with a given BUSCO partition
             and parses the output table for complete BUSCOs to generate a hints
             file for AUGUSTUS/BRAKER/GALBA.
Author: Katharina J. Hoff
Email: katharina.hoff@uni-greifswald.de
Date: November 27th, 2023

Copyright (C) 2023, Katharina J. Hoff, University of Greifswald

This program is free software; you can redistribute it and/or modify
it under the terms of the Artistic License.
"""

# install instructions for sepp that will be needed for containers:
# git clone https://github.com/smirarab/sepp.git
# sudo apt-get install default-jre
# python3 setup.py config -c
# sudo python3 setup.py install
# other than that, the script depends on compleasm.py and its dependencies


import argparse
import re
import shutil
import os
import csv
import subprocess
from inspect import currentframe, getframeinfo

''' Function that runs a subprocess with arguments '''


def run_simple_process(args_lst):
    try:
        # bed files need sorting with LC_COLLATE=C
        myenv = os.environ.copy()
        myenv['LC_COLLATE'] = 'C'
        print("Trying to execute the following command:")
        print(" ".join(args_lst))
        result = subprocess.run(
            args_lst, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=myenv)
        print("Suceeded in executing command.")
        if(result.returncode == 0):
            return(result)
        else:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Return code of subprocess was " +
                  str(result.returncode) + str(result.args))
            quit(1)
    except subprocess.CalledProcessError as grepexc:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed executing: ",
              " ".join(grepexec.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        quit(1)


def extract_tx_ids_from_tsv(tsv_file):
    busco_ids = {}
    with open(tsv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            if (row['Status'] == 'Single' or row['Status'] == 'Duplicated') and (row['Frameshift events'] == '0'):
                busco_ids[row['Best gene']] = True
    return busco_ids

def read_and_filter_gff(gff_file, id_dict):
    gff_list = []
    try:
        with open(gff_file, 'r') as gff_handle:
            for line in gff_handle:
                if not line.startswith("#"):
                    if re.search(r'Target=', line):
                        bid = re.search(r'Target=(\S+)', line).group(1)
                        if bid in id_dict.keys():
                            gff_list.append(line)
    except IOError:
        print("Could not open file ", gff_file)
    return gff_list

def miniprot_to_hints(gff_lines):
    hints_lines = []
    for line in gff_lines:
        line = line.rstrip()
        line_fields = line.split("\t")
        if re.search(r'\tCDS\t.*Parent=([^;]+)', line):
            if (int(line_fields[4])-3) - (int(line_fields[3])+3) + 1 > 0:
                regex_result = re.search(r'\tCDS\t.*Parent=([^;]+)', line)
                hint = line_fields[0] + "\t" + "c2h" + "\tCDSpart\t" + str(int(line_fields[3])+3) + "\t" + str(int(line_fields[4])-3) + "\t" + "1" + "\t" + line_fields[6] + "\t" + line_fields[7] + "\t" + "src=M;grp=" + regex_result.group(1) + ";pri=4\n"
                hints_lines.append(hint)
    return hints_lines



def main():
    parser = argparse.ArgumentParser(description="Run compleasm and generate a gtf file with complete BUSCO genes.")

    # Mandatory input arguments
    parser.add_argument("-g", "--genome", required=True, help="Genome file in fasta format")
    parser.add_argument("-d", "--database", required=True, help="BUSCO database to use")
    parser.add_argument("-p", "--compleasm", required=False, help="Location of compleasm binary incl. binary name")
    parser.add_argument("-t", "--threads", required=False, help="Number of threads to use, default is 1")
    parser.add_argument("-o", "--output", required=True, help="Output file name, file is in GTF format.")
    parser.add_argument("-s", "--scratch_dir", required=False, help="Temporary directory for compleasm output, default is current directory, must be writable.")
    parser.add_argument("-L", "--library_path", required=False, help="Path to pre-downloaded BUSCO lineage data (compleasm mb_downloads directory). If provided, compleasm will use this instead of downloading.")
    args = parser.parse_args()

    if args.compleasm is None:
        # Try to find compleasm.py in PATH
        found = shutil.which('compleasm.py') or shutil.which('compleasm')
        if found:
            args.compleasm = found
        elif os.environ.get('COMPLEASM_PATH'):
            args.compleasm = str(os.environ.get('COMPLEASM_PATH')) + '/compleasm.py'
        elif os.path.exists('/opt/compleasm_kit/compleasm.py'):
            args.compleasm = '/opt/compleasm_kit/compleasm.py'
        else:
            raise FileNotFoundError("compleasm is not in PATH, COMPLEASM_PATH is not set, and /opt/compleasm_kit/compleasm.py does not exist")
    if args.compleasm is not None:
        # check whether provided compleasm is executable
        if not os.access(args.compleasm, os.X_OK):
            raise FileNotFoundError("compleasm is not executable")

    # Ensure database uses _odb12 suffix (compleasm only supports odb12).
    # Convert e.g. eukaryota_odb12 -> eukaryota_odb12, or add _odb12 if missing.
    if re.search(r'_odb\d+$', args.database):
        args.database = re.sub(r'_odb\d+$', '_odb12', args.database)
    else:
        args.database = args.database + "_odb12"

    # apply compleasm to genome file with run_subprocess and the database
    if args.scratch_dir is None:
        args.scratch_dir = "compleasm_genome_out"
    compleasm_cmd = [args.compleasm, 'run', '-l', args.database, '-a', args.genome, '-t', str(args.threads), '-o', args.scratch_dir]
    if args.library_path is not None:
        compleasm_cmd.extend(['--library_path', args.library_path])
    run_simple_process(compleasm_cmd)

    # Detect the actual lineage directory created by compleasm.
    # compleasm may auto-upgrade the lineage (e.g. eukaryota_odb12 -> eukaryota_odb12),
    # so we cannot assume the output directory matches the requested lineage name.
    actual_lineage_dir = os.path.join(args.scratch_dir, args.database)
    if not os.path.isdir(actual_lineage_dir):
        # Search for any _odb* directory in scratch_dir
        import glob
        lineage_base = re.sub(r'_odb\d+$', '', args.database)
        candidates = sorted(glob.glob(os.path.join(args.scratch_dir, lineage_base + '_odb*')))
        if candidates:
            actual_lineage_dir = candidates[-1]  # pick highest version
            print(f"[INFO] Lineage was upgraded: requested {args.database}, found {os.path.basename(actual_lineage_dir)}")
        else:
            # Fall back: any directory containing full_table.tsv
            for root, dirs, files in os.walk(args.scratch_dir):
                if 'full_table.tsv' in files:
                    actual_lineage_dir = root
                    print(f"[INFO] Found full_table.tsv in {root}")
                    break

    # parse compleasm output table for complete BUSCOs without frameshifts
    busco_ids = extract_tx_ids_from_tsv(actual_lineage_dir + '/full_table.tsv')

    # filter the miniprot alignments for those that have no frame shifts
    gff_lines = read_and_filter_gff(actual_lineage_dir + '/miniprot_output.gff', busco_ids)

    # convert the miniprot lines to CDSpart hints
    hints_lines = miniprot_to_hints(gff_lines)

    # gff_lines is not compatible with getAnnoFastaFromJoingenes, we need to fix this!
    try:
        with open(args.output, "w") as out_handle:
            for line in hints_lines:
                out_handle.write(line)
    except IOError:
        print("Failed to open file", args.out)

    # print the BUSCO scores for genome level statistics to STDOUT
    print("The following BUSCOs were found in the genome:")
    # compleasm puts summary.txt either at scratch_dir/summary.txt or
    # scratch_dir/{lineage}/summary.txt depending on version
    summary_found = False
    for summary_candidate in [
        args.scratch_dir + '/summary.txt',
        args.scratch_dir + '/' + args.database + '/summary.txt',
    ]:
        if os.path.exists(summary_candidate):
            try:
                with open(summary_candidate, 'r') as busco_summary:
                    for line in busco_summary:
                        line = line.strip()
                        print(line)
                summary_found = True
                break
            except IOError:
                pass
    if not summary_found:
        # Try to find any summary.txt in the scratch dir
        import glob
        summaries = glob.glob(args.scratch_dir + '/**/summary.txt', recursive=True)
        if summaries:
            with open(summaries[0], 'r') as busco_summary:
                for line in busco_summary:
                    print(line.strip())
        else:
            print("No summary.txt found in", args.scratch_dir)

    # Note: scratch_dir is NOT deleted here — Snakemake's collect_results
    # rule handles cleanup. The summary.txt is needed as a rule output.

if __name__ == "__main__":
    main()
