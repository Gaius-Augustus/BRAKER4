#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
from inspect import currentframe, getframeinfo
import shutil
import re
from datetime import datetime
from collections import defaultdict

__author__ = "Katharina J. Hoff"
__copyright__ = "Copyright 2023. All rights reserved."
__credits__ = "Huang Neng"
__license__ = "Artistic License 2.0, in part Apache License (see notes in code about functions copied & modified from compleasm)"
__version__ = "1.0.1"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "production"

argparser = argparse.ArgumentParser(description = 'Find or build the best gene set generated with BRAKER and ' + 
                                    'TSEBRA minimizing missing BUSCOs using compleasm. Will only compute a ' +
                                    'new BRAKER gene set if the percentage of missing BUSCOs exceeds 5%.')
argparser.add_argument('-m', '--tmp_dir', type=str, required = True,
                          help = 'Temporary directory where intermediate files will be written')
argparser.add_argument('-c', '--compleasm_bin', type=str, required = False,
                          help = 'Location of compleasm.py on system')
argparser.add_argument('-d', '--input_dir', type=str, required = True,
                            help = 'Output directory of BRAKER')
argparser.add_argument('-y', '--tsebra', type=str, required = False,
                            help = 'Location of tesbra.py on system')
argparser.add_argument('-f', '--getanno', type=str, required = False,
                            help = 'Location of getAnnoFastaFromJoingenes.py on system')
argparser.add_argument('-g', '--genome', type=str, required = True,
                            help = 'Genome FASTA file')
argparser.add_argument('-t', '--threads', type=str, required = False, default = 1,
                       help = 'Number of threads to use for running compleasm')
argparser.add_argument('-p', '--busco_db', type=str, required = True,
                            help = 'BUSCO lineage for running compleasm')
argparser.add_argument('-n', '--missing_busco_threshold', type=int, default = 20,
                       required = False, help = "Threshold for the percentage of missing BUSCOs that " +
                       "decides until which point the BUSCOs in augustus and genemark gene sets will " +
                       "simply be added on top of the braker.gtf gene set, upper boundary.")
argparser.add_argument('-a', '--diff_to_braker', type=int, default = 5,
                        required = False, help = "Difference in percentage of missing BUSCOs between " +
                        "braker.gtf and the other gene sets that is allowed for the BUSCOs in the " +
                        "augustus and genemark gene sets to be added on top of the braker.gtf gene set.")

argparser.add_argument('-L', '--library_path', type=str, required = False,
                       help = 'Path to pre-downloaded BUSCO lineage data (compleasm mb_downloads directory). '
                              'If provided, compleasm will use this instead of downloading.')
argparser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))

args = argparser.parse_args()

# Functions copied from compleasm by Huang Neng under the Apache License 2.0
# modifications are noted in the comments

def load_score_cutoff(scores_cutoff_file): # changed arguments
    cutoff_dict = {}
    try:
        with open(scores_cutoff_file, "r") as f:
            for line in f:
                line = line.strip().split()
                try:
                    taxid = line[0]
                    score = float(line[1])
                    cutoff_dict[taxid] = score
                except IndexError:
                    raise Error("Error parsing the scores_cutoff file.")
    except IOError:
        raise Error("Impossible to read the scores in {}".format(scores_cutoff_file))
    return cutoff_dict

def load_length_cutoff(lengths_cutoff_file): # changed arguments
    cutoff_dict = {}
    try:
        with open(lengths_cutoff_file, "r") as f:
            for line in f:
                line = line.strip().split()
                try:
                    taxid = line[0]
                    sigma = float(line[2])
                    length = float(line[3])
                    if sigma == 0.0:
                        sigma = 1
                    cutoff_dict[taxid] = {}
                    cutoff_dict[taxid]["sigma"] = sigma
                    cutoff_dict[taxid]["length"] = length
                except IndexError:
                    raise Error("Error parsing the lengths_cutoff file.")
    except IOError:
        raise Error("Impossible to read the lengths in {}".format(lengths_cutoff_file))
    return cutoff_dict

def parse_hmmsearch_output(hmm_out_dir, score_cutoff_dict, length_cutoff_dict): # turned this into a function
    busco_tx = {}
    for hmmsearch_output in os.listdir(hmm_out_dir):
        outfile = os.path.join(hmm_out_dir, hmmsearch_output)
        try:
            if outfile.endswith(".out"): # changed that summary.txt is ignored
                with open(outfile, 'r') as fin:
                    coords_dict = defaultdict(list)
                    for line in fin:
                        if line.startswith('#'):
                            continue
                        line = line.strip().split()
                        target_name = line[0]
                        query_name = line[3]
                        hmm_score = float(line[7])
                        hmm_from = int(line[15])
                        hmm_to = int(line[16])
                        assert hmm_to >= hmm_from
                        if hmm_score < score_cutoff_dict[query_name]:
                            # failed to pass the score cutoff
                            continue
                        coords_dict[target_name].append((hmm_from, hmm_to))
                    for tname in coords_dict.keys():
                        coords = coords_dict[tname]
                        interval = []
                        coords = sorted(coords, key=lambda x: x[0])
                        for i in range(len(coords)):
                            hmm_from, hmm_to = coords[i]
                            if i == 0:
                                interval.extend([hmm_from, hmm_to, hmm_to - hmm_from])
                            else:
                                try:
                                    assert hmm_from >= interval[0]
                                except:
                                    raise Error("Error parsing the hmmsearch output file {}.".format(outfile))
                                if hmm_from >= interval[1]:
                                    interval[1] = hmm_to
                                    interval[2] += hmm_to - hmm_from
                                elif hmm_from < interval[1] <= hmm_to:
                                    interval[2] += hmm_to - interval[1]
                                    interval[1] = hmm_to
                                elif hmm_to < interval[1]:
                                    continue
                                else:
                                    raise Error("Error parsing the hmmsearch output file {}.".format(outfile))
                        busco_tx[tname] = 0 # change that targets to with BUSCO hits are stored
        except IOError:
            print("Error opening file {}".format(outfile))
            sys.exit(1)
    return busco_tx

# end of functions under Apache License 2.0

def parse_compleasm(file):
    """
    Parse the compleasm statistics from the specified file.

    Args:
        file (str): Path to the file containing compleasm statistics.

    Returns:
        float: percentage of missing BUSCOs.

    Raises:
        SystemExit: If the file cannot be opened.

    """
    missing = 0
    stat_pattern = r'M:(\d+\.\d+)\%, \d+'
    try:
        with open(file, "r") as f:
            for line in f:
                if re.search(stat_pattern, line):
                    missing = float(re.search(stat_pattern, line).group(1))
    except IOError:
        print("ERROR: Could not open file: " + file)
        sys.exit(1)
    return missing


def find_input_files(args):
    """
    Find all the required input files for the process.

    Args:
        args (argparse.Namespace): Command-line arguments.

    Returns:
        dict: Dictionary containing the file paths.

    """
    file_paths = {}
    file_paths["braker_aa"] = check_file(os.path.join(args.input_dir, "braker.aa"))
    file_paths["augustus_aa"] = check_file(os.path.join(args.input_dir, "Augustus", "augustus.hints.aa"))
    if not file_paths["augustus_aa"]:
        file_paths["augustus_aa"] = check_file(os.path.join(args.input_dir, "augustus.hints.aa"))
    file_paths["genome"] = check_file(args.genome)  # required to generate genemark protein file
    file_paths["braker_gtf"] = check_file(os.path.join(args.input_dir, "braker.gtf"))
    file_paths["augustus_gtf"] = check_file(os.path.join(args.input_dir, "Augustus", "augustus.hints.gtf"))
    if not file_paths["augustus_gtf"]:
        file_paths["augustus_gtf"] = check_file(os.path.join(args.input_dir, "augustus.hints.gtf")) 
    return file_paths

def find_genemark_gtf(input_dir):
    """
    Find the genemark.gtf file and corresponding training.gtf file in the specified directory.

    Args:
        input_dir (str): Path to the BRAKER working directory.

    Returns:
        tuple: Tuple containing the paths to the genemark.gtf and training.gtf files.
        tuple consists of False,False if no files are found.

    """
    genemark_directories = ["GeneMark-ETP", "GeneMark-EP", "GeneMark-ET", "GeneMark-ES"]

    for genemark_dir in genemark_directories:
        genemark_gtf_file = os.path.join(input_dir, genemark_dir, "genemark.gtf")
        training_gtf_file = os.path.join(input_dir, genemark_dir, "training.gtf")
        if not os.path.isfile(training_gtf_file):
            training_gtf_file = os.path.join(input_dir, "traingenes.gtf")
        if os.path.isfile(genemark_gtf_file):
            return check_file(genemark_gtf_file), check_file(training_gtf_file)

    return False, False


def run_compleasm(protein_files, threads, busco_db, tmp_dir):
    """
    Run compleasm on the specified protein files

    Args:
        protein_files (list): List of protein files to run compleasm on.
        threads (int): Number of threads to use.
        busco_db (str): Path to the BUSCO database.
        tmp_dir (str): Temporary directory to store the output.

    Returns:
        dict: Dictionary containing the paths to the compleasm result files.

    Raises:
        SystemExit: If there is an error in execting compleasm.

    """
    # Determine BUSCO library directory.
    # If library_path is provided, use the pre-downloaded data instead of downloading.
    if args.library_path is not None:
        lineage_dir = os.path.join(args.library_path, args.busco_db)
        if not os.path.exists(lineage_dir):
            print("ERROR: Pre-downloaded BUSCO lineage not found at " + lineage_dir)
            sys.exit(1)
    else:
        lineage_dir = "mb_downloads/" + args.busco_db
        if not os.path.exists(lineage_dir):
            # Strip any _odbNN suffix before requesting download
            db_for_download = re.sub(r'_odb\d+$', '', busco_db)
            compleasm_cmd = [args.compleasm_bin, "download", db_for_download]
            run_simple_process(compleasm_cmd)

    # read key data of BUSCO lineage
    score_cutoff_dict = load_score_cutoff(os.path.join(lineage_dir, "scores_cutoff"))
    length_cutoff_dict = load_length_cutoff(os.path.join(lineage_dir, "lengths_cutoff"))

    protein_hmmsearch_output_dict = {}  ## key: hmm protein name, value: list of aligned hmm and complete or fragment
    for profile in os.listdir(os.path.join(lineage_dir, "hmms")):
        outfile = profile.replace(".hmm", ".out")
        target_specie = profile.replace(".hmm", "")
        protein_hmmsearch_output_dict[target_specie] = []

    # run compleasm on the protein files
    for protein_file in protein_files:
        # create a tool-specific output subdirectory
        tool = re.search(r'^([^.]+)\.', os.path.basename(protein_file)).group(1)
        tool_out_dir = args.tmp_dir + "/" + tool
        compleasm_cmd = [args.compleasm_bin, "protein", "-p", protein_file, "-l", busco_db, "-t", str(args.threads), "-o", tool_out_dir]
        if args.library_path is not None:
            compleasm_cmd.extend(["--library_path", args.library_path])
        run_simple_process(compleasm_cmd)

    result_dict = {}
    for protein_file in protein_files:
        # identify the gene prediction program from the protein file name
        tool = re.search(r'^([^.]+)\.', os.path.basename(protein_file)).group(1)
        result_dict[tool] = check_file(tmp_dir + "/" + tool + "/summary.txt")

    return result_dict, score_cutoff_dict, length_cutoff_dict, protein_hmmsearch_output_dict

def run_getanno(annobin, genome_file, gtf, output_dir):
    """
    Run the getAnnoFastaFromJoingenes.py tool to process the Genemark GTF file and generate protein sequences.

    Args:
        annobin: getAnnoFastaFromJoingenes.py script.
        genome_file (str): Path to the genome file.
        gtf (str): Path to the GTF file.
        output_dir (str): Output directory to store the generated protein sequences

    """
    tool = re.search(r'^([^.]+)\.', os.path.basename(gtf)).group(1)
    cmd = [annobin, '-g', genome_file, '-f', gtf, '-o', output_dir + "/" + tool, '-s', 'True']
    run_simple_process(cmd)
    return check_file(output_dir + "/" + tool + ".aa")

def run_simple_process(args_lst):
    """
    Execute a subprocess command with the provided arguments.

    Args:
        args_lst (list): List of command-line arguments.

    Returns:
        CompletedProcess: Result of the subprocess execution.

    Raises:
        SystemExit: If the subprocess returns a non-zero exit code.

    """
    try:
        print(" ".join(args_lst))
        result = subprocess.run(
            args_lst, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if(result.returncode == 0):
            return(result)
        else:
            frameinfo = getframeinfo(currentframe())
            print('Error in file ' + frameinfo.filename + ' at line ' +
                  str(frameinfo.lineno) + ': ' + "Return code of subprocess was " +
                  str(result.returncode) + str(result.args))
            sys.exit(1)
    except subprocess.CalledProcessError as grepexc:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed executing: ",
              " ".join(grepexc.args))
        print("Error code: ", grepexc.returncode, grepexc.output)
        sys.exit(1)


def check_binary(binary, name):
    """
    Check if the provided binary is executable or available in the system's PATH.

    Args:
        binary (str): Path to the binary file or None.
        name (str): Name of the binary/tool.

    Returns:
        str: Path to the binary if found and executable.

    Raises:
        SystemExit: If the binary is not found, not executable, or not in the PATH.

    """        
    # check if binary is in PATH
    if binary is not None:
        if os.path.isfile(binary):
            # check if binary is executable
            if os.access(binary, os.X_OK):
                return binary
            else:
                print("ERROR: " + name + " binary is not executable: " + binary)
                sys.exit(1)
        else:
            print("ERROR: " + name + " binary not found at " + binary)
            sys.exit(1)
    else:
        # check if binary is in PATH
        if shutil.which(name) is not None:
            return shutil.which(name)
        else:
            print("ERROR: " + name + " binary not found in PATH")
            sys.exit(1)


def check_file(file):
    """
    Check if the specified file exists.

    Args:
        file (str): Path to the file.

    Returns:
        str: Path to the file if it exists, otherwise False

    """
    if os.path.isfile(file):
        return file
    else:
        return False


def determine_mode(path_dir):
    """
    Determine whether BRAKER run, and whether the files are complete.

    Args:
        path_dir (dict): Dictionary containing the file paths.

    Returns:
        str: Mode of the run (BRAKER).

    Raises:
        SystemExit: If the files are not complete for either of the runs

    """
    if path_dir["braker_aa"] and path_dir["braker_gtf"] and path_dir["genome"] and path_dir["augustus_aa"] and path_dir["augustus_gtf"] and path_dir["genemark_gtf"]:
        return "BRAKER"
    else:
        print("ERROR: The specified directory does not contain all required files.")
        print("These are the files that were found:")
        print(path_dir)
        print("We require the following key files for a BRAKER run: braker.aa, braker.gtf, genome.fa, augustus.hints.aa, augustus.hints.gtf, genemark.gtf")
        sys.exit(1)


def check_dir(dir):
    """
    Check if the specified directory exists.

    Returns:
        str: Path to the directory if it exists.

    Raises:
        SystemExit: If the directory is not found.

    """
    if os.path.isdir(dir):
        return dir
    else:
        print("ERROR: Directory not found: " + dir)
        sys.exit(1)

def write_busco_keep_gtf(tx_dict, in_file, out_file):
    """
    Filters and writes GTF lines based on transcript IDs.

    This function reads a GTF file specified by 'in_file', filters its contents
    based on the presence of transcript IDs in 'tx_dict', and writes the 
    filtered content to 'out_file'. It retains lines starting with '#' 
    (comments) and writes only the lines where the transcript ID is found in
      'tx_dict'.

    Args:
        tx_dict (dict): A dictionary where keys are transcript IDs to be kept.
        in_file (str): Path to the input GTF file to be read.
        out_file (str): Path to the output file where filtered lines are 
                        written.

    Raises:
        IOError: If there is an error opening 'in_file' or 'out_file'.

    Notes:
        The function exits the script with a status of 1 if an IOError is 
        encountered.
        It assumes that each line in the GTF file has a 'transcript_id' 
        attribute.
    """
    try:
        with open(in_file, 'r') as fin, open(out_file, 'w') as fout:
            for line in fin:
                if line.startswith('#'):
                    continue
                if re.search(r'transcript_id \"([^"]+)\"', line):
                    tx_id = re.search(r'transcript_id \"([^"]+)\"', line).group(1)
                    if tx_id in tx_dict:
                        fout.write(line)
    except IOError:
        print("Error opening file " + in_file + " or " + out_file)
        sys.exit(1)

def main():
    """
    Execute workflow
    """
    # Step 0: complete the busco lineage name if necessary
    # Accept any _odbNN suffix; default to _odb12 if no suffix is given.
    if not re.search(r'_odb\d+$', args.busco_db):
        args.busco_db = args.busco_db + "_odb12"
    
    # Step 1: Find all input files
    file_paths = find_input_files(args)

    # Step 2: find out which GeneMark was used: ETP or EP or ET or ES or none
    file_paths["genemark_gtf"], file_paths["training_gtf"] = find_genemark_gtf(args.input_dir)

    # Step 3: determine whether BRAKER was run and whether files are complete
    w_mode = determine_mode(file_paths)

    # Step 4: Check if all dependencies are available
    args.tsebra = check_binary(args.tsebra, "tsebra.py")
    args.getanno = check_binary(args.getanno, "getAnnoFastaFromJoingenes.py")
    args.compleasm_bin = check_binary(args.compleasm_bin, "compleasm.py")

    # Step 5: Check if temporary directory exists
    # if not, create it
    if not os.path.exists(args.tmp_dir):
        try:
            os.makedirs(args.tmp_dir)
        except OSError:
            print("ERROR: Creation of the directory %s failed" % args.tmp_dir)
            sys.exit(1)

    if w_mode == "BRAKER":
        # Step 6: Create protein sequene file for GeneMark gtf file
        file_paths["genemark_aa"] = run_getanno(args.getanno, args.genome, file_paths["genemark_gtf"], args.tmp_dir)

    # Step 7: run compleasm
    if w_mode == "BRAKER":
        protein_file_list = [file_paths["genemark_aa"], file_paths["braker_aa"], file_paths["augustus_aa"]]

    compleasm_out_dict, score_dict, length_dict, protein_hmmsearch_dict = run_compleasm(protein_file_list, args.threads, args.busco_db, args.tmp_dir)

    # Step 8: parse and compare the number of missing BUSCOs
    if w_mode == "BRAKER":
        genemark_missing = parse_compleasm(compleasm_out_dict["genemark"])
        augustus_missing = parse_compleasm(compleasm_out_dict["augustus"])
        braker_missing = parse_compleasm(compleasm_out_dict["braker"])

    
    # Step 9: Decide whether the provided final gene set is good enough
    if w_mode == "BRAKER":
        print("BRAKER is missing " + str(braker_missing) + " BUSCOs.")
        print("GeneMark is missing " + str(genemark_missing) + " BUSCOs.")
        print("Augustus is missing " + str(augustus_missing) + " BUSCOs.")
        if braker_missing <= augustus_missing and braker_missing <= genemark_missing:
            print("The BRAKER gene set " + file_paths["braker_gtf"] + " is the best one. It lacks " + str(braker_missing) + "% BUSCOs.")
            sys.exit(0)
        elif (braker_missing <= args.missing_busco_threshold or ((braker_missing-augustus_missing)<args.diff_to_braker and (braker_missing-genemark_missing)<args.diff_to_braker)) and (augustus_missing < braker_missing or genemark_missing < braker_missing):
            print("All BUSCOs present in augustus.hints.gtf and genemark.gtf will be added to the braker.gtf gene set.")
            # in this case we want to find the BUSCOs in all gene sets and merge them on top of the braker gene set
            tsebra_force = file_paths["braker_gtf"]
        elif augustus_missing <= genemark_missing:
            tsebra_force = file_paths["augustus_gtf"]
            not_tsebra_force = file_paths["genemark_gtf"]
            print("Will enforce augustus.hints.gtf and the BUSCOs from augustus and genemark gene set.")
        elif genemark_missing < augustus_missing:
            print("Will enforce genemark.gtf and the BUSCOs from augustus and genemark gene set.")
            tsebra_force = file_paths["genemark_gtf"]
            not_tsebra_force = file_paths["augustus_gtf"]
        # used to have a case for GALBA but that was removed because we implemented DIAMOND filter for GALBA, instead

    # Step 10: find the BUSCOs in augustus and genemark gene sets
    hmmoutdir = os.path.join(args.tmp_dir, "augustus", args.busco_db + "_hmmsearch_output")
    keep_aug_dict = parse_hmmsearch_output(hmmoutdir, score_dict, length_dict)
    write_busco_keep_gtf(keep_aug_dict, file_paths["augustus_gtf"], args.tmp_dir + "/augustus_keep.gtf")
    # find the BUSCOs in genemark.gtf
    hmmoutdir = os.path.join(args.tmp_dir, "genemark", args.busco_db + "_hmmsearch_output")
    keep_gm_dict = parse_hmmsearch_output(hmmoutdir, score_dict, length_dict)
    write_busco_keep_gtf(keep_gm_dict, file_paths["genemark_gtf"], args.tmp_dir + "/genemark_keep.gtf")
    # concatenate the two keep gtf files
    try:
        with open(args.tmp_dir + "/augustus_keep.gtf", 'r') as fin, open(args.tmp_dir + "/genemark_keep.gtf", 'r') as fin2, open(args.tmp_dir + "/augustus_genemark_keep.gtf", 'w') as fout:
            for line in fin:
                fout.write(line)
            for line in fin2:
                fout.write(line)
    except IOError:
        print("Error opening file " + args.tmp_dir + "/augustus_keep.gtf" + " or " + args.tmp_dir + "/genemark_keep.gtf")
        sys.exit(1)

    # Step 11: Run TSEBRA and enforce the best gene set
    if w_mode == "BRAKER":
        if tsebra_force == file_paths["braker_gtf"]:
            # if augustus_genemark_keep.gtf is not empty
            if os.stat(args.tmp_dir + "/augustus_genemark_keep.gtf").st_size > 0:
                tsebra_cmd = [args.tsebra, "-k", tsebra_force + "," + args.tmp_dir + "/augustus_genemark_keep.gtf",
                              "-o", args.tmp_dir + "/better.gtf", "-q"]
            else:
                print("Attempted to merge additional BUSCOs onto braker.gtf but there are no BUSCOs to be added.")
                print("The BRAKER gene set " + file_paths["braker_gtf"] + " will be kept. It lacks " + str(braker_missing) + "% BUSCOs.")
                sys.exit(0)
        elif file_paths['training_gtf']:
            tsebra_cmd = [args.tsebra, "-k", file_paths["training_gtf"] + "," + tsebra_force # enforcing training.gtf may seem redundant if genemark is enforced but it causes no harm, and it must be enforced if augustus is enforced
                          + "," + args.tmp_dir + "/augustus_genemark_keep.gtf",
                          "-g", not_tsebra_force, "-o", args.tmp_dir + "/better.gtf", "-q"]
        else:
            print("Warning: No training.gtf file found. TSEBRA will be run without it, you may be using an older version of BRAKER.")
            tsebra_cmd = [args.tsebra, "-k", tsebra_force + "," + args.tmp_dir + "/augustus_genemark_keep.gtf",
                      "-g", not_tsebra_force, "-o", args.tmp_dir + "/better.gtf", "-q"]
    run_simple_process(tsebra_cmd)
    # name genes/transcripts consistently
    shutil.move(args.tmp_dir + "/better.gtf", args.tmp_dir + "/better_tmp.gtf")
    # replace tsebra.py by gtf_rename.py in args.tsebra
    rename_gtf = args.tsebra.replace("tsebra.py", "rename_gtf.py")
    gtf_rename_cmd = [rename_gtf, "--gtf", args.tmp_dir + "/better_tmp.gtf", "--out", args.tmp_dir + "/better.gtf"]
    run_simple_process(gtf_rename_cmd)
    os.remove(args.tmp_dir + "/better_tmp.gtf")

    # Step 12: generate protein sequence file for the new BRAKER gene set
    better_gtf = check_file(args.tmp_dir + "/better.gtf")
    bb_aa = run_getanno(args.getanno, args.genome, better_gtf, args.tmp_dir)

    # Step 13: Run compleasm on the new gene set
    secondary_compleasm_out, score_dict, length_dict, protein_hmmsearch_dict = run_compleasm([bb_aa], args.threads, args.busco_db, args.tmp_dir)

    # Step 14: parse compleasm output and report numbers
    better_missing = parse_compleasm(list(secondary_compleasm_out.values())[0])
    if w_mode == "BRAKER":
        if better_missing < braker_missing:
            print("The new best BRAKER gene set is " + better_gtf + ".")
            print("It is missing " + str(better_missing) + "% BUSCOs.")
        else:
            # if the new gene set is not superior, produce an output that tells the user what
            # of the previously existing gene sets had the lowest percentage of missing BUSCOs
            print("WARNING: The new BRAKER gene set is not better than the original one!")
            print("The best gene set produced by the original BRAKER run is:")
            gene_sets = [file_paths["braker_gtf"], file_paths["genemark_gtf"], file_paths["augustus_gtf"]]
            # Create a dictionary to map gene set names to their respective missing values
            gene_set_values = {
                file_paths["braker_gtf"]: braker_missing,
                file_paths["genemark_gtf"]: genemark_missing,
                file_paths["augustus_gtf"]: augustus_missing
            }

            # Find the minimum value among the three
            min_value = min(gene_set_values.values())

            # Check if braker.gtf has the smallest number of missing BUSCOs or is tied for the smallest number
            if gene_set_values[file_paths["braker_gtf"]] == min_value:
                min_gene_set = file_paths["braker_gtf"]
            else:
                # Find the gene set name associated with the minimum value (excluding braker.gtf)
                for gene_set, value in gene_set_values.items():
                    if gene_set != file_paths["braker_gtf"] and value == min_value:
                        min_gene_set = gene_set
                        break
            # Print the name of the gene set with the lowest number of missing BUSCOs
            print(min_gene_set)


if __name__ == "__main__":
    main()