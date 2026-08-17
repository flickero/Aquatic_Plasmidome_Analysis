# -*- coding: utf-8 -*-
"""
Created on 21/09/2023

Author: Lucy Androsiuk
"""

import time, re, os, sys, subprocess, logging, shutil, dotenv, errno, glob, threading
from dotenv import load_dotenv
from pathlib import Path
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio import SearchIO
from Bio import SeqIO
import pandas as pd
from shutil import copyfile
from csv import writer
from subprocess import Popen,PIPE

load_dotenv()
### directions
## working directories
outing = r"../"
db_dir = r"../DBs"
resource = r"../resources"
logs = r"../logs"
Path(logs).mkdir(parents=True, exist_ok=True)
Path(db_dir).mkdir(parents=True, exist_ok=True)

SampleName = sys.argv[1]
print(SampleName)
reads = sys.argv[2]
print(reads)

### working filesv
blast_dir = os.getenv("BLAST")
spades_dir = f'{os.getenv("SPADES")}spades.py'
trim_file = r'./trim_files.csh'
megahit_bash = r'./megahit_assembly.csh'
scapp_bash = r"./scapp.csh"

all_data_log = f'{logs}/all_data_{SampleName}.log'
lib_file = f'{resource}/library_size.csv'
start_time = time.time()

formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')

def setup_logger(name, log_file, level=logging.DEBUG):
    """To setup as many loggers as you want"""
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger

def error_tracker():
    error_file = logs + "/Error.log"
    error_logger = setup_logger("error_logging", error_file)
    exc_type, exc_value, traceback = sys.exc_info()
    error_logger.exception(sys.exc_info())
    all_logger.exception(sys.exc_info())
    print("exception caught: %s" % exc_type, exc_value)

def CreateDirectory (name):
    '''Creating directory with sepcified name in the parent_dir'''
    parent_dir = outing
    path = os.path.join(parent_dir, name)
    try:
        if not os.path.exists(path):
            os.mkdir(path)
            all_logger.info("Directory '%s' created" % path)
    except:
        raise  # Re-raise the exception if it's a different error
    return path

all_out = CreateDirectory('Output')

all_logger = setup_logger('second_logger', all_data_log)
all_logger.info("Working with %s." % SampleName)
def CreateFolderDirectory (name):
    '''Creating folder in the directory'''
    parent_dir = CreateDirectory(name)
    pathFolder = os.path.join(parent_dir, SampleName)
    print(f"The directory for {SampleName} will be created.")
    if not os.path.exists(pathFolder):
        os.mkdir(pathFolder)
        all_logger.info("Directory '%s' created" % pathFolder)
    return pathFolder

def Trimming(path, r1, r2):
    '''Running Trim_Galore to trim and clean read files'''
    try:
        trim_time = time.time()
        pattern = "*fastq.gz"
        # Find all files matching the pattern
        fastq_files = glob.glob(os.path.join(path, pattern))

        if len(fastq_files) == 2:
            all_logger.info("Starting trimming and cleaning reads.")
            subprocess.call(["bash", trim_file, path, r1, r2, SampleName, logs])
            all_logger.warning("Check ERRORs file for records")
            all_logger.info("Reads are trimmed and cleaned.")
            all_logger.info("Trimming took %s seconds" % (time.time() - trim_time))
            assert os.path.isfile(trim_file), "trimming bash file is missing"
        else:
            all_logger.warning("Looks like reads are already trimmed and cleaned. Check ERRORs file for records")
        # Use glob to find files matching the pattern
        if len(glob.glob(os.path.join(path, "*.fq.gz"))) == 2:
            print(f"2 files ending with '*.fq.gz'. Removing raw read files")
            files_to_remove = glob.glob(os.path.join(path, "*.fq.gz"))

            # Iterate through the list of files and remove them
            if file_path not in files_to_remove:
                os.remove(file_path)
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def GetLibrary(path):
    for file in os.listdir(path):
        all_logger.info("***** Getting information for %s " % file)
        if re.match('.*\.fq.gz', file):
            fullpath = os.path.join(path, file)

            # Use seqkit stats to obtain library size
            try:
                seqkit_cmd = ["seqkit", "stats", fullpath]
                process = subprocess.Popen(seqkit_cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                out, err = process.communicate()

                if process.returncode == 0:
                    # Parse the output to get the library size
                    lines = out.decode('utf-8').strip().split('\n')
                    numseqs = int(lines[1].split()[3].replace(',', ""))
                    print("****** The library size is %d reads" % numseqs)
                    return numseqs
                else:
                    all_logger.warning("Error running seqkit stats. Check the following error message:")
                    all_logger.warning(err.decode('utf-8'))
                    return None
            except Exception as e:
                all_logger.warning("Error running seqkit stats: %s" % str(e))
                return None

def RemoveF(folder_path, important_files):
    files_in_dir = os.listdir(folder_path)
    for file in files_in_dir:  # loop to delete each unnecessary file in folder
        path_to_file = f'{folder_path}/{file}'
        if os.path.isfile(path_to_file) and file not in important_files:
            all_logger.info("Will not need: %s. Removing it" % file)
            os.remove(f'{folder_path}/{file}')
        elif os.path.isdir(path_to_file) and file not in important_files:
            all_logger.info("Will not need: %s. Removing it" % file)
            shutil.rmtree(f'{folder_path}/{file}')
        elif os.path.isfile(path_to_file) and file in important_files:
            all_logger.warning("Will need the file %s. Saved it in output direcoty" % file)
        else:
            all_logger.warning("No files to delete anymore")

def RunSpades():
    '''Running SPAdes to generate assembly_graphs'''
    spades_time = time.time()
    try:
        if not os.path.isfile(output_dir + "/assembly_graph.fastg"):
            subprocess.call(["python", spades_dir,
                             "-1", trim_fw,
                             "-2", trim_rv,
                             "-o", output_dir])
            all_logger.info("SPAdes running took %s seconds" % (time.time()-spades_time))
        else:
            all_logger.warning("SPAdes has already created 'assembly_graph.fastg'")
        files_in_dir = os.listdir(output_dir + '/')  # get list of files in the directory
        for file in files_in_dir:  # loop to delete each unnecessary file in folder
            important_files = ['scaffolds.fasta', 'SCAPP_res', 'spades.log', 'assembly_graph.fastg', 'metaplasmidSPAdes']
            path_to_file = f'{output_dir}/{file}'
            if os.path.isfile(path_to_file) and file not in important_files:
                all_logger.info("Will not need: %s. Removing it" % file)
                os.remove(f'{output_dir}/{file}')
            elif os.path.isdir(path_to_file) and file not in important_files:
                all_logger.info("Will not need: %s. Removing it" % file)
                shutil.rmtree(f'{output_dir}/{file}')
            elif os.path.isfile(path_to_file) and file in important_files:
                all_logger.warning("Will need the file %s. Saved it in output direcoty" % file)
            else:
                all_logger.warning("No files to delete anymore")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def RunMegaHit():
    ''' Running MegaHit to generate assembly_graphs'''
    megahit_path = os.path.join(output_dir, "megahit")
    try:
        megahit_time = time.time()
        if not os.path.isfile(output_dir + "/assembly_graph.fastg"):
            all_logger.info("Starting assembly with MegaHit.")
            subprocess.call(["bash", megahit_bash, trim_fw, trim_rv, megahit_path, output_dir, logs])
            all_logger.warning("Check ERRORs file for records")
            all_logger.info("MegaHit took %s seconds" % (time.time() - megahit_time))
        else:
            all_logger.warning("MegaHit has already created 'assembly_graph.fastg'")

        # REMOVING unnecessary files, make sure all necessary files are in important_files list
        important_files = ['log', 'final.contigs.fa', 'intermediate_contigs']
        RemoveF(megahit_path, important_files)

    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def SCAPP():
    ''' Running SCAPP to assemble plasmids'''
    SCAPP_path = os.path.join(output_dir, "SCAPP_res")
    in_file = output_dir + "/assembly_graph.fastg"
    if not os.path.exists(SCAPP_path):
        os.mkdir(SCAPP_path)
        all_logger.info("Directory ' %s ' created" % SCAPP_path)
    try:
        scapp_time = time.time()
        if not os.path.exists(SCAPP_path + "/assembly_graph.confident_cycs.fasta"):
            all_logger.info("Starting SCAPP.")
            subprocess.call(["bash", scapp_bash, in_file, SCAPP_path, trim_fw, trim_rv, logs])
            all_logger.warning("Check ERRORs file for records")
            all_logger.info("assembly_graph.confident_cycs.fasta file created")
            all_logger.info("SCAPP took %s seconds" % (time.time() - scapp_time))

        # REMOVING unnecessary files, make sure all necessary files are in important_files list
        important_files = ['assembly_graph.confident_cycs.fasta', 'logs', 'intermediate_files', 'intermediate_files/assembly_graph.cycs.fasta']
        RemoveF(SCAPP_path, important_files)

    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def MetaPlasmidSPAdes():
    ''' Running metaplasmidSPAdes to assemble plasmids'''
    MPspades_path = os.path.join(output_dir, "metaplasmidSPAdes")
    if not os.path.exists(MPspades_path):
        os.mkdir(MPspades_path)
        all_logger.info("Directory ' %s ' created" % MPspades_path)
    try:
        metaplasmid_time = time.time()
        if not os.path.isfile(output_dir + "/metaplasmidSPAdes/scaffolds.fasta"):
            subprocess.call(["python", spades_dir,
                             "--meta", "--plasmid",
                             "-1", trim_fw,
                             "-2", trim_rv,
                             "-o", MPspades_path])
            all_logger.info("MetaPlasmidSPAdes running took %s seconds" % (time.time()-metaplasmid_time))
            all_logger.warning("Something may have gone wrong here. Check file 'scaffolds.fasta'.")
        else:
            all_logger.warning("Looks like MetaPlasmidSPAdes has already found nodes. Check file 'scaffolds.fasta' and 'spades.log' for errors.")

        # REMOVING unnecessary files, make sure all necessary files are in important_files list
        important_files = ['scaffolds.fasta', 'spades.log']
        RemoveF(MPspades_path, important_files)

        # Define input and output file paths
        input_file = f"{MPspades_path}/scaffolds.fasta"
        output_file = f"{MPspades_path}/filt_scaffolds.fasta"
        min_length = 1000

        # Read the input fasta file
        sequences = list(SeqIO.parse(input_file, "fasta"))

        # Filter sequences by length
        filtered_sequences = [seq for seq in sequences if len(seq) >= min_length]

        # Write the filtered sequences to a new fasta file
        SeqIO.write(filtered_sequences, output_file, "fasta")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()

def alignNodes():
    '''Comparing candidates obtained from SCAPP with ones obtained from metaplasmidSPAdes with BLASTn'''
    query_file = f'{output_dir}/SCAPP_res/intermediate_files/assembly_graph.cycs.fasta'
    subject_file = f'{output_dir}/metaplasmidSPAdes/filt_scaffolds.fasta'
    newName = "aligned.tab"
    newFile_path = f'{output_dir}/{newName}'
    newFile = open(newFile_path, 'wt')
    try:
        unique_nodes_time = time.time()
        if not os.stat(query_file).st_size==0 and not os.stat(subject_file).st_size==0:
            blastn = blast_dir + "blastn"
            all_logger.info("Starting the comparison of SCAPP and MetaPlasmidSPAdes outputs")
            cline = NcbiblastnCommandline(blastn,
                                          query = query_file,
                                          subject = subject_file,
                                          evalue = 0.001,
                                          out = newFile_path,
                                          outfmt = '6 qseqid sseqid evalue length pident score '
                                                         'qcovs qstart qend sstart send qseq sseq')
            all_logger.info(cline)
            os.system(str(cline))
            stdout, stderr = cline()
            newFile.close()
            all_logger.info("Alignment finished and written to the file aligned.tab in your output directory")
            all_logger.info("Alignment of SCAPP to MPSPAdes output nodes took %s seconds" % (time.time() - unique_nodes_time))
        else:
            open(newFile_path, "a").close()
            all_logger.warning("SCAPP output and/or MetaPlasmidSPAdes output are empty. Could not align.")
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()
    return newFile_path

def MetaPlasmidOnly():
    ''' Getting list of candidates detected with metaplasmidSPAdes only'''
    toCompare = alignNodes()
    all_logger.info("Analyzing SCAPP and MetaPlasmidSPAdes output")
    try:
        if not os.stat(toCompare).st_size == 0:
            df = pd.read_table(toCompare, header = None)
            custom_columns = ['qseqid', 'sseqid', 'evalue', 'length', 'pident',
                              'score', 'qcovs', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq']
            df.columns = custom_columns  # assigning specific columns
            df_filtered = df[(df['pident'] >= 99.0) & (df['qcovs'] >= 99.0)]  # filtering similar
            with open(output_dir + "/metaplasmidSPAdes/filt_scaffolds.fasta") as metaPlasmid_file:
                metaPlasmid_list = []  # list for unique metaplasmids
                for record in SeqIO.parse(metaPlasmid_file, "fasta"):
                    metaPlasmid_list.append(record.id)
                RecordToAppend = df_filtered['sseqid'].values.tolist()
                metaPlasmid_dif = [rec for rec in metaPlasmid_list if
                                   rec not in RecordToAppend]  # writing unique plasmids
            all_logger.info("Unique MetaPlasmidSPAdes nodes defined")
            all_logger.info(metaPlasmid_dif)
        elif not os.path.isfile(toCompare):
            all_logger.warning("There is no alignment. SCAPP or MetaPlasmidSPAdes output is empty.")
            with open(output_dir + "/metaplasmidSPAdes/filt_scaffolds.fasta") as metaPlasmid_file:
                metaPlasmid_list = []  # list for unique metaplasmids
                for record in SeqIO.parse(metaPlasmid_file, "fasta"):
                    metaPlasmid_list.append(record.id)
                metaPlasmid_dif = [rec for rec in metaPlasmid_list]
        else:
            all_logger.warning("There is no alignment. All records in SCAPP and MetaPlasmidSPAdes are unique.")
            with open(output_dir + "/metaplasmidSPAdes/filt_scaffolds.fasta") as metaPlasmid_file:
                metaPlasmid_list = []  # list for unique metaplasmids
                for record in SeqIO.parse(metaPlasmid_file, "fasta"):
                    metaPlasmid_list.append(record.id)
                metaPlasmid_dif = [rec for rec in metaPlasmid_list]
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()
    return metaPlasmid_dif

def ReplaceMetaPlasmidID():
    metaplasmidSpades = output_dir + "/metaplasmidSPAdes/filt_scaffolds.fasta"
    new_scaffold = output_dir + "/NewScaffolds.fasta"
    try:
        with open(new_scaffold, "wt") as new_spades:
            records = SeqIO.parse(metaplasmidSpades, 'fasta')
            for record in records:
                all_logger.info("We are replacing MetaPlasmidSPAdes record IDs./n Initial record ID is: %s" % record.id)
                prefix = str(SampleName[-6:]) + "_N"
                record.description = record.description.replace("N", prefix)
                record.id = record.description
                all_logger.info("We are assigning new ID to MetaPlasmidSPAdes record./n New record ID is: %s" % record.id)
                SeqIO.write(record, new_spades, 'fasta')
    except:
        all_logger.warning("Something has gone wrong here")
        error_tracker()
    return new_scaffold

# Create a lock
file_lock = threading.Lock()

def AppendToCombined():
    ''' Generating combined fasta for candidates'''
    try:
        scapp = f'{output_dir}/SCAPP_res/intermediate_files/assembly_graph.cycs.fasta'
        combined_dir = all_out + "/CombinedOutput.fasta"

        with open(scapp) as scapp_fasta:
            with file_lock:
                with open(combined_dir, 'a') as combined_fasta:
                    for header in scapp_fasta:
                        prefix = ">" + str(SampleName[-6:]) + "_"
                        combined_fasta.write(header.replace(">", prefix))

        listToCompare = MetaPlasmidOnly()
        spades_records = SeqIO.parse(ReplaceMetaPlasmidID(), 'fasta')

        for record in spades_records:
            name_regex = re.search("N\w+.+", record.id)
            if name_regex is not None:
                name = name_regex.group(0)
                all_logger.info("Appending '%s ' to list." % name)

                if name in listToCompare:
                    with file_lock:
                        all_logger.info("Writing '%s ' to CombinedOutput." % name)
                        combined_fasta = open(combined_dir, 'a')  # Reopen file within the lock
                        SeqIO.write(record, combined_fasta, 'fasta')
                        combined_fasta.close()
                else:
                    all_logger.warning("%s is duplicate" % name)

    except Exception as e:
        all_logger.warning("Something has gone wrong here: %s" % str(e))
        error_tracker()

if not os.path.isfile(lib_file):
    header = ['Sample', 'Library_size']
    with open(lib_file, 'a') as f_object:
        writer_object = writer(f_object)
        writer_object.writerow(header)
        f_object.close()

def get_file_names(folder_path):
    try:
        # Get a list of all files in the folder
        files = os.listdir(folder_path)

        # Filter out non-files (directories, etc.)
        file_names = [file for file in files if os.path.isfile(os.path.join(folder_path, file))]

        return file_names
    except OSError as e:
        print(f"Error reading files from {folder_path}: {e}")
        return []

all_logger.info("Starting processing of TaraOceans sample number: " + SampleName)

pattern = "*.fq.gz"
# Find all files matching the pattern
fq_files_before = glob.glob(os.path.join(reads, pattern))

if len(fq_files_before) != 2:
    raw_reads = "*fastq.gz"
    # Find all files matching the pattern
    raw_read_files = glob.glob(os.path.join(reads, raw_reads))

    if len(raw_read_files) != 2:
        all_logger.warning(f"Error: Expected 2 files ending with '*fastq.gz' but found {len(raw_read_files)}")
    else:
        file_names = get_file_names(reads)
        print(file_names)
        forward = [file for file in file_names if file.endswith("_1.fastq.gz")]
        print(forward)
        reverse = [file for file in file_names if file.endswith("_2.fastq.gz")]
        print(reverse)
        # Check if exactly one file matches the "_1.fastq.gz" pattern and one matches the "_2.fastq.gz" pattern
        if len(forward) == 1 and len(reverse) == 1:
            forward = forward[0]
            reverse = reverse[0]
            file_fw = os.path.join(reads, forward)
            file_rv = os.path.join(reads, reverse)
            print("Forward: %s" % file_fw)
            print("Reverse: %s" % file_rv)

            # Call Trimming function
            Trimming(reads, file_fw, file_rv)

            # Check if two files with the pattern "*norm.fq" exist after Trimming
            fq_files_after = glob.glob(os.path.join(reads, pattern))
            if len(fq_files_after) == 2:
                # Check if two files with the pattern "*norm.fq" exist after Trimming
                f_nrom = [file for file in fq_files_after if file.endswith("_1_val_1.fq.gz")]
                r_norm = [file for file in fq_files_after if file.endswith("_2_val_2.fq.gz")]
                if len(f_nrom) == 1 and len(r_norm) == 1:
                    all_logger.info("Trimming function worked.")
                    trim_fw = f_nrom[0]
                    trim_rv = r_norm[0]
                    all_logger.info("Forward: %s" % trim_fw)
                    all_logger.info("Reverse: %s" % trim_rv)
                else:
                    all_logger.warning("Error: Trimming did not produce the expected output.")
        else:
            all_logger.warning("Error: Incorrect file naming pattern. Check for '_1.fastq.gz' and '_2.fastq.gz' suffixes.")
else:
    # Assuming you want to store the file paths in variables trim_fw and trim_rv
    #trim_fw, trim_rv = fq_files_before
    f_nrom = [file for file in fq_files_before if file.endswith("_1_val_1.fq.gz")]
    r_norm = [file for file in fq_files_before if file.endswith("_2_val_2.fq.gz")]
    trim_fw = f_nrom[0]
    trim_rv = r_norm[0]
    all_logger.info("Forward: %s" % trim_fw)
    all_logger.info("Reverse: %s" % trim_rv)
    
lib_size = GetLibrary(reads)
if lib_size is not None:
    # Continue with the rest of your code
    print("Library size:", lib_size)
    for_csv = [SampleName, lib_size]
    with open(lib_file, 'a') as f_object:
        writer_object = writer(f_object)
        writer_object.writerow(for_csv)
        f_object.close()
else:
    print("Failed to obtain library size.")


output_dir = CreateFolderDirectory('Output')
all_logger.info("***************WE START MegaHit for sample #%s****************" % SampleName)
RunMegaHit()
all_logger.info("***************MegaHit for sample #%s finished****************" % SampleName)
all_logger.info("***************WE START SCAPP for sample #%s****************" % SampleName)
SCAPP()
all_logger.info("***************SCAPP for sample #%s finished****************" % SampleName)
all_logger.info("***************WE START MetaPlasmidSPAdes for sample #%s****************" % SampleName)
MetaPlasmidSPAdes()
all_logger.info("***************MetaPlasmidSPAdes finished for sample #%s****************" % SampleName)
shutil.rmtree(reads)
all_logger.info("Input directory removed for sample #%s" % SampleName)
all_logger.info("***************APPENDING RECORDS for sample #%s TO COMBINED FILE****************" % SampleName)
AppendToCombined()
all_logger.info("***************FINISHED APPENDING RECORDS for sample #%s TO COMBINED FILE****************" % SampleName)

all_logger.info("--- %s seconds ---" % (time.time() - start_time))

