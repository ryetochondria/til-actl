#!/usr/bin/env python3
#SBATCH --job-name=array_job_test   # Job name
#SBATCH --ntasks-per-node=1         # Run a single task
#SBATCH --mem=120gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12          # Job Memory
#SBATCH --time=00:30:00             # Time limit hrs:min:sec
#SBATCH --output=star_stdout/array_%A-%a.log    # Standard output and error log
#SBATCH --array=1-50                # Array range

# Author: Francois Aguet
# Modified by Ryan Sowell 10/28/2018
# Ruddle FileName: run_STAR.slurm

from datetime import datetime
import argparse
import os
import sys
import subprocess
import gzip
import shutil
import contextlib
import multiprocessing
import re
sys.path.append(os.getcwd())

start_wallclock = datetime.datetime.now()

##### --- DEFINE FUNCTIONS --------------------------------------------------------------------------------------------------------- #####


@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

# List of samples in this experiment
# Import list from samplesheet



## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
#   Directory tree created by mk_dir.sh   #
# --------------------------------------- #
#    [parent directory for project]       #
#     |_ /data                            #
#         |                               #
#         |__ /pipelines                  #
#         |    |_______ /GTExTOPMed       #
#         |    |         |___ /genomes    #
#         |    |         |___ /scripts    #
#         |    |                          #
#         |    |_______ /TCR-MiXCR        #
#         |                               #
#         |__ /RNAseq                     #
#         |    |____ /fastqs_raw          #
#         |    |____ /outs                #
#         |                               #
#         |__ /TCRseq                     #
#              |____ /fastqs_raw          #
#              |____ /outs                #
#                                         #
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

##### =========================================================================================================================== #####
##### --- PIPELINE INPUTS ------------------------------------------------------------------------------------------------------- #####
##### =========================================================================================================================== #####

##### --- GLOBAL VARIABLES ------------------------------------------------------------------------------------------------------ #####
project_name = 'MedImmuneSRA'                                         # Project name and directory
data_type = 'RNAseq'                                                  # e.g. RNAseq, TCRseq, etc
pipeline_name = 'GTExTOPMed'                                          # Alignment and analysis pipeline to use 
                                                                      #  (should be contained in a folder w/ sub-directories genomes 
                                                                      #  and scripts)
# run STAR modifiers
testing=True
threads=12
twopassMode = 'None'


samplesheet_path = '/home/rts38/project/MedImmuneSRA/data/RNAseq/fastqs_raw'
samplesheet = open(samplesheet_path + '/samples.txt', "r")
samples = samplesheet.read().splitlines()

# Hard input of sample names
# samples = ['MEDIR201', 'MEDIR202', 'MEDIR205', 'MEDIR206', 'MEDIR207',
#            'MEDIR208', 'MEDIR209', 'MEDIR210', 'MEDIR213', 'MEDIR214',
#            'MEDIR301', 'MEDIR302', 'MEDIR303', 'MEDIR304', 'MEDIR305',
#            'MEDIR306', 'MEDIR307', 'MEDIR308', 'MEDIR309', 'MEDIR310',
#            'MEDIR311', 'MEDIR313', 'MEDIR314', 'MEDIR316', 'MEDIR317',
#            'MEDIR401', 'MEDIR402', 'MEDIR403', 'MEDIR404', 'MEDIR405',
#            'MEDIR406', 'MEDIR407', 'MEDIR408', 'MEDIR409', 'MEDIR410',
#            'MEDIR411', 'MEDIR412', 'MEDIR501', 'MEDIR502', 'MEDIR503',
#            'MEDIR504', 'MEDIR505', 'MEDIR506', 'MEDIR507', 'MEDIR508',
#            'MEDIR509', 'MEDIR510']

##### --- PATHS ----------------------------------------------------------------------------------------------------------------- #####
base_dir = os.path.join('/home/rts38/project/' + project_name + '/data')                      # Path to base directory
project_path = os.path.join(base_dir + '/' + data_type)                             # Path to parent folder of project
path_to_data_raw = os.path.join(project_path + '/fastqs_raw')                       # Path to raw_fastq files (in sample-specific directories)
star_outs = os.path.join(project_path + '/outs/STAR_outs')                          # Path to STAR output directory
rsem_outs = os.path.join(project_path + '/outs/RSEM_outs')                          # Path to RSEM output directory
pipeline = os.path.join(base_dir + '/pipelines/' + pipeline_name + '/scripts')      # Path to pipeline infrastructure
genomes = os.path.join(base_dir + '/pipelines/' + pipeline_name + '/genomes')       # Path to genome and annotation references directory


# def make_samplesheet(fastq_directory):
#     ss = 'bash' + os.path.join(project_path + '/runs/mk_samplesheet.sh')
#     subprocess.run(ss, shell=True, check=True)


# def read_samplesheet(samplesheet_path):
#     samplesheet = open(samplesheet_path + '/samples.txt', "r")
#     samples = samplesheet.read().splitlines()


##### --- PATH TO REFERENCE FILES ----------------------------------------------------------------------------------------------- #####
##### |______ Reference files generated using the GENCODE v26 build (GRCh38)


### STAR Index (ERCC genome included; HLA genes removed; built in STAR for 2x75bp)
star_index = os.path.join(genomes + '/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v26_oh75')
### RSEM Reference
rsem_reference = os.path.join(genomes + '/rsem_reference_GRCh38_gencode26_ercc')                  
### Genome Reference [GRCh38, GENCODE v26]  (ERCC genome included; HLA genes removed)
genome_fasta = os.path.join(genomes + '/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta')   
### Gene Annotations (ERCC genes included)
genes_gtf = os.path.join(genomes + '/gencode.v26.GRCh38.ERCC.genes.gtf')                          


##### --- SAMPLE INPUTS --------------------------------------------------------------------------------------------------------- #####
#### Load in list of samples
# make_samplesheet(path_to_data_raw)
# read_samplesheet(project_path + '/samples.txt')
# print('List of samples:')
# print(samples)

##### =========================================================================================================================== #####
##### --- EXECUTE STAR RUN ------------------------------------------------------------------------------------------------------ #####
##### =========================================================================================================================== #####
cd(os.path.join(project_path + '/outs'))

if not os.path.exists('star-stdout'):
    os.makedirs('star-stdout')
if not os.path.exists(star_outs):
    os.makedirs(star_outs)
if not os.path.exists(rsem_outs):
     os.makedirs(rsem_outs)

# Load STAR index into shared memory
loadSTARindex = 'STAR --genomeDir ' +star_index+ ' --genomeLoad LoadAndExit'
SLURM_ARRAY_TASK_ID = int(os.environ('SLURM_ARRAY_TASK_ID'))
if SLURM_ARRAY_TASK_ID == 0:
    subprocess.check_call(loadSTARindex, shell=True, executable='/bin/bash')
    print('Loading STAR index genome...')

# Variable assignment for current array task ID
curr_sample=str(samples[SLURM_ARRAY_TASK_ID])
#curr_sample = subprocess.check_call('`'+'sed "${SLURM_ARRAY_TASK_ID}q;d" '+samplesheet+'`', shell=True, executable=/bin/bash)
curr_path = os.path.join(path_to_data_raw, curr_sample)
curr_files = os.listdir(curr_path)
curr_outdir = os.path.join(star_outs, curr_sample)
curr_prefix = os.path.join(curr_outdir + '/' + curr_sample + '.')

# Create comma-separated list of fastq files for pair-end reads of current sample in array job
R1_fastq = ','.join(sorted(' '.join([os.path.join(curr_path, f) for f in curr_files if re.match('.*_R1_.*\.fastq\.gz', f)]).split()))
R2_fastq = ','.join(sorted(' '.join([os.path.join(curr_path, f) for f in curr_files if re.match('.*_R2_.*\.fastq\.gz', f)]).split()))
     
# Base STAR command line inputs
base_cmd = 'STAR --runMode alignReads --runThreadN ' + threads \
    + ' --genomeDir ' + star_index + ' --genomeLoad LoadAndKeep' \
    + ' --readFilesIn ' + R1_fastq + ' ' + R2_fastq \
    + ' --sjdbGTFfile ' + genes_gtf \
    + ' --twopassMode ' + twopassMode \
    + ' --readFilesCommand zcat' \
    + ' --outFileNamePrefix ' + curr_prefix #\
#    + ' --parametersFiles STAR.config'

# Additional variables for STAR 
##  move this eventually to a .config file to be called by this script
default_cmd += ' --outFilterMultimapNmax 20' \
    + ' --alignSJoverhangMin 8' \
    + ' --alignSJDBoverhangMin 1' \
    + ' --outFilterMismatchNmax 999' \
    + ' --outFilterMismatchNoverLmax 0.1' \
    + ' --alignIntronMin 20' \
    + ' --alignIntronMax 1000000' \
    + ' --alignMatesGapMax 1000000' \
    + ' --outFilterType BySJout' \
    + ' --outFilterScoreMinOverLread 0.33' \
    + ' --outFilterMatchNminOverLread 0.33' \
    + ' --limitSjdbInsertNsj 1200000' \
    + ' --outSAMstrandField intronMotif' \
    + ' --outFilterIntronMotifs None' \
    + ' --alignSoftClipAtReferenceEnds Yes' \
    + ' --quantMode TranscriptomeSAM GeneCounts' \
    + ' --outSAMtype BAM Unsorted' \
    + ' --outSAMunmapped Within' \
    + ' --chimSegmentMin 15' \
    + ' --chimJunctionOverhangMin 15' \
    + ' --chimOutType WithinBAM SoftClip' \
    + ' --chimMainSegmentMultNmax 1' \
    + ' --outSAMattributes NH HI AS nM NM ch' \
    + ' --outSAMattrRGline ID:rg1 SM: sm1'#\
#    + ' --sjdbFileChrStartEnd None'

print('--------------------------STAR RUN PARAMETERS-------------------------')
print('')
print('CURRENT ARRAY TASK ID: ' + SLURM_ARRAY_TASK_ID +' ' + "array_%A-%a.log")
print('SAMPLE NAME: ' + curr_sample)
print('NUMBER OF THREADS PER TASK: '+ threads)
print('TWO PASS MODE: ' + twopassMode)
print('')
print('----------------------------------------------------------------------')

cmd = ' '.join(base_cmd, default_cmd)

if testing:
    print(cmd)
else:
    subprocess.run(cmd, shell=True, check=True)

endSTAR_wallclock = .datetime.now()
STARruntime = endSTAR_wallclock - start_wallclock
print('STAR RUNTIME: ' + timedelta.total_seconds(STAR_runtime) + 'seconds')

# Post-processing
with cd(curr_outdir):
    # set permissions
    for r,d,f in os.walk(curr_prefix + '._STARpass1'):
        os.chmod(r, 0o755)

    # delete unneeded files
    shutil.rmtree(curr_prefix + '._STARgenome')
    if os.path.exists(curr_prefix + '._STARtmp'):
        shutil.rmtree(curr_prefix + '._STARtmp')

    # sort BAM (use samtools to get around the memory gluttony of STAR)
    print('[' + datetime.now().strftime("%b %d %H:%M:%S") + '] Sorting BAM', flush=True)
    cmd = 'samtools sort --threads ' + threads + ' -o ' + curr_prefix + '.Aligned.sortedByCoord.out.bam ' + curr_prefix + '.Aligned.out.bam'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    os.remove(curr_prefix + '.Aligned.out.bam')
    print('[' + datetime.now().strftime("%b %d %H:%M:%S") + '] Finished sorting BAM', flush=True)

    # index BAM
    print('['+datetime.now().strftime("%b %d %H:%M:%S") + '] Indexing BAM', flush=True)
    cmd = 'samtools index ' + curr_prefix + '.Aligned.sortedByCoord.out.bam'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    print('[' + datetime.now().strftime("%b %d %H:%M:%S") + '] Finished indexing BAM', flush=True)

    # sort and index chimeric BAM
    if int(chimSegmentMin)>0:
        cmd = 'samtools sort --threads ' + threads + ' -o ' + curr_prefix + '.Chimeric.out.sorted.bam ' + curr_prefix + '.Chimeric.out.sam'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        cmd = 'samtools index ' + curr_prefix + '.Chimeric.out.sorted.bam'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        os.remove(curr_prefix + '.Chimeric.out.sam')


endPostProcess_wallclock = datetime.datetime.now()
PostProcess_runtime = endPostProcess_wallclock - endSTAR_wallclock
total_elapsed = endPostProcess_wallclock - start_wallclock

print('POST_PROCESS RUNTIME: ' + datetime.timedelta.total_seconds(PostProcess_runtime) + 'seconds')
print('TOTAL ELEAPSED TIME: ' + datetime.timedelta.total_seconds(total_elapsed) + 'seconds')

datetime.datetime.now


