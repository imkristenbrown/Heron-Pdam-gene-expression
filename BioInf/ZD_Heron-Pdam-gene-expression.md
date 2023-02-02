# Heron-Pdam-gene-expression
Script Written By: Putnam, Brown, Dellaert
Last Updated: 20230202

[BaseSpace Data Info for CLI](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-examples)
Run ID: 376886519
Project ID: JA22512

# Download Data
Work will be done on URI server Andromeda

```
mkdir /data/putnamlab/KITT/hputnam/20230125_Barott_Pdam
cd /data/putnamlab/KITT/hputnam/20230125_Barott_Pdam
mkdir scripts
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20230125_Barott_Pdam

module load IlluminaUtils/2.11-GCCcore-9.3.0-Python-3.8.2

bs download project -i 376886519 -o /data/putnamlab/KITT/hputnam/20230125_Barott_Pdam
```

# QC raw files

```
nano /data/putnamlab/KITT/hputnam/20230125_Barott_Pdam/scripts/qc.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20230125_Barott_Pdam

module load FastQC/0.11.9-Java-11 
module load MultiQC/1.9-intel-2020a-Python-3.8.2

#make qc output folder
mkdir raw_qc/

#run fastqc on raw data
fastqc *.fastq.gz -o raw_qc/

#Compile MultiQC report from FastQC files
multiqc ./raw_qc
mv multiqc_report.html raw_qc/raw_qc_multiqc_report.html
mv multiqc_data raw_qc/raw_multiqc_data

echo "Initial QC of Seq data complete." $(date)
```

```
sbatch /data/putnamlab/KITT/hputnam/20230125_Barott_Pdam/scripts/qc.sh
```

### initiated QC 20230125 - slurm-216067.out

MultiQC report (.html formats) can be found in the [GitHub repository](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/raw_qc)

# TagSeq Pipeline based on [Dr. Ariana Huffmyer's pipeline](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/Genome_V3/TagSeq_BioInf_genomeV3.md) and [Dr. Sam Gurr's pipeline](https://github.com/SamGurr/SamGurr.github.io/blob/master/_posts/2021-01-07-Geoduck-TagSeq-Pipeline.md)

# Trimming and Trimmed QC

Working on URI Andromeda Server

Raw data are in:

> /data/putnamlab/KITT/hputnam/20230125_Barott_Pdam/

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq #Enter working directory
mkdir raw_data #make folder for raw data
mkdir scripts #make folder for scripts

ln -s /data/putnamlab/KITT/hputnam/20230125_Barott_Pdam/*.fastq.gz /data/putnamlab/zdellaert/Pdam-TagSeq/raw_data/ #symlink the raw reads into my directory
```

Raw data are now sym-linked in:

> /data/putnamlab/zdellaert/Pdam-TagSeq/raw_data/


This script will:

1. Generate FastQC/MultiQC for raw sequences
2. Conduct trimming and cleaning

    Settings:
    - remove adapter sequences --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    - enable polyX trimming on 3' end at length of 6 --trim_poly_x 6
    - filter by minimum phred quality score of >30  -q 30
    - enable low complexity filter -y
    - set complexity filter threshold of 50% required -Y 50

3. Generate reports for cleaned sequences.

MultiQC report (.html formats) can be found in the [GitHub repository](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/trimmed_qc)

```
nano scripts/trim_qc.sh #make script for trimming and QC, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --error=../"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/raw_data

# load modules needed
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

#make qc output folder
mkdir /data/putnamlab/zdellaert/Pdam-TagSeq/trimmed_qc/

#make processed folder for trimmed reads
mkdir /data/putnamlab/zdellaert/Pdam-TagSeq/processed/

# Make an array of sequences to trim
array1=($(ls *.fastq.gz)) 

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA) 

for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 /data/putnamlab/zdellaert/Pdam-TagSeq/processed/clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 
        fastqc /data/putnamlab/zdellaert/Pdam-TagSeq/processed/clean.${i} -o /data/putnamlab/zdellaert/Pdam-TagSeq/trimmed_qc/ # fastqc the cleaned reads
done 

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

multiqc /data/putnamlab/zdellaert/Pdam-TagSeq/trimmed_qc/ #Compile MultiQC report from FastQC files 

mv multiqc_* trimmed_qc/

mv raw_data/multiqc_report.html trimmed_qc/ #move output files to the QC directory
mv raw_data/multiqc_data trimmed_qc/ #move output files to the QC directory

echo "Cleaned MultiQC report generated." $(date)
```

```
sbatch /data/putnamlab/zdellaert/Pdam-TagSeq/scripts/trim_qc.sh
```

Exported multiQC report to my computer using scp and uploaded results to [GitHub repository](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/trimmed_qc).

```
scp  zdellaert@ssh3.hac.uri.edu:/data/putnamlab/zdellaert/Pdam-TagSeq/trimmed_qc/multiqc_report.html /Users/zoedellaert/Documents/URI/Heron-Pdam-gene-expression/BioInf/trimmed_qc/trimmed_multiqc_report.html

scp  -r zdellaert@ssh3.hac.uri.edu:/data/putnamlab/zdellaert/Pdam-TagSeq/trimmed_qc/multiqc_data /Users/zoedellaert/Documents/URI/Heron-Pdam-gene-expression/BioInf/trimmed_qc/trimmed_multiqc_data
```

### Initiated trimming and QC 20230202 sbatch job id 221857