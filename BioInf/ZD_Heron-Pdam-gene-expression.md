# Heron-Pdam-gene-expression

Script Written By: Putnam, Brown, Dellaert
Last Updated: 20230204

[BaseSpace Data Info for CLI](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-examples)
Run ID: 376886519
Project ID: JA22512

## Download Data

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

## QC raw files

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

## TagSeq Pipeline based on [Dr. Ariana Huffmyer's pipeline](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/Genome_V3/TagSeq_BioInf_genomeV3.md) and [Dr. Sam Gurr's pipeline](https://github.com/SamGurr/SamGurr.github.io/blob/master/_posts/2021-01-07-Geoduck-TagSeq-Pipeline.md)

## Trimming and Trimmed QC

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

mv multiqc_report.html trimmed_qc/ #move output files to the QC directory
mv multiqc_data trimmed_qc/ #move output files to the QC directory

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

**Data look great after trimming, but we should keep an eye out for samples RF16A and RF16C, both of which had high sequence duplication levels.**

## Download Genome: [*Pocillopora damicornis*](https://www.ncbi.nlm.nih.gov/nuccore/RCHS00000000)

Cunning et al. 2018 [Publication](https://www.nature.com/articles/s41598-018-34459-8)

Obtain reference genome assembly and gff annotation file.

```
mkdir references/ 
cd references/ 

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/095/GCF_003704095.1_ASM370409v1/GCF_003704095.1_ASM370409v1_genomic.fna.gz #download reference genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/095/GCF_003704095.1_ASM370409v1/GCF_003704095.1_ASM370409v1_genomic.gff.gz #download genome annotation file

gunzip GCF_003704095.1_ASM370409v1_genomic.fna.gz #unzip genome file
gunzip GCF_003704095.1_ASM370409v1_genomic.gff.gz #unzip gff annotation file
```

## Alignment with HISAT2

### Concatenate L001 and L002 reads for each sample before alignment

#### This code is from Dr. Kevin Wong's [pipeline](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/scripts/TagSeq/TagSeq_Analysis_HPC.md)

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq #Enter working directory
nano scripts/cat.clean.sh #make script for concatenation, enter text in next code chunk
```

Note: I added a backslash before the ls *R1_001.fastq.gz because I have an alias for ls in my bash.profile

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --account=putnamlab   
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/processed
#SBATCH --mem=100GB

module load SAMtools/1.9-foss-2018b

echo "Concatenate files" $(date)

\ls *R1_001.fastq.gz | awk -F '[_]' '{print $1"_"$2}' | sort | uniq > ID #Create list of the sample IDs in this format: "clean.RF25C_S110"

for i in `cat ./ID`;
 do cat $i\_L001_R1_001.fastq.gz $i\_L002_R1_001.fastq.gz > $i\_ALL.fastq.gz;
 done

echo "Mission complete." $(date)
```

```
sbatch /data/putnamlab/zdellaert/Pdam-TagSeq/scripts/cat.clean.sh
```

#### Initiated 20230203 sbatch job id 222008

### Alignment with HISAT2: Pdam genome

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq #Enter working directory
nano scripts/align.sh #make script for alignment, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --account=putnamlab              
#SBATCH --error="align_error" #if your job fails, the error report will be put in this file
#SBATCH --output="align_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/processed

# load modules needed
module load HISAT2/2.2.1-foss-2019b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

# index the reference genome for Pocillopora damicornis output index to working directory
hisat2-build -f /data/putnamlab/zdellaert/Pdam-TagSeq/references/GCF_003704095.1_ASM370409v1_genomic.fna ./Pdam_ref # called the reference genome (scaffolds)
echo "Referece genome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed

array=($(ls /data/putnamlab/zdellaert/Pdam-TagSeq/processed/*_ALL.fastq.gz)) # call the clean sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
 hisat2 -p 8 --dta -x Pdam_ref -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
      echo "${i} bam-ified!"
        rm ${sample_name}.sam
done

cp align_* ../scripts/script_outputs/ #copy script outputs to script_outputs folder
```

```
sbatch /data/putnamlab/zdellaert/Pdam-TagSeq/scripts/align.sh
```

#### Initiated Alignment 20230203 sbatch job id 222014

- Nodes: 1
- Cores per node: 10
- CPU Utilized: 05:59:59
- CPU Efficiency: 65.10% of 09:13:00 core-walltime
- Job Wall-clock time: 00:55:18
- Memory Utilized: 2.87 GB
- Memory Efficiency: 2.87% of 100.00 GB

#### To view mapping percentages

```
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

for i in *.bam; do
    echo "${i}" >> mapped_reads_counts_Pdam
    samtools flagstat ${i} | grep "mapped (" >> mapped_reads_counts_Pdam
done
```

Average mapping rate: 29%

- min 24.01%
- max 36.06%
- average 29.10%
- count 48

### Alternative genome to map to: [*Pocillopora acuta*](http://cyanophora.rutgers.edu/Pocillopora_acuta/)

Rutgers University Stephens et al. 2022 [Publication](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac098/6815755)

Obtain reference genome assembly and gff annotation file.

```
mkdir references/
cd references/ 

wget http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.assembly.fasta.gz

wget http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.genes.gff3.gz

gunzip Pocillopora_acuta_HIv2.assembly.fasta.gz #unzip genome file
gunzip Pocillopora_acuta_HIv2.genes.gff3.gz #unzip gff annotation file
```

### Alignment with HISAT2 to *P. acuta* genome

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq #Enter working directory
mkdir processed/aligned_Pacuta
nano scripts/align_Pacuta.sh #make script for alignment, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                 
#SBATCH --error="align_Pacuta_error" #if your job fails, the error report will be put in this file
#SBATCH --output="align_Pacuta_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/processed

# load modules needed
module load HISAT2/2.2.1-foss-2019b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

# index the reference genome for Pocillopora acuta output index to working directory
hisat2-build -f /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.assembly.fasta ./Pacuta_ref # called the reference genome (scaffolds)
echo "Referece genome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed

array=($(ls /data/putnamlab/zdellaert/Pdam-TagSeq/processed/*_ALL.fastq.gz)) # call the clean sequences - make an array 

for i in ${array[@]}; do
    sample_name=`echo $i| awk -F [.] '{print $2}'`
    hisat2 -p 8 --dta -x Pacuta_ref -U ${i} -S /data/putnamlab/zdellaert/Pdam-TagSeq/processed/aligned_Pacuta/${sample_name}.sam
        samtools sort -@ 8 -o /data/putnamlab/zdellaert/Pdam-TagSeq/processed/aligned_Pacuta/${sample_name}.bam /data/putnamlab/zdellaert/Pdam-TagSeq/processed/aligned_Pacuta/${sample_name}.sam
      echo "${i} bam-ified!"
        rm /data/putnamlab/zdellaert/Pdam-TagSeq/processed/aligned_Pacuta/${sample_name}.sam
done

cp align_Pacuta_* ../scripts/script_outputs/ #copy script outputs to script_outputs folder
```

```
sbatch /data/putnamlab/zdellaert/Pdam-TagSeq/scripts/align_Pacuta.sh
```

Moved the following files into the aligned_Pacuta folder so that bam files and reference files were together:

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq/processed #Enter working directory
mv align_Pacuta_* aligned_Pacuta/
mv Pac* aligned_Pacuta/
```

#### Initiated Alignment 20230203 sbatch job id #SBATCH 222015

- Nodes: 1
- Cores per node: 10
- CPU Utilized: 07:49:16
- CPU Efficiency: 60.97% of 12:49:40 core-walltime
- Job Wall-clock time: 01:16:58
- Memory Utilized: 4.12 GB
- Memory Efficiency: 4.12% of 100.00 GB

#### To view mapping percentages

```
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

cd /data/putnamlab/zdellaert/Pdam-TagSeq/processed/aligned_Pacuta #Enter working directory

for i in *.bam; do
    echo "${i}" >> mapped_reads_counts_Pacuta
    samtools flagstat ${i} | grep "mapped (" >> mapped_reads_counts_Pacuta
done
```

Average mapping rate: 72.48%

- min 59.41%
- max 77.98%
- average 72.48%
- count 48

## Assembly with Stringtie 2, using reads aligned to *P. acuta* genome and *P. acuta* GFF3 file

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq #Enter working directory
mkdir Stringtie2 #make directory for assembly results
nano scripts/stringtie.sh #make script for assembly, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --account=putnamlab              
#SBATCH --error="stringtie_error" #if your job fails, the error report will be put in this file
#SBATCH --output="stringtie_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/processed/aligned_Pacuta

#load packages
module load StringTie/2.1.4-GCC-9.3.0

#Transcript assembly: StringTie

array=($(ls *_ALL.bam)) #Make an array of sequences to assemble
 
for i in ${array[@]}; do #Running with the -e option to compare output to exclude novel genes. Also output a file with the gene abundances
        sample_name=`echo $i| awk -F [_] '{print $1"_"$2}'` #use only "RS11B_S86" part of name, removes text after second underscore
  stringtie -p 8 -e -B -G /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes.gff3 -A ../../Stringtie2/${sample_name}.gene_abund.tab -o ../../Stringtie2/${sample_name}.gtf ${i}
        echo "StringTie assembly for seq file ${i}" $(date)
done

cp stringtie_* ../../scripts/script_outputs/ #copy script outputs to script_outputs folder

echo "StringTie assembly COMPLETE, starting assembly analysis" $(date)
```

-p means number of threads/CPUs to use (8 here)

-e means only estimate abundance of given reference transcripts (only genes from the genome) - dont use if using splice variance aware to merge novel and ref based.

-B means enable output of ballgown table files to be created in same output as GTF

-G means genome reference to be included in the merging

```
sbatch /data/putnamlab/zdellaert/Pdam-TagSeq/scripts/stringtie.sh
```

This will make a .gtf file for each sample.

#### Initiated Assembly 20230204 sbatch job id #SBATCH 223879

- Nodes: 1
- Cores per node: 20
- CPU Utilized: 00:17:04
- CPU Efficiency: 6.14% of 04:38:00 core-walltime
- Job Wall-clock time: 00:13:54
- Memory Utilized: 164.22 MB
- Memory Efficiency: 0.08% of 200.00 GB

## Generate gene count matrix Prep DE

### First, need to download script prepDE.py and add to scripts folder

1. This can be downloaded from the [Stringtie github repository](https://github.com/gpertea/stringtie/blob/master/prepDE.py) and uploaded to andromeda into your scripts folder
2. Alternatively, copy the script from github and create your own script
3. 3rd option, copy from another user in the lab :) (see below)

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq #Enter working directory
cp /data/putnamlab/ashuffmyer/pairs-rnaseq/prepDE.py scripts
nano scripts/prepDE.sh #make script for assembly, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --account=putnamlab              
#SBATCH --error="prepDE_error" #if your job fails, the error report will be put in this file
#SBATCH --output="prepDE_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/Stringtie2

#load packages
module load Python/2.7.15-foss-2018b #Python
module load StringTie/2.1.4-GCC-9.3.0 #Transcript assembly: StringTie
module load GffCompare/0.12.1-GCCcore-8.3.0 #Transcript assembly QC: GFFCompare

#make gtf_list.txt file
ls *.gtf > gtf_list.txt

stringtie --merge -p 8 -G /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes.gff3 -o HeronPdam_merged.gtf gtf_list.txt #Merge GTFs 
echo "Stringtie merge complete" $(date)

gffcompare -r /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes.gff3 -G -o merged HeronPdam_merged.gtf #Compute the accuracy 
echo "GFFcompare complete, Starting gene count matrix assembly..." $(date)

#make gtf list text file
for filename in *.gtf; do echo $filename $PWD/$filename; done > listGTF.txt

python ../scripts/prepDE.py -g HeronPdam_gene_count_matrix.csv -i listGTF.txt #Compile the gene count matrix

cp prepDE_* ../scripts/script_outputs/ #copy script outputs to script_outputs folder

echo "Gene count matrix compiled." $(date)
```

```
sbatch /data/putnamlab/zdellaert/Pdam-TagSeq/scripts/prepDE.sh
```

#### Initiated prepDE 20230204 sbatch job id #SBATCH 223882

### Error here, likely due to the GFF3 issues mentioned in [Ariana's pipeline](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/Genome_V3/TagSeq_BioInf_genomeV3.md) and in [this github issue](https://github.com/Putnam-Lab/Lab_Management/issues/51) or [this one](https://github.com/Putnam-Lab/Lab_Management/issues/11) however my script is not creating any errors...

### My issue: no gene count matrix is created.
```
cat prepDE_error
```
```
GCC/9.3.0(24):ERROR:150: Module 'GCC/9.3.0' conflicts with the currently loaded module(s) 'GCC/7.3.0-2.30'
GCC/9.3.0(24):ERROR:102: Tcl command execution failed: conflict GCC

GCCcore/9.3.0(24):ERROR:150: Module 'GCCcore/9.3.0' conflicts with the currently loaded module(s) 'GCCcore/7.3.0'
GCCcore/9.3.0(24):ERROR:102: Tcl command execution failed: conflict GCCcore

GCCcore/8.3.0(24):ERROR:150: Module 'GCCcore/8.3.0' conflicts with the currently loaded module(s) 'GCCcore/7.3.0'
GCCcore/8.3.0(24):ERROR:102: Tcl command execution failed: conflict GCCcore

  33730 reference transcripts loaded.
  33730 query transfrags loaded.
```

# Fixing GFF3 file in hopes of fixing gene count matrix issue as described by Ariana in [her pipeline](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/Genome_V3/TagSeq_BioInf_genomeV3.md) and in [this github issue](https://github.com/Putnam-Lab/Lab_Management/issues/51)

### The Mcap and Pacuta genomes both came from the same source, so we need to adjust the GFF3 file for Pacuta in the same way.

1. Download *P. acuta* gff3 from andromeda (or directly from [Rutgers](http://cyanophora.rutgers.edu/Pocillopora_acuta/))

```
scp  -r zdellaert@ssh3.hac.uri.edu:/data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes.gff3 /Users/zoedellaert/Documents/URI/Heron-Pdam-gene-expression/BioInf
```

2. In R on your computer, use [this script](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/fix_gff_format.Rmd), modified from [Dr. Ariana Huffmyer](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/Genome_V3/fix_gff_format.Rmd) to correct the Pacuta GFF3 file. 

3. Upload fixed *P. acuta* gff3 to andromeda

```
scp /Users/zoedellaert/Documents/URI/Heron-Pdam-gene-expression/BioInf/Pocillopora_acuta_HIv2.genes_fixed.gff3.gz zdellaert@ssh3.hac.uri.edu:/data/putnamlab/zdellaert/Pdam-TagSeq/references/
```

4. Unzip fixed gff3 file in andromeda

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq/references #Enter working directory
gunzip Pocillopora_acuta_HIv2.genes_fixed.gff3.gz
```

5. Redo pipeline from Stringtie onwards:

## Assembly with Stringtie 2, using reads aligned to *P. acuta* genome and **Fixed** *P. acuta* GFF3 file

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq #Enter working directory
mkdir Stringtie2 #make directory for assembly results
nano scripts/stringtie.sh #make script for assembly, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --account=putnamlab              
#SBATCH --error="stringtie_error" #if your job fails, the error report will be put in this file
#SBATCH --output="stringtie_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/processed/aligned_Pacuta

#load packages
module load StringTie/2.1.4-GCC-9.3.0

#Transcript assembly: StringTie

array=($(ls *_ALL.bam)) #Make an array of sequences to assemble
 
for i in ${array[@]}; do #Running with the -e option to compare output to exclude novel genes. Also output a file with the gene abundances
        sample_name=`echo $i| awk -F [_] '{print $1"_"$2}'` #use only "RS11B_S86" part of name, removes text after second underscore
  stringtie -p 8 -e -B -G /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes_fixed.gff3 -A ../../Stringtie2/${sample_name}.gene_abund.tab -o ../../Stringtie2/${sample_name}.gtf ${i}
        echo "StringTie assembly for seq file ${i}" $(date)
done

cp stringtie_* ../../scripts/script_outputs/ #copy script outputs to script_outputs folder

echo "StringTie assembly COMPLETE, starting assembly analysis" $(date)
```

-p means number of threads/CPUs to use (8 here)

-e means only estimate abundance of given reference transcripts (only genes from the genome) - dont use if using splice variance aware to merge novel and ref based.

-B means enable output of ballgown table files to be created in same output as GTF

-G means genome reference to be included in the merging

```
sbatch /data/putnamlab/zdellaert/Pdam-TagSeq/scripts/stringtie.sh
```

This will make a .gtf file for each sample.

#### Initiated Assembly 20230204 sbatch job id #SBATCH 223895

- Nodes: 1
- Cores per node: 20
- CPU Utilized: 00:16:59
- CPU Efficiency: 6.32% of 04:28:40 core-walltime
- Job Wall-clock time: 00:13:26
- Memory Utilized: 167.74 MB
- Memory Efficiency: 0.08% of 200.00 GB

## Generate gene count matrix Prep DE

### First, need to download script prepDE.py and add to scripts folder

1. This can be downloaded from the [Stringtie github repository](https://github.com/gpertea/stringtie/blob/master/prepDE.py) and uploaded to andromeda into your scripts folder
2. Alternatively, copy the script from github and create your own script (see below)
3. 3rd option, copy from another user in the lab :)
        - `cp /data/putnamlab/ashuffmyer/pairs-rnaseq/prepDE.py scripts`

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq #Enter working directory
nano scripts/prepDE.py #paste in code from github page above
nano scripts/prepDE.sh #make script for assembly, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --account=putnamlab              
#SBATCH --error="prepDE_error" #if your job fails, the error report will be put in this file
#SBATCH --output="prepDE_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/Stringtie2

#load packages
module load Python/2.7.15-foss-2018b #Python
module load StringTie/2.1.4-GCC-9.3.0 #Transcript assembly: StringTie
module load GffCompare/0.12.1-GCCcore-8.3.0 #Transcript assembly QC: GFFCompare

#make gtf_list.txt file
ls *.gtf > gtf_list.txt

stringtie --merge -p 8 -G /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes_fixed.gff3 -o HeronPdam_merged.gtf gtf_list.txt #Merge GTFs 
echo "Stringtie merge complete" $(date)

gffcompare -r /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes_fixed.gff3 -G -o merged HeronPdam_merged.gtf #Compute the accuracy 
echo "GFFcompare complete, Starting gene count matrix assembly..." $(date)

#make gtf list text file
for filename in *.gtf; do echo $filename $PWD/$filename; done > listGTF.txt

python ../scripts/prepDE.py -g HeronPdam_gene_count_matrix.csv -i listGTF.txt #Compile the gene count matrix

cp prepDE_* ../scripts/script_outputs/ #copy script outputs to script_outputs folder

echo "Gene count matrix compiled." $(date)
```

```
sbatch /data/putnamlab/zdellaert/Pdam-TagSeq/scripts/prepDE.sh
```

#### Initiated prepDE 20230204 sbatch job id 223896

### Note: This did not fix the missing gene count table issue. However, fixing the gff3 file is still likely needed so I will rewrite a smoother workflow above.




#### To do once csv file is created:

Export gene count matrix report to my computer using scp and upload results to [GitHub repository](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/HeronPdam_gene_count_matrix.csv).

```
scp  zdellaert@ssh3.hac.uri.edu:/data/putnamlab/zdellaert/Pdam-TagSeq/Stringtie2/HeronPdam_gene_count_matrix.csv /Users/zoedellaert/Documents/URI/Heron-Pdam-gene-expression/BioInf/HeronPdam_gene_count_matrix.csv
```