# Heron-Pdam-gene-expression
Script Written By: Putnam, Brown, Dellaert
Last Updated: 20230125

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