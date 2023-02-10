# Alternative version of Stringtie.

## Full pipeline with re-estimation and merge: Assembly with Stringtie 2, using reads aligned to *P. acuta* genome and **Fixed** *P. acuta* GFF3 file

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq #Enter working directory
mkdir FullP_Stringtie2 #make directory for assembly results
nano scripts/FullP_stringtie.sh #make script for assembly, enter text in next code chunk
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
#SBATCH --error="FullP_stringtie_error" #if your job fails, the error report will be put in this file
#SBATCH --output="FullP_stringtie_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/processed/aligned_Pacuta

#load packages
module load StringTie/2.1.4-GCC-9.3.0

#Transcript assembly: StringTie using full pipeline without -e and -B in initial assembly

array=($(ls *_ALL.bam)) #Make an array of sequences to assemble
 
for i in ${array[@]}; do #Running without the -e and -B
        sample_name=`echo $i| awk -F [_] '{print $1"_"$2"_"$3}'`
        stringtie -p 8 -G /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes_fixed.gff3 -A ../../FullP_Stringtie2/${sample_name}.gene_abund.tab -o ../../FullP_Stringtie2/${sample_name}.gtf ${i}
        echo "StringTie assembly for seq file ${i}" $(date)
done

cp FullP_stringtie_* ../../scripts/script_outputs/ #copy script outputs to script_outputs folder

echo "StringTie assembly COMPLETE, starting assembly analysis" $(date)
```

-p means number of threads/CPUs to use (8 here)

-e means only estimate abundance of given reference transcripts (only genes from the genome) - dont use if using splice variance aware to merge novel and ref based. #removing from this step for full pipeline, see Figure [in docs](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#de).

-B means enable output of ballgown table files to be created in same output as GTF #removing from this step for full pipeline, see Figure [in docs](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#de).

-G means genome reference to be included in the merging

```
sbatch /data/putnamlab/zdellaert/Pdam-TagSeq/scripts/FullP_stringtie.sh
```

This will make a .gtf file for each sample.

#### Initiated Assembly 20230206 sbatch job id #SBATCH 224097

- Nodes: 1
- Cores per node: 20
- CPU Utilized: 00:19:41
- CPU Efficiency: 7.52% of 04:21:40 core-walltime
- Job Wall-clock time: 00:13:05
- Memory Utilized: 129.48 MB
- Memory Efficiency: 0.06% of 200.00 GB

## Generate gene count matrix Prep DE

### First, need to download script prepDE.py and add to scripts folder

1. This can be downloaded from the [Stringtie github repository](https://github.com/gpertea/stringtie/blob/master/prepDE.py) and uploaded to andromeda into your scripts folder
2. Alternatively, copy the script from github and create your own script (see below)
3. 3rd option, copy from another user in the lab :)
        - `cp /data/putnamlab/ashuffmyer/pairs-rnaseq/prepDE.py scripts`

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq #Enter working directory
nano scripts/prepDE.py #paste in code from github page above
nano scripts/FullP_prepDE.sh #make script for assembly, enter text in next code chunk
```

**Adding in re-estimation step, this is the only time that stringtie --merge is recommended, it is unnecessary (and unused) in the pipeline above**

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --account=putnamlab              
#SBATCH --error="FullP_prepDE_error" #if your job fails, the error report will be put in this file
#SBATCH --output="FullP_prepDE_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/FullP_Stringtie2

#load packages
module load GCCcore/9.3.0 #I needed to add this to resolve conflicts between loaded GCCcore/7.3.0 and GCCcore/9.3.0
module load Python/2.7.15-foss-2018b #Python
module load StringTie/2.1.4-GCC-9.3.0 #Transcript assembly: StringTie
module load GffCompare/0.12.1-GCCcore-8.3.0 #Transcript assembly QC: GFFCompare

#make gtf_list.txt file
ls *.gtf > gtf_list.txt

stringtie --merge -e -p 8 -G /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes_fixed.gff3 -o HeronPdam_merged.gtf gtf_list.txt #Merge GTFs 
echo "Stringtie merge complete" $(date)

gffcompare -r /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes_fixed.gff3 -G -o merged HeronPdam_merged.gtf #Compute the accuracy and examine how the transcripts compare with the reference annotation (optional)
echo "GFFcompare complete, Starting gene count matrix assembly..." $(date)

#Note: the merged part is actually redundant and unnecessary unless we perform the original stringtie step without the -e function and perform
#re-estimation with -e after stringtie --merge, but will redo the pipeline later and confirm that I get equal results.

#stringtie –e –B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf ERR188044_chrX.bam

cd /data/putnamlab/zdellaert/Pdam-TagSeq/processed/aligned_Pacuta 

array=($(ls *_ALL.bam)) #Make an array of sequences to assemble

for i in ${array[@]}; do #Running with the -e and -B options for re-estimation
        sample_name=`echo $i| awk -F [_] '{print $1"_"$2"_"$3}'`
        stringtie -p 8 -e -B -G /data/putnamlab/zdellaert/Pdam-TagSeq/FullP_Stringtie2/HeronPdam_merged.gtf -A /data/putnamlab/zdellaert/Pdam-TagSeq/FullP_Stringtie2/${sample_name}.gene_abund.merged.tab -o /data/putnamlab/zdellaert/Pdam-TagSeq/FullP_Stringtie2/${sample_name}.merged.gtf ${i}
        echo "StringTie assembly for seq file ${i}" $(date)
done

cd /data/putnamlab/zdellaert/Pdam-TagSeq/FullP_Stringtie2

#make gtf list text file
for filename in *bam.merged.gtf; do echo $filename $PWD/$filename; done > listGTF.txt

#this could need to be .bam.merged.gtf, what does prepDE want after the merged case?

python ../scripts/prepDE.py -g HeronPdam_gene_count_matrix_FullPipeline.csv -i listGTF.txt #Compile the gene count matrix

cp FullP_prepDE_* ../scripts/script_outputs/ #copy script outputs to script_outputs folder

echo "Gene count matrix compiled." $(date)
```

```
sbatch /data/putnamlab/zdellaert/Pdam-TagSeq/scripts/FullP_prepDE.sh
```

### initiated prepDE 20230206 sbatch job id 224130

- Nodes: 1
- Cores per node: 20
- CPU Utilized: 00:22:29
- CPU Efficiency: 6.44% of 05:49:00 core-walltime
- Job Wall-clock time: 00:17:27
- Memory Utilized: 352.55 MB
- Memory Efficiency: 0.17% of 200.00 GB

### Export gene count matrix report to my computer using scp and upload results to [GitHub repository](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/tree/master/BioInf/TagSeq_Output/Different_Pipeline_Versions/HeronPdam_gene_count_matrix_FullPipeline.csv).

```
scp  zdellaert@ssh3.hac.uri.edu:/data/putnamlab/zdellaert/Pdam-TagSeq/FullP_Stringtie2/HeronPdam_gene_count_matrix_FullPipeline.csv /Users/zoedellaert/Documents/URI/Heron-Pdam-gene-expression/BioInf/HeronPdam_gene_count_matrix_FullPipeline.csv
```
