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
#module load GCCcore/9.3.0 #I needed to add this to resolve conflicts between loaded GCCcore/7.3.0 and GCCcore/9.3.0
module load Python/2.7.15-foss-2018b #Python
module load StringTie/2.1.4-GCC-9.3.0 #Transcript assembly: StringTie
module load GffCompare/0.12.1-GCCcore-8.3.0 #Transcript assembly QC: GFFCompare

#make gtf_list.txt file
ls *.gtf > gtf_list.txt

stringtie --merge -p 8 -G /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes_fixed.gff3 -o #HeronPdam_merged.gtf gtf_list.txt #Merge GTFs 
#echo "Stringtie merge complete" $(date)

gffcompare -r /data/putnamlab/zdellaert/Pdam-TagSeq/references/Pocillopora_acuta_HIv2.genes_fixed.gff3 -G -o merged #HeronPdam_merged.gtf #Compute the accuracy and compare with reference annotation
echo "GFFcompare complete, Starting gene count matrix assembly..." $(date)

#Re-Estimate transcript abundances using HeronPdam_merged.gtf, as in example: "stringtie –e –B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf ERR188044_chrX.bam"

cd /data/putnamlab/zdellaert/Pdam-TagSeq/processed/aligned_Pacuta

array=($(ls *_ALL.bam)) #Make an array of sequences to assemble
 
for i in ${array[@]}; do 
        
        stringtie -e -p 8 -G /data/putnamlab/zdellaert/Pdam-TagSeq/Stringtie2/HeronPdam_merged.gtf -o /data/putnamlab/zdellaert/Pdam-TagSeq/Stringtie2/${i}.merge.gtf ${i}
        
        echo "StringTie re-estimation of abundance for seq file ${i}" $(date)
done

#make gtf list text file

cd /data/putnamlab/zdellaert/Pdam-TagSeq/Stringtie2

for filename in *merge.gtf; do echo $filename $PWD/$filename; done > listGTF.txt

python ../scripts/prepDE.py -g HeronPdam_gene_count_matrix.csv -i listGTF.txt #Compile the gene count matrix

cp prepDE_* ../scripts/script_outputs/ #copy script outputs to script_outputs folder

echo "Gene count matrix compiled." $(date)