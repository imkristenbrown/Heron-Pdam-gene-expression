# BLAST-ing our genes to the Biomineralization Toolkit


On andromeda

Download protein file from genome

```
cd /data/putnamlab/zdellaert/Pdam-TagSeq/references

wget http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.genes.pep.faa.gz

gunzip Pocillopora_acuta_HIv2.genes.pep.faa.gz

cd ..

mkdir blast
cd blast
```

On personal computer

```
scp  /Users/zoedellaert/Documents/URI/Heron-Pdam-gene-expression/BioInf/data/Biomineralization_Toolkit_FScucchia/Biomineralization_Toolkit_FScucchia.fasta zdellaert@ssh3.hac.uri.edu:/data/putnamlab/zdellaert/Pdam-TagSeq/blast/
Biomineralization_Toolkit_FScucchia.fasta
```

On andromeda

```
nano Biomineralization_blast.sh
```

```
#!/bin/bash
#SBATCH --job-name="Pacuta_TRP_blast"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --error="blast_out_error"
#SBATCH --output="blast_out"
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/zdellaert/Pdam-TagSeq/blast/
#SBATCH --nodes=1 --ntasks-per-node=20

module load BLAST+/2.9.0-iimpi-2019b

makeblastdb -in ../references/Pocillopora_acuta_HIv2.genes.pep.faa -out Pacuta_prot -dbtype prot

blastp -query Biomineralization_Toolkit_FScucchia.fasta -db Pacuta_prot -out Biomineralization_blast_results.txt -outfmt 0

blastp -query Biomineralization_Toolkit_FScucchia.fasta -db Pacuta_prot -out Biomineralization_blast_results_tab.txt -outfmt 6 -max_target_seqs 1
```


```
sbatch Biomineralization_blast.sh
```

Errors:



On personal computer:

```
scp  zdellaert@ssh3.hac.uri.edu:/data/putnamlab/zdellaert/Pdam-TagSeq/blast/Biomineralization_blast_results.txt /Users/zoedellaert/Documents/URI/Heron-Pdam-gene-expression/BioInf/output

scp  zdellaert@ssh3.hac.uri.edu:/data/putnamlab/zdellaert/Pdam-TagSeq/blast/Biomineralization_blast_results_tab.txt /Users/zoedellaert/Documents/URI/Heron-Pdam-gene-expression/BioInf/output
```