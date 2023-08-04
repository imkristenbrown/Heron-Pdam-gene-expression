# Confirming species ID for Pocillopora gene expression analysis
[NCBI Fasta Download](https://github.com/kblin/ncbi-acc-download)

```
pip install ncbi-acc-download

cd Heron-Pdam-gene-expression/BioInf/species_id
```

#Pacuta IDs
```
ncbi-acc-download --format fasta JX994073 JX624999 KF583928  KF583935 EU374226 FJ424111 KJ720240 KM215075 KJ690906 KP698585 KX538985 KX538986 KM215104 KJ720241
```


#Pdam IDs 
```
ncbi-acc-download --format fasta JX994077 KJ720219 JX994086 JX624991 KF583925 KF583946 EU374235 KJ720218 KP698587 KM215098 KX538982 KJ720226 KX538983 KF583950 JX994087 KJ720235 JX985618 KJ690905 JX625025 KX538984
```
```
cat *.fa > pacuta_pdam.fasta
```
add Brown et al consensus seqeunce = [Brown et al 2022](https://doi.org/10.1098/rspb.2022.0941)
GenBank Accession numbers OP296503 – OP296521

add Strand consensus seqeunce = [Stephens et al 2022](10.1093/gigascience/giac098) P acuta genome sample 

Trim to all bases in all samples

Align in Genious Prime with MUSCLE

Tree with Genious Prime defaults

outcome = Pdamicornis species ID for [Brown et al 2022](https://doi.org/10.1098/rspb.2022.0941)

#Pdam Pacuta species tree

![Tree](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/species_id/pocillopora_pdam_pacuta_tree.png?raw=true)

## 4 new samples added by ZD, [PCR performed 04112023](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/Barott-and-Brown-Pdam-mtORF-SpeciesID/) and Sanger Sequencing on 04192023 of samples RF16A, RF24B, RS2C, and RS11D


Using same protocol, generated consensus sequence for these 4 samples (showed 100% alignment with the Brown et al consensus seqeunce of the other samples) and aligned to the Pdam, Pacuta, Brown et al consensus sequence, and Strand consensus sequence as detailed above.


#Pdam Pacuta species tree, with all samples

![Tree](https://github.com/imkristenbrown/Heron-Pdam-gene-expression/blob/master/BioInf/species_id/041923%20mtORF/pocillopora_ZD_mtORF_tree.png?raw=true)
