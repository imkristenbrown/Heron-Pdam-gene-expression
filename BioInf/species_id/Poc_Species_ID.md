https://github.com/kblin/ncbi-acc-download

pip install ncbi-acc-download

cd Heron-Pdam-gene-expression/BioInf/species_id

#Pacuta IDs
ncbi-acc-download --format fasta JX994073 JX624999 KF583928  KF583935 EU374226 FJ424111 KJ720240 KM215075 KJ690906 KP698585 KX538985 KX538986 KM215104 KJ720241



#Pdam IDs 
ncbi-acc-download --format fasta JX994077 KJ720219 JX994086 JX624991 KF583925 KF583946 EU374235 KJ720218 KP698587 KM215098 KX538982 KJ720226 KX538983 KF583950 JX994087 KJ720235 JX985618 KJ690905 JX625025 KX538984

cat *.fa > pacuta_pdam.fasta

add Brown et al consensus seqeunce

add Strand et al consensus seqeunce

Trim to all bases in all samples

Align in Genious Prime with MUSCLE

Tree with Genious Prime defaults

outcome = Pdamicornis species ID

[Tree]()