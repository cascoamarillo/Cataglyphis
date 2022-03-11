# Cataglyphis

## Scripts to look at TEdenovo with REPET package v2.5 (Quesneville et al., 2005; Flutre et al., 2011).

## Code for https://doi.org/10.1101/2022.01.07.475286

## Transposons

#TEdenovo using REPET v2 #==============================================================================

#TEdenovo.py from REPET v2 WORKS with phyton2:

#module load python/2.7.15

#You can run TEdenovo using the sequence provided as example "Chispanica.fa", or replace with your fasta file (assembly) of interest.

#RENAME the project name in the configuration file TEdenovo.cfg #[project] #project_name: Chispanica

#Edit TEdenovo.cfg and ./runTEdenovo in order to adapt it to your personal situation. More in https://urgi.versailles.inra.fr/Tools/REPET/TEdenovo-tuto

#The output is a MCL-filtered library of TEs located in {NAME}_Blaster_GrpRecPil_Struct_Map_TEclassif_Filtered_MCL

#==============================================================================

On the FASTA files you use: * There is a 15-character limit on the name of the file (and the project name) * Make sure there is a line break at least once every 60 bases * Make sure there are no illegal characters in the FASTA headers

##LINE BREAK (https://bioinf.shenwei.me/seqkit/) seqkit seq -w 60 Sp.fasta > Sp.fna

##RENAME fasta headers awk '/^>/{print ">Sp_ctg" ++i; next}{print}' < Sp.fna > Sp.fst

##ALL CAPITAL LETTERS awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' Sp.fst > Sp.fa

## Repeat Classification of TE libary

#Run RepeatClassifier (Repeat Modeler) to generate a .classified file

RepeatClassifier -consensi ${TElib}

## Extract #Unknown sequences

##There are different ways. One would be:

###Extract fasta headers
grep "^>" ${TElib}.classified > ${TElib}.classified.IDs.txt

###Sort by only #Unknown
${TElib}.classified.IDs.txt > ${TElib}.classified.IDs-Unknown.txt

###Extract fasta sequences for #Unknown and "Known"
seqkit grep --pattern-file ${TElib}.classified.IDs-Unknown.txt ${TElib}.classified > ${TElib}.classified.Unknown.fa
seqkit grep --pattern-file ${TElib}.classified.IDs-Known.txt ${TElib}.classified > ${TElib}.classified.Known.fa

### Compare "Unknown" sequences to Bacteria, mitochondria and CDS Cardiocondyla obscurior sets

####Bacteria
#ncbi online search can be used if the query file is not big (query length size to 1,000,000 or less), using in "Organism": bacteria (taxid:2). I.e. split original query file.
#Another way:
#Use https://github.com/ncbi/sra-tools
#and work with blast in remote. First, link WGS of interest (Bacteria taxID = 2)
taxid2wgs.pl -title "Bacteria WGS" -alias_file bacteria-wgs 2
#Use sra-blastn to search:
sra-blastn -query ${TElib}.classified.Unknown.fa -db bacteria-wgs -outfmt '6 qseqid sseqid pident evalue staxids stitle' -evalue 0.001 -num_alignments 1 -out Chis_LibTEs_Unknown_blast_bacteria.csv

####Cardiocondyla obscurior mitochondria set
#Genbank accession # KX951753
#https://www.ncbi.nlm.nih.gov/nuccore/KX951753
blastn -query ${TElib}.classified.Unknown.fa -db mt-Cobs_KX951753.1.fa -outfmt '6 qseqid sseqid pident evalue staxids stitle' -num_alignments 1 -out Chis_LibTEs_Unknown_blast_CobsCDS.csv -evalue 0.01

####Cardiocondyla obscurior CDS set
#Available in http://hymenopteragenome.org/cardiocondyla/?q=genome_consortium_datasets
#CDS https://elsiklab-data.missouri.edu/data/hgd/cds_fasta/Cardiocondyla_obscurior_Cobs_1.4_OGSv1.4_cds.fa.gz
blastn -query ${TElib}.classified.Unknown.fa -db Cobs_1.4_OGS_CDS.fa -outfmt '6 qseqid sseqid pident evalue staxids stitle' -num_alignments 1 -out Chis_LibTEs_Unknown_blast_CobsCDS.csv -evalue 0.01

## Final TE library: merge ${TElib}.classified.Unknown.fa with unmatched ${TElib}.classified.Unknown.fa

