# Cataglyphis hispanica

## Scripts to look at TEdenovo with REPET package v2.5, filter out potential contamination, Repeat Masking and diversity plot

## Code for:
Darras, Hugo; De Souza Araujo, Natalia; Baudry, Lyam; Guiglielmoni, Nadège; Lorite, Pedro; Marbouty, Martial; Rodriguez, Fernando; Arkhipova, Irina; Koszul, Romain; Flot, Jean-François; Aron, Serge. Chromosome-level genome assembly and annotation of two lineages of the ant Cataglyphis hispanica: stepping stones towards genomic studies of hybridogenesis and thermal adaptation in desert ants. Peer Community Journal, Volume 2 (2022), article no. e40. doi : 10.24072/pcjournal.140. https://peercommunityjournal.org/articles/10.24072/pcjournal.140/

## Transposons with TEdenovo

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

## RUN runTEdenovo

## RUN TEcleanup

## RUN RepeatMasking
