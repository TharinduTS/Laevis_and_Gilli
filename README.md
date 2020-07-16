# Laevis_and_Gilli

##  Index a reference genome, align and sort fastq files, and check the depth per base for a defined region. (had to run seperately without using loops for some fq files to work)

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load bwa
module load samtools/1.10
bwa index $1
bwa mem $1 $2 -t 4 | samtools view -Shu - | samtools sort - -o $2sorted.bam
samtools index $2sorted.bam
samtools depth -r chr1L:1-10000000 $2sorted.bam > $2sorted.bam.depth_per_base

```
# Do the same with bam files instead fq

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load bwa
module load samtools/1.10
samtools index $1
samtools depth -r chr1L:1-10000000 $1> $1sorted.bam.depth_per_base
```
# Concat files, Create a concat file, plot depth for all the samples

```r
library (ggplot2)
#set working directory here
working_directory<-"~/Desktop/OneDrive - McMaster University/for lab and research/Tharindu on Mac/lab/structure_lplots/depth_files"

setwd(working_directory)
#remove scientific notation
options(scipen=999)

#*******create concat file****

#set path with raw files here
raw_files_are_in<-working_directory

#create output directory
dir.create(paste(working_directory,"/out",sep = ""))

#save concat file in
save_path<-paste(working_directory,"/out",sep = "")

require(dplyr)
require(data.table)


# read file path raw
all_paths_raw_pre <-
  list.files(path = raw_files_are_in,
             full.names = TRUE,recursive = FALSE)
all_paths_raw<-all_paths_raw_pre[!file.info(all_paths_raw_pre)$isdir]

#create first data frame
x=1

all_data<-read.table(file=all_paths_raw[1])
colnames(all_data)<-c("chr","pos","dep")
all_data["sample_name"]=basename(all_paths_raw[1])
all_data$sample_no<-x

#add data
for (i in all_paths_raw[2:length(all_paths_raw)]) {
  x=x+1
  
  current_table<-read.table(i)
  setnames(current_table,c("chr","pos","dep"))
  current_table["sample_name"]=basename(i)
  current_table$sample_no<-x
  all_data<-rbind(all_data,current_table)
  
}
write.table(all_data,file =paste(save_path,"/concat",sep=''))

#*****end of concat******

#set new working directory
setwd(save_path)

all_data<-read.table(file = "concat")

p <- ggplot(data = all_data, aes(x=sample_no, y=dep))  +
  ylim(0,60)+
  ylab("Depth")+
  xlab("sample_number")+
  geom_boxplot(aes(fill=sample_name))
ggsave(filename = "all_samples.pdf",plot = p,height = 10,width = 20)
```

# ANGSD

## Load ANGSD
```bash
module load nixpkgs/16.09
module load intel/2018.3
module load angsd/0.929
```
## Set paths to the programs and the data 
(in server)

## NB this must be done every time you open a new terminal
```bash
AA=/scratch/premacht/xlaevis_and_xgilli/ANGSD
```
## Set path to a bam file list with several bam files
```bash
BAMFOLDER=$AA/../all_bam_files
```
# Make input data using ANGSD

## Make a file that contains the paths of all the bam files
```bash
find $BAMFOLDER |  grep bam$ > all.files
```
## Calculate GLs from all the BAM files using ANGSD
```bash
angsd -bam all.files -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -minMapQ 20 -minQ 20 -doCounts 1 -doDepth 1 -setMinDepth 2 -setMaxDepth 100  -minInd 15 -minMaf 0.05 -doGlf 2 -out all -P 1
```
Description
```txt
-bam all.files : tells ANGSD that the bam files to calculate GL from are listed in the file all.files
 -GL 2 : tells ANGSD to use the GATK genotype likelihood model
 -doMajorMinor 1 : tells ANGSD to infer the major and minor alleles
 -doMaf 1 : tells ANGSD to calculate minor allele frequencies (needed by two of the options below: -SNP_pval and -minMaf)
 -SNP_pval 2e-6 : tells ANGSD to use a p-value threshold of 2e-6 for calling SNPs
 -minMapQ 20 : tells ANGSD what to require as minimum mapping quality (quality filter)
 -minQ 20 : tells ANGSD what to require as minimum base quality (quality filter). One third of the sample size was used.
 -doCounts 1 :	Count the number A,C,G,T. All sites, All samples) - needed to use -doDepth
 -doDepth 1 : find depth distribution for every sample and for all samples jointly - needed to use depth limits
 -setMinDepth 2 :	If total depth is smaller then site is removed from analysis.
 -setMaxDepth 100 : If total depth is larger then site is removed from analysis.
 -minInd 25 : tells ANGSD to only output GLs for loci with data for at least 25 individuals for each site (quality filter)
 -minMaf 0.05 : tells ANGSD to only output GLS for loci with a minimum minor allele frequency of 0.05 (quality filter)
 -doGlf 2 : tells ANGSD to write the final genotype likelihoods into a file in beagle format
 -out all : tells ANGSD to call all output files it generate "all" followed by the appropriate suffix e.g. "all.beagle.gz"
 -P 1 : tells ANDSG use 1 thread (up to 100% CPU)
```
# Download and install NGSadmix

```bash
wget https://raw.githubusercontent.com/ANGSD/angsd/master/misc/ngsadmix32.cpp

g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
```
# Take NGSadmix output

```bash
../NDSadmix/NGSadmix -likes all.beagle.gz -K 3 -minMaf 0.05 -seed 1 -o all
```
# Run 10 reps for K values 1to 5

```bash
for run in `seq 10`;   do   mkdir run_$run ; cd run_$run/;   for K in `seq 5`;     do ../../NDSadmix/NGSadmix -likes ../all.beagle.gz -K $K -P 10 -o $K\_outfiles -minMaf 0.05;   done;   cd ../; done
```

# Load python
```bash
module load python/3.5.4
```

# Copy clumpp_input_maker_new.py to the directory. Following is clumpp_input_maker_new.py

```python
"""

This will take files from multiple runs of a particular K (from, say, Structure, or NGSAdmix), create a clumpp input file, and recommend the clumpp algorithm to choose (by calculating the D stat from the Clumpp manual).


Ex.

python3.7 clumpp_input_maker.py -in test/run*/2*qopt -type ngsadmix -out test/forClumpp

## produces test/forClumpp_clumpp.param & test/forClumpp_clumpp.in

"""


def get_cmdline_args():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-in", dest='input_files', nargs = '+',
						help = "pattern for input file name. Something unique"
						"to grab all names. should be able to search through"
						"all subdirectories holding the replicates (run_1/,"
						"run_2/, etc.): E.g., pattern: run*/2_*qopt ",
						required = True
						)
	parser.add_argument("-type", "--cluster_type",
						help = "was this and 'ngsadmix' or 'structure' run?",
						required = True
						)
	# parser.add_argument("-param", "--write_clumpp_param", help = "write a parameter file of this name for clumpp based on recommendations here.")
	parser.add_argument("-out", "--base_output_name",
						help = "base name to add to the output files for"
						 "clumpp"
						 )
	args = parser.parse_args()
	return(args)

def read_input(name):
	"""
	take a file name, read it in.
	return - list of lines
	"""
	input = list()
	IN = open(name)
	for line in IN.readlines():
		input.append(line)
	return(input)


def collect_input_data(files, filetype):
	"""
	read in the input files and collect the data
	"""
	data = list()
	for n in files:
		data.append(read_input(n))

	return(data, len(files))

def get_num_individuals(data):
	"""
	determines the number of individuals in the analysis
	"""
	return(len(data[0]))

def get_k(data):
	"""
	determines the number of clusters
	"""
	indv1 = data[0][0]
	indv1 = indv1.strip().split()# "\t"
	return(len(indv1)) ## THIS ASSUMES IT IS ONLY CLUSTER ASSIGNMENT DATA HERE

def compute_D(C, K, R, N = 1000):
	"""
	C: # individuals,
	K: # clusters,
	R: # replicates,
	N: number of replicates to test
	"""
	from math import factorial

	try: # b/c some values, the factorial will be too large for python to handle
		T_fullsearch = (factorial(K) ** (R-1)) * (R * (R-1)/2) * K * C # make sure this math works
		D_fullsearch = T_fullsearch * C * 1

		T_greedy = (factorial(K)) * (R * (R-1)/2) * K * C
		D_greedy = T_greedy * C * N
	except:
		return(3)

	# T_largeK = (R * (R − 1)/2) * (K ** 2) * C

	if D_fullsearch <= 10 ** 13:
		print("\n use fullsearch, M = 1")
		return(1)
	elif D_greedy <= 10 ** 13:
		print("\n use Greedy, M = 2, ideally a greedy_option = 1 would be set, but that may be too slow (all possible orders tested), greedy_option = 2 could be used with high R (e.g., 1000). But if K gets large, may need greedy_option = 3 and high R (e.g., 1000).")
		return(2)
	else:
		print("\n use largeKGreedy, M = 3, again, ideally greedy_option 1, but most likely it will be 2 or 3. Use large R for 2 or 3.")
		return(3)

def print_clumpp_input(data, basename = 'forClumpp', type = 'ngsadmix'):
	"""
	make input data file for clumpp
	"""
	import re
	outfile = basename + "_clumpp.in"

	## need to add the starting 5 columns to this including an individual identifier, which will not exist for NGSADMIX.
	with open(outfile, 'w') as out_fh:
		for set in data:
			# for indv in set:
			for i in range(1, len(set)+1):
				if re.search("ngsadmix", type, re.I):
					# create set of indvidual IDs and other expected columns
					line = [str(i), i, "(0)", 99, ":"]
					for l in line:
						out_fh.write("{}\t".format(l))
					out_fh.write("{}".format(set[i-1]))
					# print(set[i-1], end = '')
				else:
					print(indv, end = '')
					out_fh.write("{}".format(set[i-1]))
			out_fh.write("\n")


def print_param_file(C, K, R, M, basename):
	"""
	prints the parameter file for clumpp
	"""
	print("\nSetting GREEDY_OPTION 2, and testing 1000 repeats of input order. Unless M = 1, in which case these will be ignored.\n")

	param_outf = basename + "_clumpp.param"
	header = ["DATATYPE\t0","INDFILE\t" + basename + "_clumpp.in",
				"POPFILE\tnotneeded.popfile","OUTFILE\t" + basename + "_clumpp.out",
				"MISCFILE\tnotneeded.miscfile", "K\t" + str(K), "C\t" + str(C) , "R\t" + str(R),
				"M\t" + str(M), "W\t0", "S\t2", "GREEDY_OPTION\t2", "REPEATS\t1000","PERMUTATIONFILE\tNOTNEEDED.permutationfile",
				"PRINT_PERMUTED_DATA\t0",
				"PERMUTED_DATAFILE\tk4.perm_datafile",
				"PRINT_EVERY_PERM\t0\t",
				"EVERY_PERMFILE\tnotneeded.every_permfile",
				"PRINT_RANDOM_INPUTORDER\t0",
				"RANDOM_INPUTORDERFILE\tk4.random_inputorderfile",
				"OVERRIDE_WARNINGS\t0",
				"ORDER_BY_RUN\t0"]
	with open(param_outf, 'w') as out_fh:
		for l in header:
			out_fh.write("{}\n".format(l))

# DATATYPE 0
# INDFILE k2_forClumpp
# POPFILE notneede.popfile
# OUTFILE k2.outfile
# MISCFILE NOTNEEDED.miscfile
# K 2
# C 69
# R 15
# M 2
# W 0
# S 2
# GREEDY_OPTION 2
# REPEATS 1000
# "PERMUTATIONFILE\tNOTNEEDED.permutationfile",
# "PRINT_PERMUTED_DATA\t0",
# "PERMUTED_DATAFILE\tk4.perm_datafile",
# "PRINT_EVERY_PERM\t0\t",
# "EVERY_PERMFILE\tnotneeded.every_permfile",
# "PRINT_RANDOM_INPUTORDER\t0",
# "RANDOM_INPUTORDERFILE\tk4.random_inputorderfile",
# "OVERRIDE_WARNINGS\t0",
# "ORDER_BY_RUN\t0"





def main():
	import re


	args = get_cmdline_args()

	if not re.search("ngsadmix", args.cluster_type, re.I):
		exit("\n\nonly for ngsadmix data right now. Assumes the input cluster files only have cluster assignment on each line, no other information (Q-files).\n\n")

	cluster_data, R = collect_input_data(args.input_files, args.cluster_type)

	C = get_num_individuals(cluster_data)

	K = get_k(cluster_data)

	if re.search("ngsadmix", args.cluster_type, re.I):
		print_clumpp_input(
							cluster_data,
							args.base_output_name,
							args.cluster_type
							)

	M = compute_D(C, K, R) ## not actually doing anything about this. add that.

	print_param_file(C, K, R, M, args.base_output_name)

if __name__ == '__main__':
	main()
```
# Create CLUMPP input

```bash
for k in `seq 5` ; do python3 clumpp_input_maker_new.py -in run_*/${k}_*qopt -type ngsadmix -out k$k ; done
```

# Install CLUMPP

```bash

wget https://rosenberglab.stanford.edu/software/CLUMPP_Linux64.1.1.2.tar.gz

gunzip CLUMPP_Linux64.1.1.2.tar.gz

tar xvf CLUMPP_Linux64.1.1.2.tar
```

# Create CLUMPP file

```bash
for f in k*param ; do ../CLUMPP_Linux64.1.1.2/CLUMPP $f ; done

```

Then download these all runs folders and use it for R

#**************************************

# Make lnlks file and bamnames using following script using an excel file with relavent columns. 

```r
#make lnlks
setwd("~/Desktop/OneDrive - McMaster University/for lab and research/Tharindu on Mac/lab/structure_lplots")
require(tidyverse)

#make empty df
my_df<-setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("rep","k","lnlk"))

#list all runs
all_runs<-list.dirs(path = "./runs", full.names = TRUE, recursive = FALSE)
#list all files
all_files<-list.files(path = all_runs,recursive = TRUE,full.names = TRUE,pattern = ".log")

#make empty df
my_df<-setNames(data.frame(matrix(ncol = 3, nrow = length(all_files))), c("rep","k","lnlk"))
#add 0 and sort correctly
#for (i in 1:9) {all_files<-gsub(pattern = paste(i,"/",sep = ""),replacement = paste("0",i,"/",sep = ""),all_files)
  
#}
#all_files<-sort(all_files)



y=1

for (file in all_files) {
  

current_file<-read_lines(file = file)

#observe lnlk
lg_line<-grep(pattern = "best like=",current_file,value = TRUE) %>% str_split(pattern = "\\ ",simplify = TRUE) %>% str_split(pattern = "=") 
a<-lg_line[[2]]
lnlk<-a[[2]]

#observe k
k_val<-basename(file) %>% str_split(pattern = "\\_",simplify = TRUE)
k<-k_val[[1]]

#observe run
c<-path.expand(file)%>% str_split(pattern = "run_",simplify = TRUE) %>% str_split(pattern = "/") 
d<-c[[2]]
rep<-d[[1]]

#add data to dataframe
my_df$rep[y]=rep
my_df$k[y]=k
my_df$lnlk[y]=lnlk

y<-y+1

}

#order dataframe
ordered_df<-my_df[
  with(my_df, order(as.integer(rep), k)),
  ]

#save lnlks
write.table(ordered_df,"lnlks_allRuns.txt",sep="\t",row.names=FALSE,quote = FALSE)

#make bam filenames
require("readxl")
library("readxl")
locality_info<-read_excel("./locality_info/full_summary.xlsx")
name<-subset.data.frame(locality_info,select = 5)
location<-subset.data.frame(locality_info,select =3)
id<-subset.data.frame(locality_info,select =2)
species<-subset.data.frame(locality_info,select =4)

locality_df<-data_frame(name,location,id,species)

loc_info<-cbind(name,location,id,species)
write.table(loc_info,"bam_names.txt",sep="\t",row.names=FALSE,quote = FALSE)
```
# Example excel file

```text
IND	ID	LOCATION	SPECIES	FILE
IND1	 BJE266	vic-Lwiro	laevis	 BJE266_sorted.bam
IND2	 XGUAE_92	la-kleinmond	laevis	 XGUAE_92_sorted.bam
IND3	 XGUAE_65	la-CoGH	laevis	 XGUAE_65_sorted.bam
IND4	 BJE265	vic-Lwiro	laevis	 BJE265_sorted.bam
IND5	 XGL713_123	gi-kleinmond	gilli	 XGL713_123_sorted.bam
IND6	 XGUAE_71	la-kleinmond	laevis	 XGUAE_71_sorted.bam
IND7	 XL_CPT1	la-CoGH	laevis	 XL_CPT1_sorted.bam
IND8	 XGUAE_42	gi-CoGH	gilli	 XGUAE_42_sorted.bam
IND9	 XGUAE_59	la-CoGH	laevis	 XGUAE_59_sorted.bam
IND10	 XGL713_180	gi-kleinmond	gilli	 XGL713_180_sorted.bam
IND11	 XG92	gi-CoGH	gilli	 XG92_sorted.bam
IND12	 XGUAE_43	gi-CoGH	gilli	 XGUAE_43_sorted.bam
IND13	 BJE3545	la-Beaufort West_ Green clade	laevis	 BJE3545_sorted.bam
IND14	 BJE1488	vic-Lwiro	laevis	 BJE1488_sorted.bam
IND15	 BJE1489	vic-Lwiro	laevis	 BJE1489_sorted.bam
IND16	 BJE3540	la-Beaufort West_ Lemoenfontein Game Lodge_ Kotara Dam	laevis	 BJE3540.fqsorted.bam
IND17	 BJE3574	la-Victoria West; Jas Fontein Farm	laevis	 BJE3574.fqsorted.bam
IND18	 XGL713_177	gi-kleinmond	gilli	 XGL713_177_sorted.bam
IND19	 XGUAE_70	la-kleinmond	laevis	 XGUAE_70_sorted.bam
IND20	 BJE261	la-Bukavu	laevis	 BJE261_sorted.bam
IND21	 BJE3581	la-30 km south of Kimberley_ Debonaire farm	laevis	 BJE3581.fqsorted.bam
IND22	 BJE3536	la-Beaufort West_ Lemoenfontein Game Lodge_ Kotara Dam	laevis	 BJE3536.fqsorted.bam
IND23	 XGUAE_36	gi-CoGH	gilli	 XGUAE_36_sorted.bam
IND24	 XGUAE_72	la-kleinmond	laevis	 XGUAE_72_sorted.bam
IND25	 XGUAE_44	gi-CoGH	gilli	 XGUAE_44_sorted.bam
IND26	 XGUAE_124	la-kleinmond	laevis	 XGUAE_124_sorted.bam
IND27	 BJE267	vic-Lwiro	laevis	 BJE267_sorted.bam
IND28	 XL_CPT2	la-CoGH	laevis	 XL_CPT2_sorted.bam
IND29	 BJE3608	la-30 km south of Kimberley_ Debonaire farm	laevis	 BJE3608.fqsorted.bam
IND30	 BJE264	vic-Lwiro	laevis	 BJE264_sorted.bam
IND31	 XG153	gi-CoGH	gilli	 XG153_sorted.bam
IND32	 XGUAE_97	la-kleinmond	laevis	 XGUAE_97_sorted.bam
IND33	 BJE3639	la-Niewoudtville_ purple clade	laevis	 BJE3639_sorted.bam
IND34	 XL_CPT4	la-CoGH	laevis	 XL_CPT4_sorted.bam
IND35	 CoGH105	la-Cape of Good Hope Nature ReserveBobs hole	laevis	 CoGH105.fqsorted.bam
IND36	 XL_CPT3	la-CoGH	laevis	 XL_CPT3_sorted.bam
IND37	 XGL713_232	la-Yellow clade Krombridge	laevis	 XGL713_232_sorted.bam
IND38	 jonk_02	la-Jonkershoek	laevis	 jonk_02.fqsorted.bam
IND39	 XG12_07	gi-CoGH	gilli	 XG12_07_sorted.bam
IND40	 XGL713_181	gi-kleinmond	gilli	 XGL713_181_sorted.bam
IND41	 XLJONK_14	la-Blue clade	laevis	 XLJONK_14_sorted.bam
IND42	 BJE263	vic-Lwiro	laevis	 BJE263_sorted.bam
IND43	 JMEC006	la-Pirie Trout Hatchery 2	laevis	 JMEC006.fqsorted.bam
IND44	 XGL713_179	gi-kleinmond	gili	 XGL713_179_sorted.bam
IND45	 XGUAE_93	la-kleinmond	laevis	 XGUAE_93_sorted.bam
IND46	 JMEC003	la-Rooikrantz Dam	laevis	 JMEC003.fqsorted.bam

```

# Create Admix plots

```r
######
# plot Ks and lnlks
######

#setwd
setwd("~/Desktop/OneDrive - McMaster University/for lab and research/Tharindu on Mac/lab/structure_lplots")
#import all locations to sort
#****************************#
# you can observe sort list from here and might have to paste output for levels below

library(tidyverse)
require("readxl")
library("readxl")
library(dplyr)
library(cowplot)

full_table<-read_excel("./locality_info/full_summary.xlsx")
location_list<-(distinct(full_table[3]))
#used desc to take vic to top
sort_list<-arrange(.data = location_list,desc(LOCATION)) %>%toString %>% noquote
cat(sort_list)
#************************#

## ---- Functions -----
get_k <- function(file, k_pos = 1) {
  basename(file) %>% str_split('[_*]|K') %>%
    unlist() %>% pluck(k_pos)
}

import_names_file <- function(file){
  ids <- read_tsv(file, col_names = T)
  names(ids) <- c("indv","location","ID","species")
  ids
}

import_clumpp <- function(file, ids) {
  # need file of names in the name order as the bams, which should be same
  # as the clumpp order.
  clust <- read_table(file, col_names = F)
  clust <- read_table(file, col_names = F) %>%
    mutate(k = paste('K', ncol(clust[,c(6:ncol(clust))]), sep = ''))

  clust <- bind_cols(clust, ids)

  clust %>%
    
    gather(key = 'clust', value = 'prop', -c(X1, X2,X3,X4,X5, k, indv, location, ID,species)) %>%
    mutate(location = factor(location,
                             #change oreder here
                             levels = c("vic-Lwiro","la-Bukavu","la-Victoria West; Jas Fontein Farm", "la-Rooikrantz Dam", "la-30 km south of Kimberley_ Debonaire farm","la-Pirie Trout Hatchery 2",  "la-Beaufort West_ Lemoenfontein Game Lodge_ Kotara Dam","la-Beaufort West_ Green clade","la-Yellow clade Krombridge","la-Cape of Good Hope Nature ReserveBobs hole","la-Jonkershoek","la-Niewoudtville_ purple clade","la-Blue clade","la-kleinmond","la-CoGH", "gi-kleinmond", "gi-CoGH"),
                             labels = c("vic-Lwiro","la-Bukavu","la-Victoria West; Jas Fontein Farm", "la-Rooikrantz Dam", "la-30 km south of Kimberley_ Debonaire farm","la-Pirie Trout Hatchery 2",  "la-Beaufort West_ Lemoenfontein Game Lodge_ Kotara Dam","la-Beaufort West_ Green clade","la-Yellow clade Krombridge","la-Cape of Good Hope Nature ReserveBobs hole","la-Jonkershoek","la-Niewoudtville_ purple clade","la-Blue clade","la-kleinmond","la-CoGH", "gi-kleinmond", "gi-CoGH")
    )
    )

}
#without labs
plot_k_wtout_labs <- function(clumpp) {
  ggplot(clumpp, aes(x = indv, y = prop, fill = clust)) +
    geom_col(width = 1) +
    facet_grid( ~ location+species, switch = "x"
                , scales = "free_x"
                , space = 'free') +
    theme_bw() +
    theme(
      legend.position = 'none',
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      strip.background.x = element_blank(),
      
      #uncomment following if you want file names
      #axis.text.x = element_text(angle = 90),
      
      
      axis.text.x  = element_blank(),
      
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      #remove individual x labelling for all plots
      strip.text.x = element_blank()
    ) +
    # scale_fill_brewer(type = "qual",palette = "Paired") +
    scale_fill_manual(
      #assign different color scale for different run
      values =get(paste("pal",i,sep = "",collapse = NULL)),
      limits = names(pal)
    )+
    xlab("") +
    scale_y_continuous(breaks = c(0, 0.5, 1.0))
  
  #print(1)
 # x<-x+1
}


#final plot with labs
plot_k_wt_labs <- function(clumpp) {
  ggplot(clumpp, aes(x = indv, y = prop, fill = clust)) +
    geom_col(width = 1) +
    facet_grid( ~ location+species, switch = "x"
                , scales = "free_x"
                , space = 'free') +
    theme_bw() +
    theme(
      legend.position = 'none',
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      strip.background.x = element_blank(),
      
      #uncomment following if you want file names
      #axis.text.x = element_text(angle = 90),
      
     
      axis.text.x  = element_blank(),
      
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      #remove individual x labelling for all plots
      #strip.text.x = element_blank()
    ) +
   # scale_fill_brewer(type = "qual",palette = "Paired") +
    scale_fill_manual(
      values = pal,
      limits = names(pal)
    )+
    xlab("") +
    scale_y_continuous(breaks = c(0, 0.5, 1.0))
 
}

## --- end of functions


## -- Handle lnlk
lnlks <- read_tsv('lnlks_allRuns.txt', col_names = T)
names(lnlks) <- c("rep","k",'lnlk')

quants <-lnlks %>%
  group_by(k) %>%
  summarise(med = median(lnlk),
            upper = quantile(lnlk, 0.975),
            lower = quantile(lnlk, 0.025)
  )


lnlk_multi_plot <-
  ggplot(filter(quants, k < 6), aes(x = k, y = med, group = k)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 1) +
    geom_point(size = 5)  +
    theme_classic() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16)
    ) +
    xlab(expression(paste("Number of Clusters (", italic("K"),")"))) +
    ylab(expression(italic(log-likelihood)))


## -- handle cluster plots

sample_names <- import_names_file("bam_names.txt")
## In the same order as samples were given to NgsAdmix (list of bams order)
## E.g., (tab seperated)
# AMNH17273_Xt_SierL_male	sierra_leone	AMNH17273
# XEN091_Xt_SierL_female	sierra_leone	XEN091
# XEN092_Xt_SierL_female	sierra_leone	XEN092
# XEN094_Xt_SierL_male	sierra_leone	XEN094


c_files <- list.files("./clumpp_files/", "out", full.names = T) %>%
  set_names(., map(., get_k))
# A directory that contains files named: k4_clumpp.out, k5_clumpp.out, etc

clumpp_dats <- map(c_files, import_clumpp, ids = sample_names)

#set colors manually
#assign colour set
left_corner<-"darkblue"
right_corner<-"lightblue"
left_middle<-"forestgreen"
right_middle<-"pink"
upper_little<-"purple"

pal1 <- c(
  "X6" = right_corner,
  "X7" = left_corner, 
  "X8" = right_middle, 
  "X9" = left_middle,
  "X10"= upper_little
)

pal2 <- c(
  "X6" = right_corner,
  "X7" = left_corner, 
  "X8" = right_middle, 
  "X9" = left_middle,
  "X10"= upper_little
)

pal3 <- c(
  "X6" = left_corner,
  "X7" = right_middle, 
  "X8" = right_corner, 
  "X9" = left_middle,
  "X10"= upper_little
)



pal <- c(
  "X6" = left_corner,
  "X7" = right_corner, 
  "X8" = left_middle, 
  "X9" = right_middle,
  "X10"= upper_little
)
plot_name_list<-1:4
for (i in 1:3) {
  plot_name_list[i] <- map(clumpp_dats[i],plot_k_wtout_labs)
  i<-i+1
}

#k_plots <- map(clumpp_dats[1:3], plot_k_wtout_labs)

#last plot with labels
plot_name_list[4] <- map(clumpp_dats[4], plot_k_wt_labs)

library(gridExtra)
multi_k_plots<-grid.arrange(
  grobs = plot_name_list,
  ncol=1
  
  #widths = c(2, 2),
  #heights=c(2,1)
)
ggsave("K2-5_ngsadmix_multi.pdf", multi_k_plots,width = 30, height = 10)


## -- combine lnlk and cluster plots
join_plots <-
  cowplot::plot_grid(multi_k_plots, lnlk_multi_plot,
                   labels = "AUTO", ncol = 2
  )
xx<-ggdraw(add_sub(join_plots, "Label", vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=4.5))

ggsave("k_lnlk_plots_multi.pdf", join_plots,
        useDingbats=FALSE, width = 30, height = 10)



```

#****************END******************************

# optional


# Creating VCF with mapping quality 0f 20

following was used to get the bam file list seperated by spaces and was added in the script as file list. $1 is reference genome.
```bash
ls ../all_bam_files/ | tr "\n" " "
```
## Creating VCF

```bash

#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load bwa
module load samtools/1.10
module load nixpkgs/16.09
module load intel/2018.3
module load bcftools/1.10.2

samtools mpileup -q20 -d8000 -ugf $1 BJE1488_sorted.bam BJE1489_sorted.bam BJE261_sorted.bam BJE263_sorted.bam BJE264_sorted.bam BJE265_sorted.bam BJE266_sorted.bam BJE267_sorted.bam BJE3536.fqsorted.bam BJE3540.fqsorted.bam BJE3545_sorted.bam BJE3574.fqsorted.bam BJE3581.fqsorted.bam BJE3608.fqsorted.bam BJE3639_sorted.bam CoGH105.fqsorted.bam JMEC003.fqsorted.bam JMEC006.fqsorted.bam jonk_02.fqsorted.bam XG12_07_sorted.bam XG153_sorted.bam XG92_sorted.bam XGL713_123_sorted.bam XGL713_177_sorted.bam XGL713_179_sorted.bam XGL713_180_sorted.bam XGL713_181_sorted.bam XGL713_232_sorted.bam XGUAE_124_sorted.bam XGUAE_36_sorted.bam XGUAE_42_sorted.bam XGUAE_43_sorted.bam XGUAE_44_sorted.bam XGUAE_59_sorted.bam XGUAE_65_sorted.bam XGUAE_70_sorted.bam XGUAE_71_sorted.bam XGUAE_72_sorted.bam XGUAE_92_sorted.bam XGUAE_93_sorted.bam XGUAE_97_sorted.bam XL_CPT1_sorted.bam XL_CPT2_sorted.bam XL_CPT3_sorted.bam XL_CPT4_sorted.bam XLJONK_14_sorted.bam | bcftools call -V indels --format-fields GQ -m -O z -O z -o Xlaevis_and_gilli_all_samples_merged_sorted.bam.vcf.gz
```


