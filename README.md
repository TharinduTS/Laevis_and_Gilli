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
# Creating unfiltered VCF

following was used to get the bam file list seperated by spaces and was added in the script as file list. $1 is reference genome.
```bash
ls ../all_bam_files/ | tr "\n" " "
```
## Creating VCF

```bash

Inactive Modules:
  1) openmpi/2.1.1

Due to MODULEPATH changes, the following have been reloaded:
  1) bwa/0.7.17     2) samtools/1.10

The following have been reloaded with a version change:
  1) gcccore/.5.4.0 => gcccore/.7.3.0
  2) icc/.2016.4.258 => icc/.2018.3.222
  3) ifort/.2016.4.258 => ifort/.2018.3.222
  4) imkl/11.3.4.258 => imkl/2018.3.222
  5) intel/2016.4 => intel/2018.3

[warning] samtools mpileup option `u` is functional, but deprecated. Please switch to using bcftools mpileup in future.
[warning] samtools mpileup option `g` is functional, but deprecated. Please switch to using bcftools mpileup in future.
Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid
[mpileup] 1 samples in 1 input files
^C
[premacht@gra-login1 all_bam_files]$ clear

[premacht@gra-login1 all_bam_files]$ vi create_unfiltered_VCF 
[premacht@gra-login1 all_bam_files]$ bash create_unfiltered_VCF ../reference_genome/Xla.v91.repeatMasked.fa

Inactive Modules:
  1) openmpi/2.1.1

Due to MODULEPATH changes, the following have been reloaded:
  1) bwa/0.7.17     2) samtools/1.10

The following have been reloaded with a version change:
  1) gcccore/.5.4.0 => gcccore/.7.3.0
  2) icc/.2016.4.258 => icc/.2018.3.222
  3) ifort/.2016.4.258 => ifort/.2018.3.222
  4) imkl/11.3.4.258 => imkl/2018.3.222
  5) intel/2016.4 => intel/2018.3

Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid
[warning] samtools mpileup option `u` is functional, but deprecated. Please switch to using bcftools mpileup in future.
[warning] samtools mpileup option `g` is functional, but deprecated. Please switch to using bcftools mpileup in future.
[mpileup] 46 samples in 46 input files
^C
[premacht@gra-login1 all_bam_files]$ vi create_unfiltered_VCF 

#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load bwa
module load samtools/1.10
module load nixpkgs/16.09
module load intel/2018.3
module load bcftools/1.10.2

samtools mpileup -d8000 -ugf $1 BJE1488_sorted.bam BJE1489_sorted.bam BJE261_sorted.bam BJE263_sorted.bam BJE264_sorted.bam BJE265_sorted.bam BJE266_sorted.bam BJE267_sorted.bam BJE3536.fqsorted.bam BJE3540.fqsorted.bam BJE3545_sorted.bam BJE3574.fqsorted.bam BJE3581.fqsorted.bam BJE3608.fqsorted.bam BJE3639_sorted.bam CoGH105.fqsorted.bam JMEC003.fqsorted.bam JMEC006.fqsorted.bam jonk_02.fqsorted.bam XG12_07_sorted.bam XG153_sorted.bam XG92_sorted.bam XGL713_123_sorted.bam XGL713_177_sorted.bam XGL713_179_sorted.bam XGL713_180_sorted.bam XGL713_181_sorted.bam XGL713_232_sorted.bam XGUAE_124_sorted.bam XGUAE_36_sorted.bam XGUAE_42_sorted.bam XGUAE_43_sorted.bam XGUAE_44_sorted.bam XGUAE_59_sorted.bam XGUAE_65_sorted.bam XGUAE_70_sorted.bam XGUAE_71_sorted.bam XGUAE_72_sorted.bam XGUAE_92_sorted.bam XGUAE_93_sorted.bam XGUAE_97_sorted.bam XL_CPT1_sorted.bam XL_CPT2_sorted.bam XL_CPT3_sorted.bam XL_CPT4_sorted.bam XLJONK_14_sorted.bam | bcftools call -V indels --format-fields GQ -m -O z -O z -o Xlaevis_and_gilli_all_samples_merged_sorted.bam.vcf.gz

```


