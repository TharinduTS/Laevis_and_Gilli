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

use following to get the bam file list seperated by spaces
```bash
ls ../all_bam_files/ | tr "\n" " "
```
