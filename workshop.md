# STACKS WORKSHOP 
*29th of April 2019, Lund*



## General uppmax code
+ log in
```
ssh janneswa@rackham.uppmax.uu.se
```
+ using the nano texteditor
```
nano textfile
```
+ submit a job
```
sbatch script.sh
```
+ progress of 1 job
```
scontrol show jobid -dd 'insert jobnumber'
```
+ work on interactive node
```
interactive -A snic2017-7-126
```
+ location cleaned reads
```
/proj/sllstore2017011/ischnura_analysis/results/decloned_reads
```
+ location reference genome ischnura elegans
```
/proj/sllstore2017011/genome_10k_trimmed.fasta
```
+ location popmap
```
/proj/sllstore2017011/ischnura_analysis/results/denovo/all_popmap
```
+ example batch script
```
#!/bin/bash -l

#SBATCH -A snic2017-7-126
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -J Stacks_pipeline_denovo

module load bioinfo-tools
module load Stacks

denovo_map.pl -M 2 -n 3 -T 16 -o ./stacks_outputM2n3 --popmap ./all_popmap --samples ../decloned_reads --paired
### de novo optimisation M and n parameter
```
## Explore raw files

```
zcat sample.1.fq.gz | head
zcat sample.2.fq.gz | head
```

## Cleaning the data
+ demultiplexing

enter the barcodes into a textfile (here called isch_barcodes_lane_1)

```
#!/bin/bash -l
#SBATCH -A snic2017-7-126
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 40:00:00
#SBATCH -J cleaning_of_reads_1

module load bioinfo-tools
module load Stacks

process_radtags -P -p ../raw/lane_1_b -b ../info/isch_barcodes_lane_1 -o ../results/clean_reads -e sbfI -i gzfastq -E phred33 -r --inline_null -c -q
```

+ decloning

Decloning can be done like this

```
clone_filter -1 ../results/clean_reads/${sample}.1.fq.gz -2 ../results/clean_reads/${sample}.2.fq.gz -o ../results/decloned_reads -i gzfastq
```

or through a loop

```
#!/bin/bash -l

#SBATCH -A snic2017-7-126
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 20:00:00
#SBATCH -J declone_reads

module load bioinfo-tools
module load Stacks

#First of all, we create a variable for saving the sample identifiers
samples="sample1 sample2 sample 3"

#Then we iterate over samples to retrieve all decloned reads in the respective folder
echo ${samples} | tr " " "\n" | while read sample; do clone_filter -1 ../results/clean_reads/${sample}.1.fq.gz -2 ../results/clean_reads/${sample}.2.fq.gz -o ../results/decloned_reads -i gzfastq; echo \
${sample} processed; done
```
## Aligning data against a reference genome

+ First we need to index the reference genome
```
#SBATCH -A snic2017-7-126
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 2-00:00:00
#SBATCH -J indexing_of_reference

module load bioinfo-tools
module load bowtie2

bowtie2-build /proj/snic2017-7-126/private/assembly_50k.fasta ischnura_index
```
+ Mapping samples
```
#SBATCH -A snic2017-7-126
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 5:00:00
#SBATCH -J mapping_of_samples

module load bioinfo-tools
module load bowtie2
module load samtools

#First of all, we create a variable for saving the sample identifiers
samples="sample1 sample2 sample 3"

#Then we map the reads iterating for all this sample files and covert to the bam format
echo ${samples} | tr " " "\n" | while read sample; do bowtie2 -N 1 --mp 6,2 -I 100 -X 1000 -p 8 -x \
ischnura_index_10K_trimmed -1 ../results/decloned_reads/${sample}.1.fq.gz -2 ../results/decloned_reads/${sample}.2.fq.gz | samtools view -bS -o ../results/alignment_10K_trimmed/${sample}.bam ; echo ${sample} processed; done
```

## Running the pipeline with reference
+ make a population map file

general structure of a popmap file:

``` 
% cat popmap 
indv_01<tab>fw 
indv_02     fw 
indv_03     fw 
indv_04     oc 
indv_05     oc 
indv_06     oc
```  

+ running the pipleline

```
#!/bin/bash -l

#SBATCH -A snic2017-7-126
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 22:00:00
#SBATCH -J Stacks_pipeline_hybrid

module load bioinfo-tools
module load Stacks

ref_map.pl -T 16 --samples ../results/sorted_bam_10K_trimmed/ -o ../results/reference/stacks_output_ref --popmap ../subset/popmap_files/all
```

## Running the pipeline denovo

+ parameter optimalisation

Understanding the de novo parameters: http://catchenlab.life.illinois.edu/stacks/param_tut.php

+ -r 0.8 =loci found in 80% of the population or, r80 loci
+ popmap for circa 20 samples
+ Paris et al. 2017: We found that the highest polymorphism of r80 loci resulted from n = M, n = M or n = M + 1
+ Do for M=1 to M=6

```
denovo_map.pl -M 1 -n 1 -T 16 -o ./opt/M1_1eachpopST --popmap ./popmap_1eachpop_speciestype --samples ../decloned_reads --paired -X "populations:-r 0.8"
```

+ Count number of loci found
```
awk '{print $1}' populations.hapstats.tsv | sort | uniq | wc -l
```

+ For M values with highest number of loci, also change n to M-1 and M+1

+ running the de novo pipeline
```
denovo_map.pl -M 2 -n 3 -T 16 -o ./stacks_outputM2n3 --popmap ./all_popmap --samples ../decloned_reads --paired
```

## De novo first, then integrate with reference

as done in https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775

when (part of) samples do not reach >80% mapping rate, worthwhile comparing this method with the reference pipeline (pers. comm. Catchen, 2018)

### mapping
+ Map catalog to reference genome

bwa gives best results

```
#!/bin/bash -l

#SBATCH -A snic2017-7-126
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -J Stacks_pipeline_denovo

module load bioinfo-tools
module load Stacks
module load bwa
bwa mem -t 16 ./index/assembly_50k.fasta ./stacks_outputM2n3/catalog.fa.gz  > ./stacks_outputM2n3/catalog.sam
```

+ Statistics using SAMSTAT:
```
module load bioinfo-tools
module load samtools
module load samstat
samstat /proj/sllstore2017011/ischnura_analysis/results/denovo/stacks_outputM2n3/catalog.sam
```
+ Copy file to home directory NOTE: Do from Home terminal - not in uppmax
```
scp -r janneswa@rackham.uppmax.uu.se:/proj/sllstore2017011/ischnura_analysis/results/denovo/stacks_outputM3n3/catalog.sam.samstat.html ~/Documents/MarieCurie/project_RAD/intermediate/denovo
```

### integrate
+ sam to bam
```
module load bioinfo-tools
module load Stacks
module load python3
module load samtools

cp ./stacks_outputM2n3/catalog.sam ./bwa/catalog.sam
samtools view -Sb  ./bwa/catalog.sam  > ./bwa/catalog.bam
```
+ integrate mapping to de novo catalog
```
module load bioinfo-tools
module load Stacks
module load python3
module load samtools

stacks-integrate-alignments -P ./stacks_outputM2n3 -B ./bwa/catalog.bam -O ./integrate
```

+ copy integrate output files into de novo output folder
+ ready for populations

## Populations
+ -r =0.75    # locus has to be detected in 75% of individuals in population
+ -p=   minimum number of populations a locus must be present in to process a locus
+	--min_maf=0.05     # minimum minor allele frequency
+	--max_obs_het=0.70   # maximum accepted heterozygosity
```
populations -P ./../stacks_output_ref --popmap ./../all_popmap -O pops.maf05.r75.het70  -r 0.75 --min_maf 0.05  --max_obs_het 0.7 --genepop --plink --structure --vcf
```

+ select one random SNP per RAD tag
```
populations -P ./../stacks_output_ref --popmap ./../all_popmap -O pops.maf05.r75.het70  -r 0.75 --min_maf 0.05  --max_obs_het 0.7 --genepop --plink --structure --vcf --write_random_snp
```

## DOWNSTREAM ANALYSIS: STRUCTURE
+ take structure of random SNP per locus!

+ delete first line of file
```
sed -i '1d' populations.structure
```
+ pops from charachter to number
```
sed -i -e 's/hervias/1/g' populations.structure
sed -i -e 's/lascanas/2/g' populations.structure
sed -i -e 's/perdiguero/3/g' populations.structure
sed -i -e 's/triciovillar/4/g' populations.structure
sed -i -e 's/valbornedo/5/g' populations.structure
sed -i -e 's/xuno/6/g' populations.structure
sed -i -e 's/arles/7/g' populations.structure
sed -i -e 's/arreo/8/g' populations.structure
sed -i -e 's/belgium/9/g' populations.structure
sed -i -e 's/laigual/10/g' populations.structure
sed -i -e 's/liverpool/11/g' populations.structure
sed -i -e 's/maraixdorx/12/g' populations.structure
sed -i -e 's/menorca/13/g' populations.structure
sed -i -e 's/sweden/14/g' populations.structure
sed -i -e 's/doninos/15/g' populations.structure
sed -i -e 's/laxe/16/g' populations.structure
sed -i -e 's/louro/17/g' populations.structure
sed -i -e 's/algarve/18/g' populations.structure
sed -i -e 's/algeria/19/g' populations.structure
sed -i -e 's/cachadas/20/g' populations.structure
sed -i -e 's/sanmateo/21/g' populations.structure
sed -i -e 's/genei/22/g' populations.structure
sed -i -e 's/fountaineae/23/g' populations.structure
sed -i -e 's/saharensis/24/g' populations.structure
```

# edit mainparams
## count markers
awk -F"\t" '{print NF;exit}' populations.structure
7498-2=7496
## count inds
wc -l populations.structure
419-1=418
418/2=209

630
243

