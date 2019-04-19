# STACKS WORKSHOP 
29th of April 2019, Lund

## General uppmax code
+ progress of 1 job
```
scontrol show jobid -dd 'insert jobnumber'
```

## Cleaning the data



## de novo integrate with reference

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

### parameter optimalisation

+ -r 0.8 =loci found in 80% of the population or, r80 loci
+ popmap for circa 20 samples
+ Paris 2017: We found that the highest polymorphism of r80 loci resulted from n = M, n = M or n = M + 1
+ Do for M=1 to M=6

```
denovo_map.pl -M 1 -n 1 -T 16 -o ./opt/M1_1eachpopST --popmap ./popmap_1eachpop_speciestype --samples ../decloned_reads --paired -X "populations:-r 0.8"
```

+ Count number of loci found
```
awk '{print $1}' populations.hapstats.tsv | sort | uniq | wc -l
```

+ For M values with highest number of loci, also change n to M-1 and M+1

### de novo analysis
```
denovo_map.pl -M 2 -n 3 -T 16 -o ./stacks_outputM2n3 --popmap ./all_popmap --samples ../decloned_reads --paired
```
### mapping
+ Map catalog to reference genome
+ bwa gives best results
+ index is latest assembly - having scaffolds more than 50k
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

### populations
+ min_samples=0.75    # minimum per-population percentage of samples and
+	min_maf=0.05     # minimum minor allele frequency
+	max_obs_het=0.70   # maximum accepted heterozygosity
```
populations -P ./../stacks_outputM2n3 --popmap ./../all_popmap -O pops.maf05.r75.het70  -r 0.75 --min_maf 0.05  --max_obs_het 0.7 --genepop --plink --structure --vcf -f p_value --fstats -k
```
