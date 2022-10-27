#############################################################
## Reverse Calling
## melodyjparker14@gmail.com - Oct 22
## GATK RNA-seq Best Practives Pipeline:
## [https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl](https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl]
#############################################################

###################################
# 0 - WORKFLOW
###################################
# - MarkDuplicates
# - SplitNCigarReads
# - BaseRecalibrator
# - ApplyBQSR
# - ScatterIntervalList
# - HaplotypeCaller
# - MergeVCFs
# - VariantFilteration

###################################
# 1 - LOAD MODULES
###################################
module load GATK/4.2.5.0-GCCcore-11.2.0-Java-11
BWA/0.7.17-GCCcore-11.2.0
module load Java/11.0.2(11)  # Java/13.0.2(13) availible but gatk uses 11 so stick with 11?
module load VCFtools/0.1.16-GCC-8.3.0
module load SAMtools/1.13-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load picard/2.23.0-Java-11

###################################
# 2 - VIEW DATA
###################################
samtools view <GTEX_filename.bam> | head -n 10
wc -l <GTEX_file.bam>
samtools view -H <GTEX_filename.bam>
# From the header of our .bam files we can see that STAR_2.5.3a was used for the alignment
# v8 uses STAR 2.5.3a, and GRCH38 v26 (https://www.gtexportal.org/home/releaseInfoPage)
# Our .bam file is, therefore, v8.
# Check for all of our other .bam files
for f in GT*.bam; do samtools view -H "$f" | if grep -q "STAR_2.5.3a"; then echo "$f"; fi ; done
# They are all v8

###################################
# 3 - SORT FILES
###################################
# Files should have already been sorted during the alignment stage. If not, use SAMtools to sort.
samtools sort file.bam -o file_sorted.bam

###################################
# 4 - INDEX BAM FILES
###################################
parallel  samtools index ::: *.bam
# If this doesn't work, make a loop
for f in *.bam; do samtools index "$f" "$f".bai ; done

###################################
# 5 - DOWNLOAD REFERENCE GENEOME 
###################################
mkdir ref_genome
cd ref_genome
#  https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md
# Ideally, use these:
wget https://storage.cloud.google.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta.fai
wget https://storage.cloud.google.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
wget https://storage.cloud.google.com/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.dict
# We do not have access. Alternatively, use these:
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
cd ..

###################################
# 6 - DOWNLOAD KNOWN SITES
###################################
mkdir known_sites
cd known_sites
# Known sites snp
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
# Known sites snp index
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
# Indels file
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
# Indels index 
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
cd ..

###################################
# 7 - START MAIN CODE
###################################
# Call genotypes on ALL postitions 
for f in GTEX*.bam
do 
	# Use Picard to mark dupicates
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT="$f" OUTPUT=dedup_"$f" M="${f%.Aligned*}"_metrics.txt
  
	# Use Picard to add or change groups (this was not needed in our case)
	# java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=dedup_"$f" O=RG_dedup_"$f" RGID=GBR RGLB=lib1 RGPL=ILLUINA SORT_ORDER=coordinate RGPU=unit1 RGSM="${f%.chrom20*}"
	
	# Use SplitNCigar to split alignment overlapping exon/intron junction and rescaled mapping quality
  	gatk SplitNCigarReads -R ref_genome/Homo_sapiens_assembly38.fasta -I dedup_"$f" -O snc_dedup_"$f"
  
  	# Use GATK to create a table for BQSR 
	gatk --java-options -Xmx4G BaseRecalibrator -R ref_genome/Homo_sapiens_assembly38.fasta -I snc_dedup_"$f" --known-sites known_sites/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz -O recal_data_"${f%.Aligned*}".table
	
  	# Use GATK to apply quality scores
	gatk --java-options -Xmx4G ApplyBQSR -R ref_genome/Homo_sapiens_assembly38.fasta -I snc_dedup_"$f" --bqsr-recal-file recal_data_"${f%.Aligned*}".table -O BQSR_snc_dedup_"$f"
	
  	# Scatter interval list
  	# Skip this step because we are not dealing with specific intervals
  	# more info: https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists
  
  	# Use GATK haplotype caller (make gvcf files)
	gatk --java-options -Xmx4G HaplotypeCaller -R ref_genome/Homo_sapiens_assembly38.fasta -I BQSR_snc_dedup_"$f" -O BQSR_snc_dedup_"${f%.bam}".g.vcf.gz -bamout BQSR_snc_dedup_"${f%.bam}".out.bam -ERC GVCF
done

# Merge VFCs
# CombineGVCFs, MergeVcfs and GatherVcfs are other options

# Create new folders
mkdir SNPS
mkdir INDELS

# Make a string including the GVCFs for GenomicsDB command
files=(*.g.vcf.gz)
printf -v joined '%s -V ' "${files[@]}"
string="${joined% -V }"

# Make a list of chromosomes then loop through them 
# I chose to include chrM (mitochondria)
bcftools query -f '%CHROM\n' <GTEX_file.g.vcf> | uniq | head -25 > chrom_list.txt
mapfile -t chrom_list < chrom_list.txt

for f in ${chrom_list[@]}
do 
	# Create a new GenomicsDB datastore from one or more GVCFs
	gatk --java-options -Xmx4G GenomicsDBImport -R ref_genome/Homo_sapiens_assembly38.fasta -V $string --genomicsdb-workspace-path GTEX_database_"$f" --intervals "$f"
	
	# Perform joint genotyping on all individual samples
	gatk --java-options -Xmx4G GenotypeGVCFs -R ref_genome/Homo_sapiens_assembly38.fasta -V gendb://GTEX_database_"$f" -O GTEX_"$f".vcf.gz
	
	# Separate snps from indels (ready for filtering)
	gatk --java-options "-Xmx4G" SelectVariants -V GTEX_"$f".vcf.gz -select-type SNP -O SNPS/GTEX_"$f"_snps.vcf.gz
	gatk --java-options "-Xmx4G" SelectVariants -V GTEX_"$f".vcf.gz -select-type INDEL -O INDELS/GTEX_"$f"_indels.vcf.gz

	# VariantFiltration
	# Label snps and indels to be filtered using gatk guidelines
	# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
	
	# Hard filtering for SNPs
	gatk VariantFiltration \
    		-V SNPS/GTEX_"$f"_snps.vcf.gz \
    		-filter "QD < 2.0" --filter-name "QD2" \
    		-filter "QUAL < 30.0" --filter-name "QUAL30" \
    		-filter "SOR > 3.0" --filter-name "SOR3" \
    		-filter "FS > 60.0" --filter-name "FS60" \
    		-filter "MQ < 40.0" --filter-name "MQ40" \
    		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    		-O SNPS/marked_filters_GTEX_"$f"_snps.vcf.gz
	
	# Hard filtering for indels
	gatk VariantFiltration \ 
    		-V INDELS/GTEX_"$f"_indels.vcf.gz \ 
    		-filter "QD < 2.0" --filter-name "QD2" \
    		-filter "QUAL < 30.0" --filter-name "QUAL30" \
    		-filter "FS > 200.0" --filter-name "FS200" \
    		-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \ 
   	 	-O INDELS/marked_filters_GTEX_"$f"_indels.vcf.gz

	# Use BCFtools to apply these filters
	# https://link.springer.com/protocol/10.1007/978-1-0716-2293-3_13#Bib1
	# https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b
	bcftools view --apply-filters .,PASS SNPS/marked_filters_GTEX_"$f"_snps.vcf.gz | bgzip -c > SNPS/filtered_GTEX_"$f"_snps.vcf.gz
	bcftools view --apply-filters .,PASS INDELS/marked_filters_GTEX_"$f"_indels.vcf.gz | bgzip -c > INDELS/filtered_GTEX_"$f"_indels.vcf.gz
done
