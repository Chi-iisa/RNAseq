#!/bin/bash

#/usr/bin/time -v

#Carpeta ubicacion Fastq
Path_service='/mnt/ionpgm_storage/bioinformatics/'
FastQ_folder=$Path_service'test_exoma/XXX'

Illumina_adapter='/mnt/ionpgm_storage/bioinformatics/services_and_collaborations/adapters'
Genome_GRCH38_index='/mnt/ionpgm_storage/bioinformatics/references/GRCh38.p13.genome.fa'
Bed_File_edit='/mnt/ionpgm_storage/bioinformatics/references/Twist_Exome_Core_Covered_Targets_hg38.bed'

Genome_fa="/mnt/ionpgm_storage/bioinformatics/references/GRCh38.p13.genome.fa"
db_snp_gatk="/mnt/ionpgm_storage/bioinformatics/references/db_snp_vcf/Homo_sapiens_assembly38.dbsnp138.vcf"

#Gatk command
gatk_command="/home/bigan/miniconda3/share/gatk4-4.3.0.0-0/gatk-package-4.3.0.0-local.jar"
#Annovar command
ann_command="/mnt/ionpgm_storage/bioinformatics/annovar"
hummandb_ann="/mnt/ionpgm_storage/bioinformatics/annovar/humandb"

#DDBB VariantRecalibration
HAPMAP="/mnt/ionpgm_storage/bioinformatics/references/VariantRecal/hapmap_3.3.hg38.vcf.gz"
ONMI="/mnt/ionpgm_storage/bioinformatics/references/VariantRecal/1000G_omni2.5.hg38.vcf.gz"
PHASE="/mnt/ionpgm_storage/bioinformatics/references/VariantRecal/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
DBSNP="/mnt/ionpgm_storage/bioinformatics/references/VariantRecal/Homo_sapiens_assembly38.dbsnp138.vcf.gz"


Fastq_folder=$1

#Create temporal folders

mkdir -p $FastQ_folder"/temp"

FastQCreport_folder=$FastQ_folder"/temp/FastQCreport"

mkdir -p $FastQCreport_folder

FastQCreport_PRE=$FastQCreport_folder"/PRE"
FastQCreport_POST=$FastQCreport_folder"/POST"

mkdir -p $FastQCreport_PRE
mkdir -p $FastQCreport_POST

FastQtrimmed=$FastQ_folder"/temp/FastQtrimmed"
#FastQtrimmed=$FastQ_folder"/FastQtrimmed"

mkdir -p $FastQtrimmed

BAM_folder=$FastQ_folder"/temp/BAM"
mkdir -p $BAM_folder

VCF_folder=$FastQ_folder"/temp/VCF"
mkdir -p $VCF_folder

ANN_folder=$FastQ_folder"/temp/ANN"
mkdir -p $ANN_folder

for file in $(ls $FastQ_folder)
do
       let i=$i+1
       #Filter fastq.gz samples
       if [[ "$file" == *".fastq.gz" ]]; then
               R="$(cut -d "_" -f3 <<< $file)"

               #Get the R1 samples
               if [[ "$R" == "R1" ]]; then

                      	SampleR1=$file

                       	SampleR2=$(echo -e $(cut -d "_" -f1,2 <<<  $SampleR1)"_R2_001.fastq.gz")

                       	Samplename="$(cut -d "_" -f1 <<< $file)"
                       	#1) Quality control samples
                       	fastqc -t 12 $FastQ_folder/$SampleR1 $FastQ_folder/$SampleR2 -o $FastQCreport_PRE
                       	#2) Trimming with trimmomatic Adapters illumina and 13 beginning bases
                       	/usr/bin/time -v trimmomatic PE $FastQ_folder/$SampleR1 $FastQ_folder/$SampleR2 \
                       	$FastQtrimmed/$Samplename"_FW_PD.fq.gz" $FastQtrimmed/$Samplename"_FW_SE.fq.gz" \
                       	$FastQtrimmed/$Samplename"_RV_PD.fq.gz" $FastQtrimmed/$Samplename"_RV_SE.fq.gz" \
                       	ILLUMINACLIP:$Illumina_adapter"/TruSeq3-PE.fa":2:30:10 HEADCROP:13 \
                       	LEADING:3 TRAILING:3 MINLEN:70
                       	SLIDINGWINDOW:4:15

                       	#Merge SE and REMOVE SE files
		       	cat  $FastQtrimmed/$Samplename"_FW_SE.fq.gz" $FastQtrimmed/$Samplename"_RV_SE.fq.gz" > $FastQtrimmed/$Samplename"_MERGE_SE.fq.gz"
		       	rm $FastQtrimmed/$Samplename"_FW_SE.fq.gz" $FastQtrimmed/$Samplename"_RV_SE.fq.gz"

                       	#3)Quality control trimmed fastq
                       	fastqc -t 12 $FastQtrimmed/$Samplename"_FW_PD.fq.gz" $FastQtrimmed/$Samplename"_RV_PD.fq.gz" \
                       	$FastQtrimmed/$Samplename"_MERGE_SE.fq.gz" -o $FastQCreport_POST
				                          
                       	#4) Mapping with BWA-MEM
                       	#Map PE
                       	bwa mem -t 4 $Genome_GRCH38_index $FastQtrimmed/$Samplename"_FW_PD.fq.gz" $FastQtrimmed/$Samplename"_RV_PD.fq.gz" \
                       	 | samtools sort -O BAM -o $BAM_folder/$Samplename"_PE.bam"
                       	#Map SE
                        bwa mem -t 4 $Genome_GRCH38_index $FastQtrimmed/$Samplename"_MERGE_SE.fq.gz" | samtools sort -O BAM -o $BAM_folder/$Samplename"_SE.bam"
		      
		       	#5) Mege BAM
                       	samtools merge $BAM_folder/$Samplename"_merge.bam" $BAM_folder/$Samplename"_PE.bam" $BAM_folder/$Samplename"_SE.bam"

                      	#6) Remove duplicates
                      	samtools rmdup -S $BAM_folder/$Samplename"_merge.bam" $BAM_folder/$Samplename"_nodup.bam"

                      	#7) Select Bed panel
                      	bedtools intersect -abam $BAM_folder/$Samplename"_nodup.bam" -b $Bed_File_edit > $BAM_folder/$Samplename"_bed.bam"

		      	#8) Filter BAM by quality Q20
			#Aqui se eliminan los reads multimapeados
			samtools view -bq 20 $BAM_folder/$Samplename"_bed.bam" > $BAM_folder/$Samplename"_q20.bam"
                      	#9)AddOrReplaceReadGroups
                      	java -jar $gatk_command AddOrReplaceReadGroups I=$BAM_folder/$Samplename"_q20.bam" \
                      	O=$BAM_folder/$Samplename"replace.bam" RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
                      
                      	#9)BaseRecalibrator
                      	java -jar $gatk_command BaseRecalibrator -I $BAM_folder/$Samplename"replace.bam" -R $Genome_fa \
                        --known-sites $db_snp_gatk -O $BAM_folder/$Samplename"recal_data.table"
                      
                      	java -jar $gatk_command ApplyBQSR -R $Genome_fa -I $BAM_folder/$Samplename"replace.bam" \
                        --bqsr-recal-file $BAM_folder/$Samplename"recal_data.table" \
                        -O $BAM_folder/$Samplename"recalibrate.bam"

                      	#10) Get VCF HaplotypeCaller
                      	java -Xmx4G -jar $gatk_command HaplotypeCaller -R $Genome_fa -I $BAM_folder/$Samplename"recalibrate.bam" \
                      	-L $Bed_File_edit -O $VCF_folder/$Samplename"gatk.vcf.gz"
                      
			#11)VariantRecalibration NO HardFilter
			java -Xmx4G -jar $gatk_command VariantRecalibrator -R $Genome_fa $VCF_folder/$Samplename"gatk.vcf.gz" \ 
			--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
			--resource:omni,known=false,training=true,truth=false,prior=12.0 $ONMI \
			--resource:1000G,known=false,training=true,truth=false,prior=10.0 $PHASE \
			--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
 			-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP \
			-O $VCF_folder/$Samplename"output_snp.recal" --tranches-file $VCF_folder/$Samplename"output_snp.tranches"

			java -Xmx4G -jar $gatk_command ApplyVQSR $VCF_folder/$Samplename"gatk.vcf.gz" \
			--recal-file $VCF_folder/$Samplename"output_snp.recal" --tranches-file $VCF_folder/$Samplename"output_snp.tranches" \
			 --truth-sensitivity-filter-level 99.7 --create-output-variant-index true -mode SNP -O $VCF_folder/$Samplename"snp.recalibrated.vcf.gz"
			
			#11)Filter Variants
			gatk VariantFiltration -V $VCF_folder/$Samplename"gatk.vcf.gz" \
			--filter-expression "QD < 2.0" --filter-name QDlessthan2 \
			--filter-expression "FS > 30.0" --filter-name FSgreaterthan30 \
			--filter-expression "MQ < 30.0" --filter-name MQlessthan30 \
			--filter-expression "MQRankSum < -12.5" --filter-name MQRankSumlessthannegative12.5 \
			--filter-expression "ReadPosRankSum < -8.0" --filter-name ReadPosRankSumlessthannegative8 \
			--filter-expression "SOR > 3.0" --filter-name SORgreaterthan3 \
			-O $VCF_folder/$Samplename"gatk_filter.vcf.gz"
			
			#12) Delete BAM files
			rm $BAM_folder/$Samplename"replace.bam" $BAM_folder/$Samplename"_q20.bam" $BAM_folder/$Samplename"_bed.bam" \
			$BAM_folder/$Samplename"_nodup.bam" $BAM_folder/$Samplename"_PE.bam" $BAM_folder/$Samplename"_SE.bam"
                      	#13a) Annotation SnpEff
			snpEff -Xmx8G GRCh38.p13 $VCF_folder/$Samplename"gatk.vcf.gz" | gzip -c > $VCF_folder/$Samplename"gatk_ann.vcf.gz"
			#13b) Annotation Annovar
			#decompress vcf
			gunzip -k $VCF_folder/$Samplename"gatk.vcf.gz"
			perl $ann_command"/table_annovar.pl" $VCF_folder/$Samplename"gatk.vcf" $hummandb_ann -buildver hg38 \
			-out $ANN_folder/$Samplename -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
			-operation g,r,f,f,f -nastring . -vcfinput -polish

			#delete vcf decompress
			rm $VCF_folder/$Samplename"gatk.vcf"
               fi

       fi
done

#Filtrar anotaciones annovar por frecuencia poblacional y variantes no sinonimas
#awk -F'\t' 'NR==1 || $12<0.01 && $9 !="synonymous SNV"' 4973673.hg38_multianno.txt > 4973673.hg38_filter_syn_fq001.txt
