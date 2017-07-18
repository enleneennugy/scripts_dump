#!/bin/bash

#SBATCH --job-name=bwa-samtools
#SBATCH --account=NN9357K
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=16
#source /cluster/bin/jobsetup
module purge   # clear any inherited modules



module load samtools 
module load python2
module load bedtools 
module load bcftools
module load vcftools
module load mafft
module load emboss
module load picard-tools
module load bwa
alias seqtk="/usit/abel/u1/vincem/nobackup/seqtk/seqtk"

#DEPENDENCIES 
	#multiqc 	(Mapping)
	#qualimap 	(Mapping)
	#seqtk 		(Phase alleles)


echo -e "\n\n####------------------------------------------------------------------------------###"
echo "##--------------   Parameters   -------------------------------------------------##"
echo -e "###------------------------------------------------------------------------------####\n\n"
#change each of the parameter 

#number of processors 
PROC=16

#fasta file with the reference markers
REFERENCE=markers.fasta

#sample name list prefix. The fastq file should me name ID*_R1_*.fastq: 9_S133_L006_R1_001.fastq.gz

	#100_S224
	#101_S225
	#102_S226
	#103_S227
	#104_S228
	#105_S229
	#106_S230
	#107_S231
	#108_S232
	#109_S233
SAMPLE_NAME=list

#marker name list. 
	#269772781
	#270811630
	#270668954
	#270238663
	#270410733
	#270607648
	#270834294
	#270257485
MARKER_NAME=markers_name.list

#fullpthe where the demultiplex data are. 
WORKING_DIRECTORY=/work/users/vincem/170602_J00146.A.Project_Nielsen-Libs2-2017-05-18/0_row_data/



echo "number of cpu: "$PROC
echo "Reference file markers :" $REFERENCE
echo "Name of the samples: " $SAMPLE_NAME
echo "Name of the markers :" $MARKER_NAME
echo "Working directory :" $WORKING_DIRECTORY
	
rm -r 0_trimmed/ 1_mapping/ 2_consensus/ temp/
	
	cd $WORKING_DIRECTORY
	mkdir 0_trimmed/
	mkdir 1_mapping/
	cp $REFERENCE 0_trimmed/
	cp $SAMPLE_NAME 0_trimmed/

	### Index the reference ###
	echo -e "\n\nIndex the reference..."$SAMPLE	
	cd 0_trimmed/
	bwa index $REFERENCE -p markers


 
	for SAMPLE in $(cat $SAMPLE_NAME)
	do 
	cd $WORKING_DIRECTORY
echo -e "\n\n####------------------------------------------------------------------------------###"
echo        "##-------------- PRocessing sample "$SAMPLE"   ------------------------------------##"
echo -e "###------------------------------------------------------------------------------####"
echo -e "\n\n####------------------------------------------------------------------------------###"
echo "##--------------       Trimming   -------------------------------------------------##"
echo -e "###------------------------------------------------------------------------------####\n"

	R1=$SAMPLE"*_R1_*.fastq.gz"
	R2=$SAMPLE"*_R2_*.fastq.gz"
	R1paired=$SAMPLE"_R1_paired.fastq.gz"
	R1unpaired=$SAMPLE"_R1_unpaired.fastq.gz"
	R2paired=$SAMPLE"_R2_paired.fastq.gz"
	R2unpaired=$SAMPLE"_R2_unpaired.fastq.gz"
	java -jar /usit/abel/u1/vincem/nobackup/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads $PROC $R1 $R2  $R1paired $R1unpaired $R2paired $R2unpaired HEADCROP:10 ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:10:20 MINLEN:100

	mv $SAMPLE*_paired.fastq.gz 0_trimmed/
	mv $SAMPLE*_unpaired.fastq.gz 0_trimmed/




echo -e "\n\n####------------------------------------------------------------------------------###"
echo "##--------------       Mapping   -------------------------------------------------##"
echo -e "###------------------------------------------------------------------------------####\n\n"

	cd 0_trimmed 
	R1paired=$SAMPLE"_R1_paired.fastq.gz"
	R1unpaired=$SAMPLE"_R1_unpaired.fastq.gz"
	R2paired=$SAMPLE"_R2_paired.fastq.gz"
	R2unpaired=$SAMPLE"_R2_unpaired.fastq.gz"
	PE=$SAMPLE"_PE"
	merge=$SAMPLE"_merge.bam"
	mapped=$SAMPLE"_mapped.bam"
	
	### Mapping the trim files ###
	echo -e "\n\nMapping the trim files..."$SAMPLE
	bwa mem -t $PROC markers $R1paired $R2paired | samtools view -@ $PROC -q 20 -bT $REFERENCE - | samtools sort -@ $PROC -  > $PE.bam
	bwa mem -t $PROC markers $R1unpaired | samtools view -@ $PROC -q 20 -bT $REFERENCE - | samtools sort -@ $PROC -  > $R1unpaired.bam
	bwa mem -t $PROC markers $R2unpaired | samtools view -@ $PROC -q 20 -bT $REFERENCE - | samtools sort -@ $PROC -  > $R2unpaired.bam
	samtools merge -f -@ $PROC ../1_mapping/$mapped $PE.bam $R1unpaired.bam $R2unpaired.bam
	rm $PE.bam $R1unpaired.bam $R2unpaired.bam 
	
echo -e "\n\n####------------------------------------------------------------------------------###"
echo "##--------------  Remove duplicate ------------------------------------------------##"
echo -e "###------------------------------------------------------------------------------####\n\n"
	cd ../1_mapping/
	j=$SAMPLE"_noduplicate.bam"
	k=$SAMPLE"_metrics.txt"
	picard.jar MarkDuplicates I=$SAMPLE"_mapped.bam" O=$j REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT M=$k
	
	FOLDER=$SAMPLE"_qc"
	OUT=$SAMPLE".html"
	~/nobackup/qualimap_v2.2.1/qualimap bamqc -bam $j -outdir $SAMPLE"_qc" -outformat html -outfile $OUT -nt $PROC
	cd ..


echo -e "\n\n####------------------------------------------------------------------------------###"
echo "##--------------   Phase alleles   ------------------------------------------------##"
echo -e "###------------------------------------------------------------------------------####\n\n"
	mkdir temp
	mkdir 2_consensus


	### Phasing alleles ###
	echo -e "\n\nPhasing alleles"


	samtools phase -AF -b $SAMPLE.allele "1_mapping/"$SAMPLE"_noduplicate.bam" > log;

	### Sort phased bam files ###
	echo -e "\n\nSort phased bam files..."
	samtools sort -@ $PROC $SAMPLE.allele.0.bam -o $SAMPLE.allele.0.sorted.bam
	samtools sort -@ $PROC $SAMPLE.allele.1.bam -o $SAMPLE.allele.1.sorted.bam


	### Generating the vcf file ###
	echo -e "\n\nGenerating the vcf file..."
	samtools mpileup -uv -f $REFERENCE $SAMPLE.allele.0.sorted.bam > $SAMPLE.allele.0.vcf
	samtools mpileup -uv -f $REFERENCE $SAMPLE.allele.1.sorted.bam > $SAMPLE.allele.1.vcf

	### Generating the bed file ###
	echo -e "\n\nGenerating the bed file..."
	grep DP=0 $SAMPLE.allele.0.vcf | awk '{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}' > $SAMPLE.allele.0.bed
	grep DP=0 $SAMPLE.allele.1.vcf | awk '{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}' > $SAMPLE.allele.1.bed


	### Generating the fasta file ###
	echo -e "\n\nGenerating the fasta file..."
	bcftools call -c $SAMPLE.allele.0.vcf | vcfutils.pl vcf2fq | ~/nobackup/seqtk/seqtk seq -a - > $SAMPLE.allele.0.fasta
	bcftools call -c $SAMPLE.allele.1.vcf | vcfutils.pl vcf2fq | ~/nobackup/seqtk/seqtk seq -a - > $SAMPLE.allele.1.fasta


	### Mask indel from the fasta files ###
	echo -e "\n\nMask indel from the fasta file..."
	bedtools maskfasta -fi $SAMPLE.allele.0.fasta -fo $SAMPLE.indel.0.fasta -bed $SAMPLE.allele.0.bed
	bedtools maskfasta -fi $SAMPLE.allele.1.fasta -fo $SAMPLE.indel.1.fasta -bed $SAMPLE.allele.1.bed

	cat $SAMPLE.indel.0.fasta | sed '/^>/ s/$/_0/g'| sed "s/>/>_/g" | sed "s/>/>$SAMPLE/g" > $SAMPLE.indel.0.fasta.rename;
	cat $SAMPLE.indel.1.fasta | sed '/^>/ s/$/_1/g'| sed "s/>/>_/g" | sed "s/>/>$SAMPLE/g" > $SAMPLE.indel.1.fasta.rename; 



echo -e "\n\n####------------------------------------------------------------------------------###"
echo "##-------------- Consensus alleles ------------------------------------------------##"
echo -e "###------------------------------------------------------------------------------####\n\n"



	### Aligned alleles ###
	echo -e "\n\nAligned alleles..."
	while read GENE;
	do
	echo -e "\n\n" $GENE
	cat $SAMPLE.indel.*.fasta.rename | awk -v var=$SAMPLE"_"$GENE '$0 ~ var {print $0;getline;print $0}' | mafft --quiet --maxiterate 1000 --ep 0.4 --thread $PROC  - > temp/$SAMPLE.$GENE.fasta ;
	/cluster/software/emboss/bin/consambig -sequence temp/$SAMPLE.$GENE.fasta -outseq 2_consensus/$SAMPLE.$GENE.cons.fasta -name $SAMPLE.$GENE;
	done < $MARKER_NAME ; 

	echo -e "\n\ncleaning...\n\n"
	rm $SAMPLE*bed $SAMPLE*.rename $SAMPLE*.fasta $SAMPLE*.vcf $SAMPLE*sorted.bam log
	cd 1_mapping/
	rm $SAMPLE*mapped.bam
	rm 0_trimmed/$REFERENCE 0_trimmed/$SAMPLE_NAME 0_trimmed/markers* 

	done;
	cd ..
	pwd
	echo -e "\n\ncleaning...\n\n"
	rm -r temp/
	


	#### Build multiqc report ###
	#cd 1_mapping/
	#echo -e "\n\nBuild multiqc report..."
	#/usit/abel/u1/vincem/nobackup/MultiQC/scripts/multiqc .


echo -e "\n\n####------------------------------------------------------------------------------###"
echo "##-------------- Alignment of genes -----------------------------------------------##"
echo -e "###------------------------------------------------------------------------------####\n\n"

	mkdir 3_aligned_genes
	### Generate alignment
	echo -e "\n\nGenerate alignment...that could take a while..."

	for GENE in $(cat $MARKER_NAME); do cat 2_consensus/*."$GENE".cons.fasta | mafft --quiet --maxiterate 1000 --ep 0.4 --thread $PROC - > 3_aligned_genes/"$GENE".alignment.fasta; done;


echo -e "\n\n####------------------------------------------------------------------------------###"
echo "##--------------    Gene   Tree     -----------------------------------------------##"
echo -e "###------------------------------------------------------------------------------####\n\n"
	
	### Discard samples 
	echo -e "\n\nDiscard samples..."
	#List of samples with low capture rate

	### Fill with the missing sample
	echo -e "\n\nFill with  mithessing sample..."
	#generate empty sequences for the missing samples

	### Discard markers
	echo -e "\n\nDiscard markers..."
	#list of markers with many missing data. ratio=nb of N\-/length, discard matrice over 20

	### Rename with the species name/voucher. 

	### Make Phylip matrices 

	### Make the gene trees (array job)
	echo -e "\n\nMake the gene trees..."









