#directory paths
cwd=$(pwd)

#specify multimapper handling
# set to uniq for only uniq mapped reads
# set to phased for only phased mapped reads# set to all for all reads
multimapperHandling="uniq"


#change here your directory paths
ngsDir="${cwd}/TestData/NGS/H.sapiens"
scriptDir="${cwd}/"
genomeDir="${cwd}/TestData/Genomes/H.sapiens"
workDir="${cwd}/TestAnalysis/"
adapterFile="${cwd}/TestData/NGS/adapter.fa"

#program paths
#change here your program paths
#the workflow was tested at the following tools versions (see comments)
bbduk="java -ea -Xmx1g -cp C:/usr/local/bbmap/current/ jgi.BBDukF" 	#BBMap version 36.14
fastqc="fastqc"                                           		#v0.11.4
samtools="samtools"                                       		#1.3 (using htslib 1.3)
tRNAscanSE="tRNAscan-SE"                                  		#1.3.1 (January 2012)
bedtools="bedtools"                                       		#v2.25.0
segemehl="segemehl.x"                                     		#0.2.0-418
picardJar="picard.jar"  				 		#2.2.1
picard="java -jar $picardJar"
gatkJar="GenomeAnalysisTK.jar" 						#3.5-0-g36282e4
gatk="java -jar $gatkJar"		

# variables
#change here the names of the variable
genomeName="hg38.genomic"
tRNAName="hg38.tRNAscan"




## pre- and post-quality control, adapter and quality trimming using BBduk
cd ${ngsDir}

## adapter and quality trimming using BBduk
for i in $(ls ${ngsDir}/*.fastq.gz)
do
  bi=$(basename $i .fastq.gz)
  $fastqc -q $i
  $bbduk in=$i  out=${ngsDir}/${bi}_trimmed.fastq ref=${adapterFile} mink=8 ktrim=r k=10 rcomp=t hdist=1 qtrim=rl trimq=25 minlength=50 maxlength=100 2> ${ngsDir}/${bi}.bbduk.log
done

## pre- and post-quality control
for i in $(ls ${ngsDir}/*trimmed.fastq)
do
  bi=$(basename $i _trimmed.fastq)
  gzip $i
  $fastqc -q ${ngsDir}/${bi}_trimmed.fastq.gz
done



###genome preparation
cd ${genomeDir}

##genome indexing using samtools
$samtools faidx ${genomeName}.fa

## scan for tRNA nuclear
$tRNAscanSE -b -q -o ${tRNAName}.nuc.csv  ${genomeName}.fa

## scan for mitochondrial tRNA, consider: tRNAscanSE finds only 21 mt tRNA
cat ${genomeName}.fa | perl -lane 'BEGIN{$c=0;}if(m/^>chrM$/){$c=1}elsif(m/^>/){$c=0;}print if $c' > ${genomeName}.chrM.fa
$tRNAscanSE -b -q -O -o ${tRNAName}.chrM.csv ${genomeName}.chrM.fa

grep -v chrM ${tRNAName}.nuc.csv > ${tRNAName}.nuc_mod.csv
cat ${tRNAName}.nuc_mod.csv ${tRNAName}.chrM.csv > ${tRNAName}.csv

##convert tRNAscan tab file into bed12 entry
perl ${scriptDir}/tRNAscan2bed12.pl ${tRNAName}.csv ${tRNAName}.bed12

##mask found tRNAs genomic
$bedtools maskfasta -fi ${genomeName}.fa -fo ${genomeName}.masked.fa -mc N -bed ${tRNAName}.bed12


###create pre-tRNA library
##add 50 nt 5' and 3' flanking regions
perl ${scriptDir}/modBed12.pl ${tRNAName}.bed12 ${tRNAName}_pre-tRNAs.bed12

##remove introns, make fasta from bed12
$bedtools getfasta -name -split -s -fi ${genomeName}.fa -bed ${tRNAName}_pre-tRNAs.bed12 -fo ${tRNAName}_pre-tRNAs.fa

##add pre-tRNAs as extra chromosoms to the genome (get the artificial genome)
cat ${genomeName}.masked.fa ${tRNAName}_pre-tRNAs.fa > ${genomeName}_artificial.fa

##indexing artificial genome
$samtools faidx ${genomeName}_artificial.fa
$segemehl -x ${genomeName}_artificial.idx -d ${genomeName}_artificial.fa



###create mature tRNA library
##remove introns, make fasta from bed12
$bedtools getfasta -name -split -s -fi ${genomeName}.fa -bed ${tRNAName}.bed12 -fo ${tRNAName}.fa

##add CCA tail to tRNA chromosomes
##remove pseudogenes
perl ${scriptDir}/addCCA.pl ${tRNAName}.fa ${tRNAName}_mature.fa



###mature tRNA clustering
##only identical tRNAs were clustered
perl ${scriptDir}/clustering.pl ${tRNAName}_mature.fa ${tRNAName}_cluster.fa ${tRNAName}_clusterInfo.fa

##indexing tRNA cluster
$samtools faidx ${tRNAName}_cluster.fa
$segemehl -x ${tRNAName}_cluster.idx -d ${tRNAName}_cluster.fa
$picard CreateSequenceDictionary R=${tRNAName}_cluster.fa O=${tRNAName}_cluster.dict



###pre-mapping agains artificial genome
mkdir -p ${workDir}/mapping
cd ${workDir}/mapping

for n in $(ls ${ngsDir}/*trimmed.fastq.gz)
do
    bn=$(basename $n _trimmed.fastq.gz)

    $segemehl --silent --evalue 500 --differences 3 --maxinterval 1000 --accuracy 80 --index ${genomeDir}/${genomeName}_artificial.idx --database ${genomeDir}/${genomeName}_artificial.fa --nomatchfilename ${bn}_unmatched.fastq --query $n -o ${bn}.sam
    gzip ${bn}_unmatched.fastq

    ##remove all reads mapping at least once to the genome
    perl ${scriptDir}/removeGenomeMapper.pl ${genomeDir}/${tRNAName}_pre-tRNAs.fa ${bn}.sam ${bn}_filtered.sam

    ##remove pre-tRNA reads, keep only mature tRNA reads
    perl ${scriptDir}/removePrecursor.pl  ${genomeDir}/${tRNAName}_pre-tRNAs.bed12 ${bn}_filtered.sam $n > ${bn}_filtered.fastq
    gzip ${bn}_filtered.fastq
done


###post-processing
mkdir -p ${workDir}/postprocessing
cd ${workDir}/postprocessing

for n in $(ls ${workDir}/mapping/*filtered.fastq.gz)
do
    bn=$(basename $n _filtered.fastq.gz)

    #post-mapping against cluster
    $segemehl --silent --evalue 500 --differences 3 --maxinterval 1000 --accuracy 85 --index ${genomeDir}/${tRNAName}_cluster.idx --database ${genomeDir}/${tRNAName}_cluster.fa --nomatchfilename ${bn}_unmatched.fastq --query $n | $samtools view -bS | $samtools sort -T ${bn} -o ${bn}.bam
    gzip ${bn}_unmatched.fastq

    ##preparing bam file for indel realignment
    #indexing
    $samtools index ${bn}.bam
    $picard BuildBamIndex I=${bn}.bam O=${bn}.bai

    #add read groups to bam file
    $picard AddOrReplaceReadGroups I=${bn}.bam O=${bn}.mod.bam RGPL=RNASeqReadSimulator RGLB=Simlib RGPU=unit1 RGSM=36bam

    #indexing
    $samtools index ${bn}.mod.bam
    $picard BuildBamIndex I=${bn}.mod.bam O=${bn}.mod.bai

    #modify mapping quality to 60 (otherwise all were removed)
    $gatk -T PrintReads -R ${genomeDir}/${tRNAName}_cluster.fa -I ${bn}.mod.bam -o ${bn}.temp.bam -rf ReassignMappingQuality -DMQ 60
    mv -f ${bn}.temp.bam ${bn}.mod.bam
    rm -f ${bn}.mod.bai ${bn}.mod.bam.bai

    #indexing
    $samtools index ${bn}.mod.bam
    $picard BuildBamIndex I=${bn}.mod.bam O=${bn}.mod.bai

    ##realignment
    $gatk -T RealignerTargetCreator -R ${genomeDir}/${tRNAName}_cluster.fa -I ${bn}.mod.bam -o ${bn}.temp.intervals
    $gatk -T IndelRealigner -R ${genomeDir}/${tRNAName}_cluster.fa -I ${bn}.mod.bam -targetIntervals ${bn}.temp.intervals -o ${bn}.realigned.bam
    rm -f ${bn}.temp.intervals

    ##filter multimapped reads
    if [ "$multimapperHandling" == "uniq" ]; then
      $samtools view -h ${bn}.realigned.bam | grep -P 'NH:i:1\D'\|'^@' | $samtools view -bS | $samtools sort -T ${bn} -o ${bn}.mmHandled.bam
    elif [ "$multimapperHandling" == "phased" ]; then
      samtools sort -n -T ${bn} -O sam ${bn}.realigned.bam  > ${bn}.nSorted.sam
      perl ${scriptDir}/multimapperPhasing.pl -ed 0 -id 0 -verbose 0 -sam ${bn}.nSorted.sam -out ${bn}.nSorted.phased.sam
      samtools view -bS ${bn}.nSorted.phased.sam | samtools sort -T ${bn} -o ${bn}.mmHandled.bam
    elif [ "$multimapperHandling" == "all" ]; then
      cp ${bn}.realigned.bam ${bn}.mmHandled.bam
    else
      echo "Unkown parameter for multimapperHandling; set to 'uniq', 'all', or 'phased'";
      exit;
    fi


    #indexing
    $samtools index ${bn}.mmHandled.bam
    $picard BuildBamIndex I=${bn}.mmHandled.bam O=${bn}.mmHandled.bai

    ##modification site calling
    $gatk -R ${genomeDir}/${tRNAName}_cluster.fa -T UnifiedGenotyper -I ${bn}.mmHandled.bam -o ${bn}.GATK.vcf -stand_call_conf 50.0
    grep -i -v lowqual ${bn}.GATK.vcf > ${bn}.GATK_filtered.vcf
done
