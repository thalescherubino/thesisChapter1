#!/bin/sh

#  deNovoAssemblyDoctorate.sh
#
#
#  Created by THCRibeiro on 31/08/2020.
#

screen -r doc
###################################
#PRJEB24850 will be called "A"#####
#coffee bean ripening##############
#transcriptome of the lower canopy#
###################################

###########
#red stage#
###########

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286369

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286368 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286367 &

for seq in `ls *fastq`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta


rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *fastq -t 24

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done

rm *fastq
fastqc *gz -t 24 &

nano sampleFiles.txt

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats


###########
#yellow stage#
###########

cd ../yellowStage

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286365 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' -v --split-files ERR2286366


for seq in `ls *fastq`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *fastq -t 24 &

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done

rm *fastq

nano sampleFiles.txt


fastqc *gz -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats




#####

###########
#green stage#
###########

cd ../greenStage/

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286364

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286363

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286362


for seq in `ls *fastq`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *fastq -t 24 &

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done

rm *fastq

nano sampleFiles.txt


fastqc *gz -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats

###################################
#PRJEB24137 will be called "B"#####
#coffee bean ripening##############
#transcriptome of the upper canopy#
###################################

###########
#red stage#
###########

cd ../../
mkdir B
cd B
mkdir redStage
cd redStage

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231843

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231842

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231841


for seq in `ls *fastq`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *fastq -t 24 &

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done

rm *fastq

nano sampleFiles.txt


fastqc *gz -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats

###########
#yellow stage#
###########


~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231840

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231839

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231838



for seq in `ls *fastq`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *fastq -t 24 &


for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done

rm *fastq

nano sampleFiles.txt


fastqc *gz -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats

###########
#green stage#
###########
cd greenStage/


~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231836

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231835




for seq in `ls *fastq`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *fastq -t 24 &

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done

rm *fastq

nano sampleFiles.txt


fastqc *gz -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats


###################################
#PRJNA609253 will be called "C"#####
#Coffe leaves under differnt temperatures#
###################################

###########
#Opt temperature#
###########

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816486

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816485

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816484

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816483

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816482

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816481

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816480

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816479

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816468

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816467


for seq in `ls *fastq`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *fastq -t 24 &

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done

rm *fastq

nano sampleFiles.txt


fastqc *gz -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats


###########
#Wat temperature#
###########

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816470 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816469 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816471 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816478 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816477 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816476 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816475 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816474 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816473 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816472 &


for seq in `ls *fastq`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *fastq -t 24 &

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done

rm *fastq

nano sampleFiles.txt


fastqc *gz -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats


###################################
#PRJEB32533 will be called "D"#####
#Coffe seeds                      #
###################################


~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348595 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348589 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348587 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348582 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348581 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348577 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348574 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348573 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348571 &



for seq in `ls *fastq`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *fastq -t 24 &

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 24 $file `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *fq
done

rm *fastq

nano sampleFiles.txt


fastqc *gz -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats

###################################
#PRJEB15539 will be called "E"#####
#Coffea roots without N sources   #
###################################


~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1744633 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1739875 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1739874 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1739873 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1735800 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1735799 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1735798 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1735797 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1735796 &



for seq in `ls *fastq`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *fastq -t 24 &

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 24 $file `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *fq
done

rm *fastq

nano sampleFiles.txt


fastqc *gz -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats

###################################
#Field      will be called "F"#####
#Due to the size of the libraries##
#I will use  ACAUÂ First###########
###################################

for seq in `ls *gz`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv all.adapter.1.fasta all.adp.fasta

fastqc *fastq -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats

#Catuaí

for seq in `ls *gz`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv all.adapter.1.fasta all.adp.fasta

fastqc *fastq -t 24 &

Trinity --seqType fq --CPU 24 --full_cleanup --max_memory 140G --min_contig_length 50 --samples_file sampleFiles.txt --no_normalize_reads


/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity.stats

#remove transcripts with at leat 95% of similarity

~/bin/cdhit/cd-hit-est -i trinity_out_dir.Trinity.fasta -o filtered.Trinity.fasta -c 0.95 -n 10  -T 0 -M 0

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl filtered.Trinity.fasta > filtered.Trinity.stats


####################################
####Predict protein coding genes####
####################################

for l1 in $(ls -d *)
do
    echo $l1
    cd $l1
    for l2 in $(ls -d *)
    do
        echo $l2
        cd $l2
        #get only the longest isoform for each "gene"
    /usr/local/bioinformatic/miniconda2/opt/trinity-2.8.5/util/misc/get_longest_isoform_seq_per_trinity_gene.pl filtered.Trinity.fasta > ~/doc/$l1.$l2.transcripts.fasta
        cd ..
    done
    cd ..
done

for file in `ls *fasta`
do
    echo $file
    name=`basename -s transcripts.fasta $file`
    echo $name
    sed -i 's/TRINITY_/'"$name"'/g' $file
done

cat *fasta > all.transcriptomes.fasta

mkdir all.transcriptomes

mv all.transcriptomes.fasta all.transcriptomes

rm *fasta

cd all.transcriptomes

/usr/local/bioinformatic/miniconda2/bin/TrinityStats.pl all.transcriptomes.fasta > all.transcriptomes.stats


busco -i all.transcriptomes.fasta -l eudicots_odb10 -o all.transcriptomes.busco.result -m transcriptome -c 24


#select a balance between number of sequences and similarity
~/bin/cdhit/cd-hit-est -i all.transcriptomes.fasta -o 0.95.colapsed.all.transcriptomes.fasta -c 0.95 -n 10  -T 0 -M 0

busco -i 0.95.colapsed.all.transcriptomes.fasta -l eudicots_odb10 -o 0.95.busco.result -m transcriptome -c 24


grep Complete full_table.tsv > ../../completeExpetedGenes.tab

cd ../../

#mult line to single line in Fasta

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < 0.95.colapsed.all.transcriptomes.fasta > test.0.95.fasta

for transcript in `awk '{print $3}' completeExpetedGenes.tab`
do
echo $transcript
grep $transcript -A 1 test.0.95.fasta >> filteredBuscoCompleteTranscripts.fasta
done


busco -i filteredBuscoCompleteTranscripts.fasta -l eudicots_odb10 -o test.complete.0.95 -m transcriptome -c 24 --update-data


#we will have to run each of them to find the best trade off between uniqueness
~/bin/cdhit/cd-hit-est -i all.transcriptomes.fasta -o 0.90.colapsed.all.transcriptomes.fasta -c 0.90 -n 7  -T 0 -M 0


busco -i 0.90.colapsed.all.transcriptomes.fasta -l eudicots_odb10 -o 0.90.busco.result -m transcriptome -c 24

~/bin/cdhit/cd-hit-est -i all.transcriptomes.fasta -o 0.85.colapsed.all.transcriptomes.fasta -c 0.85 -n 6  -T 0 -M 0


busco -i 0.85.colapsed.all.transcriptomes.fasta -l eudicots_odb10 -o 0.85.busco.result -m transcriptome -c 24



~/bin/cdhit/cd-hit-est -i 0.85.colapsed.all.transcriptomes.fasta -o 0.80.colapsed.all.transcriptomes.fasta -c 0.80 -n 5 -T 0 -M 0


busco -i 0.80.colapsed.all.transcriptomes.fasta -l eudicots_odb10 -o 0.80.busco.result -m transcriptome -c 24

#we lost 3% of completedeness using .8 and reduced duplications to 38%. Lets try to keep the single copies and reduce even more the duplicated.

#separe comple from duplicated.

grep Complete 0.80.busco.result/run_eudicots_odb10/full_table.tsv > complete.0.80.tab
grep Duplicated 0.80.busco.result/run_eudicots_odb10/full_table.tsv > duplicated.0.80.tab

for transcript in `awk '{print $3}' complete.0.80.tab`
do
echo $transcript
grep $transcript -A 1 0.80.colapsed.all.transcriptomes.fasta >> 0.8.CompleteTranscripts.fasta
done


for transcript in `awk '{print $3}' duplicated.0.80.tab`
do
echo $transcript
grep $transcript -A 1 0.80.colapsed.all.transcriptomes.fasta >> 0.8.DuplicatedTranscripts.fasta
done

#~/bin/cdhit/cd-hit-est -i 0.8.DuplicatedTranscripts.fasta -o 0.8.fromDup0.8.colapsed.all.transcriptomes.fasta -c 0.8 -n 5  -T 0 -M 0

#busco -i 0.95.colapsed.all.transcriptomes.fasta -l eudicots_odb10 -o 0.95.busco.result -m transcriptome -c 24

#it donsent matter if we subset the duplicated ones, cd-hit-est will always keep the same number of sequences if we subset 0.8 from a dulicated pool of 0.8 filtered

busco -i 0.8.CompleteTranscripts.fasta -l eudicots_odb10 -o 0.8.complete.busco.result -m transcriptome -c 24

#when run busco only in complete subset we get the same result as before. 54%

#######################

#Predict protein coding genes

/usr/bin/perl ~/bin/TransDecoder/TransDecoder.LongOrfs -m 16 -t 0.80.colapsed.all.transcriptomes.fasta

/usr/bin/perl ~/bin/TransDecoder/TransDecoder.Predict --single_best_only -t 0.80.colapsed.all.transcriptomes.fasta
#--no_refine_starts


#filter only complete sequences
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < 0.80.colapsed.all.transcriptomes.fasta.transdecoder.pep > temp.fasta

grep complete -A 1 temp.fasta > 0.8.finalProteins.fasta

grep len: 0.8.finalProteins.fasta | awk '{print $5}' > len.txt

sed -i 's/len://' len.txt

#######################
 mkdir subCanephora
cd subCanephora

#/usr/local/bioinformatic/Augustus/scripts/new_species.pl --species=Ccanephora
# export PATH=$PATH:/home/tcherubino/bin/gth-1.7.1-Linux_x86_64-64bit/bin
#export BSSMDIR="$HOME/bin/gth-1.7.1-Linux_x86_64-64bit/bin/bssm"
#export GTHDATADIR="$HOME/bin/gth-1.7.1-Linux_x86_64-64bit/bin/gthdata"
/home/tcherubino/bin/BRAKER/scripts/startAlign.pl  --CPU=24 --genome=../genome/c.c.fasta --prot=../all.transcriptomes/0.8.finalProteins.fasta --prg=gth


#Convert align_gth/gth.concat.aln to GTF format:
/home/tcherubino/miniconda3/bin/gth2gtf.pl gth.concat.aln bonafide.gtf


#Compute a flanking region length for training gene structures in GenBank flatfile format

/home/tcherubino/miniconda3/bin/computeFlankingRegion.pl bonafide.gtf


#Total length gene length (including introns): 16197306. Number of genes: 8275. Average Length: 1957.37836858006
#The flanking_DNA value is: 978 (the Minimum of 10 000 and 978)

/home/tcherubino/miniconda3/bin/gff2gbSmallDNA.pl bonafide.gtf ../../genome/c.c.fasta 978 bonafide.gb


/home/tcherubino/miniconda3/bin/randomSplit.pl bonafide.gb 1600
mv bonafide.gb.test test.gb
mv bonafide.gb.train train.gb

 /home/tcherubino/miniconda3/bin/etraining --species=subCanephora train.gb &> etrain.out


#nano /home/tcherubino/miniconda3/config/species/Ccanephora/Ccanephora_parameters.cfg

#Frequency of stop codons:
#tag: 1434 (0.246)
#taa: 1965 (0.337)
#tga: 2431 (0.417)


/home/tcherubino/miniconda3/bin/augustus --species=subCanephora test.gb > test.out

grep -c LOCUS bonafide.gb
#7864
grep -c "Variable stopCodonExcludedFromCDS set right" etrain.out
#1095

################################################
#GENERATING TRAINING GENE STRUCTURES FROM ESTs##
################################################
#Alternate protocol 2

#I now believe that I should not use only the CDS because this will remove the UTRs and this will **** with the UTR prediction.... So some tinkering will be done to copy the complete EST using the IDs from complete CDSs

#filter only complete sequences
grep complete 0.80.colapsed.all.transcriptomes.fasta.transdecoder.cds | awk '{print $1}' > complente.cds.txt


awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < 0.80.colapsed.all.transcriptomes.fasta > temp.fasta


cat complente.cds.txt |  awk -F '.p' '{print $1}' > new.complente.cds.txt

grep -A 1 -F -f new.complente.cds.txt temp.fasta > 0.8.finalEST.fasta

#grep len: 0.8.finalCds.fasta | awk '{print $5}' > len.txt

#sed -i 's/len://' len.txt

#having problems to install in the Fiocruz machine
#lets transfeer the files to LFMP and run the sutuff there...


scp 0.8.finalEST.fasta lfmp@177.105.2.190:/home/lfmp


mkdir subCanephora
cd subCanephora
cp ~/0.8.finalEST.fasta ./
mv ~/c.c.fasta ./

#Clean your transcript sequences for PASA:
/home/lfmp/bin/PASApipeline.v2.4.1/bin/seqclean 0.8.finalEST.fasta

#Prepare a configuration file for running PASA. Start by copying the pasa.alignAssembly.Template.txt file from pasa_conf to your working directory:

cp /home/lfmp/bin/PASApipeline.v2.4.1/pasa_conf/pasa.alignAssembly.Template.txt  ./alignAssembly.config

#Edit the file alignAssembly.
#
# database settings
#DATABASE=ESTsubCanephora_pasa

sudo /home/lfmp/bin/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g c.c.fasta -t 0.8.finalEST.fasta.clean --ALIGNERS blat --CPU 8

#Call open reading frames (ORFs) in PASA assemblies as follows:

~/bin/PASApipeline.v2.4.1/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta ESTsubCanephora_pasa.assemblies.fasta --pasa_transcripts_gff3 ESTsubCanephora_pasa.pasa_assemblies.gff3
ls -lh


#Create a list of complete ORFs:
grep complete  ESTsubCanephora_pasa.assemblies.fasta.transdecoder.cds | sed 's/>//' | awk '{print $1}' > complete.orfs

#| sed 's/.*//'


#Find the complete ORFs in the file sample_mydb_pasa.assemblies.fasta.transdecoder.gff3, keeping only the exon and CDS entries; rename exon entries to have the same identifier as CDS entries:

grep -F -f complete.orfs  ESTsubCanephora_pasa.assemblies.fasta.transdecoder.genome.gff3 | grep -P "(\tCDS\t|\texon\t)" | perl -pe 's/cds\.//; s/\.exon\d+//; ' >  trainingSetComplete.gff3

mv trainingSetComplete.gff3 trainingSetComplete.temp.gff3

cat trainingSetComplete.temp.gff3 | perl -pe 's/\t\S*(asmbl_ \d+).*/\t$1/' | sort -n -k 4 | sort -s -k 9 | sort -s -k 1,1 > trainingSetComplete.gff3

cp trainingSetComplete.gff3 bonafide.gtf


#Follow steps 3 and 4 of Alternate Protocol 1.

/home/lfmp/bin/Augustus/scripts/computeFlankingRegion.pl bonafide.gtf


/home/lfmp/bin/Augustus/scripts/gff2gbSmallDNA.pl bonafide.gtf c.c.fasta 1323 bonafide.gb

#Now lets transfeer the bonafide.gb to Fiocruz!

scp bonafide.gb tcherubino@bagre.cpqrr.fiocruz.br:/home/tcherubino/doc/subCanephora/EST

/home/tcherubino/miniconda3/bin/randomSplit.pl bonafide.gb 1800

mv bonafide.gb.test test.gb
mv bonafide.gb.train train.gb

/home/tcherubino/miniconda3/bin/etraining --species=subCanephora train.gb &> etrain.out

/home/tcherubino/miniconda3/bin/augustus --species=subCanephora test.gb > test.out

#SUPPORT PROTOCOL 5 - UTR TRAINING
#UTRparameters are not required for gene prediction withAUGUSTUS.When integrating expression data, however, the accuracy ofAUGUSTUS may benefit fromUTRparameters (reduction of false positive CDS feature prediction in long UTRs).


#SUPPORT PROTOCOL 5 - UTR TRAINING

#I will try in the folder used to predict based on EST
cd /home/tcherubino/doc/RNAseqBasedPrediction/ccp.subgenome

mv  test.out no.UTR.out


mv bonafide.gb.test test.gb
mv bonafide.gb.train train.gb


/home/tcherubino/miniconda3/bin/optimize_augustus.pl  --cpus=24 --species=subCanephora --kfold=24 train.gb --trainOnlyUtr=1 --AUGUSTUS_CONFIG_PATH=/home/tcherubino/miniconda3/config/ --metapars=/home/tcherubino/miniconda3/config/species/subCanephora/subCanephora_metapars.utr.cfg | tee optimizeReport.UTR.out


/home/tcherubino/miniconda3/bin/augustus --species=subCanephora  test.gb --UTR=on --print_utr=on > test.utr.out

#optimize other parameters than UTR

/home/tcherubino/miniconda3/bin/optimize_augustus.pl  --cpus=24 --species=subCanephora --kfold=24 train.gb --trainOnlyUtr=0 --AUGUSTUS_CONFIG_PATH=/home/tcherubino/miniconda3/config/ | tee optimizeReport.Final.out

/home/tcherubino/miniconda3/bin/augustus --species=subCanephora  test.gb --UTR=on --print_utr=on > final.test.utr.out

/home/tcherubino/miniconda3/bin/augustus --species=subCanephora  test.gb --UTR=off > final.test.no.utr.out


################################


#a bit deeper now. Lets run protocol 1, train with RNAseq hep
#AKA basic protocol 1.

mkdir RNAseqBasedPrediction

hisat2-build canephoraChromossomes.fa canephoraChromossomes

hisat2-build eugenioidesChromossomes.fa eugenioidesChromossomes

hisat2-build unplacedContigs.fa unplacedContigs

hisat2-build GCF_003713225.1_Cara_1.0_genomic.fna GCF_003713225.1_Cara_1.0
###################################
#PRJEB24850 will be called "A"#####
#coffee bean ripening##############
#transcriptome of the lower canopy#
###################################

###########
#red stage#
###########
~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286369

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286368 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286367 &


#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta

mv rmdup.all.adp.fasta all.adp.fasta

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./rmdup.all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done
rm *fastq

for library in $( ls *1P.fq.gz )
do
#hisat2 ../genome/canephoraChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Canephora.sam --summary-file `basename -s 1P.fq.gz $library`.Canephora.alignment.report

#hisat2 ../genome/eugenioidesChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Egenioides.sam --summary-file `basename -s 1P.fq.gz $library`.Egenioides.alignment.report

#hisat2 ../genome/unplacedContigs -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.unplacedContigs.sam --summary-file `basename -s 1P.fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Cara_1.0.sam --summary-file `basename -s 1P.fq.gz $library`.Cara_1.0.alignment.report

#rm *.Cara_1.0.sam

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Cara_1.0.sam > `basename -s 1P.fq.gz $library`.Cara_1.0.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Canephora.sam > `basename -s 1P.fq.gz $library`.Canephora.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Egenioides.sam > `basename -s 1P.fq.gz $library`.Egenioides.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.unplacedContigs.sam > `basename -s 1P.fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz

###########
#yellow stage#
###########


~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286365 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' -v --split-files ERR2286366

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./rmdup.all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done
rm *fastq

for library in $( ls *1P.fq.gz )
do
#hisat2 ../genome/canephoraChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Canephora.sam --summary-file `basename -s 1P.fq.gz $library`.Canephora.alignment.report

#hisat2 ../genome/eugenioidesChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Egenioides.sam --summary-file `basename -s 1P.fq.gz $library`.Egenioides.alignment.report

#hisat2 ../genome/unplacedContigs -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.unplacedContigs.sam --summary-file `basename -s 1P.fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Cara_1.0.sam --summary-file `basename -s 1P.fq.gz $library`.Cara_1.0.alignment.report

#rm *.Cara_1.0.sam

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Cara_1.0.sam > `basename -s 1P.fq.gz $library`.Cara_1.0.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Canephora.sam > `basename -s 1P.fq.gz $library`.Canephora.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Egenioides.sam > `basename -s 1P.fq.gz $library`.Egenioides.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.unplacedContigs.sam > `basename -s 1P.fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz



###########
#green stage#
###########

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286364

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286363

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2286362


for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./rmdup.all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done
rm *fastq

for library in $( ls *1P.fq.gz )
do
#hisat2 ../genome/canephoraChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Canephora.sam --summary-file `basename -s 1P.fq.gz $library`.Canephora.alignment.report

#hisat2 ../genome/eugenioidesChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Egenioides.sam --summary-file `basename -s 1P.fq.gz $library`.Egenioides.alignment.report

#hisat2 ../genome/unplacedContigs -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.unplacedContigs.sam --summary-file `basename -s 1P.fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Cara_1.0.sam --summary-file `basename -s 1P.fq.gz $library`.Cara_1.0.alignment.report

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Cara_1.0.sam > `basename -s 1P.fq.gz $library`.Cara_1.0.bam


#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Canephora.sam > `basename -s 1P.fq.gz $library`.Canephora.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Egenioides.sam > `basename -s 1P.fq.gz $library`.Egenioides.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.unplacedContigs.sam > `basename -s 1P.fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz


###################################
#PRJEB24137 will be called "B"#####
#coffee bean ripening##############
#transcriptome of the upper canopy#
###################################

###########
#red stage#
###########

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231843

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231842

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231841


for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./rmdup.all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done
rm *fastq

for library in $( ls *1P.fq.gz )
do
#hisat2 ../genome/canephoraChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Canephora.sam --summary-file `basename -s 1P.fq.gz $library`.Canephora.alignment.report

#hisat2 ../genome/eugenioidesChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Egenioides.sam --summary-file `basename -s 1P.fq.gz $library`.Egenioides.alignment.report

#hisat2 ../genome/unplacedContigs -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.unplacedContigs.sam --summary-file `basename -s 1P.fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Cara_1.0.sam --summary-file `basename -s 1P.fq.gz $library`.Cara_1.0.alignment.report



samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Cara_1.0.sam > `basename -s 1P.fq.gz $library`.Cara_1.0.bam


#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Canephora.sam > `basename -s 1P.fq.gz $library`.Canephora.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Egenioides.sam > `basename -s 1P.fq.gz $library`.Egenioides.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.unplacedContigs.sam > `basename -s 1P.fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz


###########
#yellow stage#
###########


~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231840

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231839

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231838



for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./rmdup.all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done
rm *fastq

for library in $( ls *1P.fq.gz )
do
#hisat2 ../genome/canephoraChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Canephora.sam --summary-file `basename -s 1P.fq.gz $library`.Canephora.alignment.report

#hisat2 ../genome/eugenioidesChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Egenioides.sam --summary-file `basename -s 1P.fq.gz $library`.Egenioides.alignment.report

#hisat2 ../genome/unplacedContigs -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.unplacedContigs.sam --summary-file `basename -s 1P.fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Cara_1.0.sam --summary-file `basename -s 1P.fq.gz $library`.Cara_1.0.alignment.report



samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Cara_1.0.sam > `basename -s 1P.fq.gz $library`.Cara_1.0.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Canephora.sam > `basename -s 1P.fq.gz $library`.Canephora.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Egenioides.sam > `basename -s 1P.fq.gz $library`.Egenioides.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.unplacedContigs.sam > `basename -s 1P.fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz


###########
#green stage#
###########


~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231836

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERR2231835

for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./rmdup.all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done
rm *fastq

for library in $( ls *1P.fq.gz )
do
#hisat2 ../genome/canephoraChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Canephora.sam --summary-file `basename -s 1P.fq.gz $library`.Canephora.alignment.report

#hisat2 ../genome/eugenioidesChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Egenioides.sam --summary-file `basename -s 1P.fq.gz $library`.Egenioides.alignment.report

#hisat2 ../genome/unplacedContigs -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.unplacedContigs.sam --summary-file `basename -s 1P.fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Cara_1.0.sam --summary-file `basename -s 1P.fq.gz $library`.Cara_1.0.alignment.report



samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Cara_1.0.sam > `basename -s 1P.fq.gz $library`.Cara_1.0.bam


#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Canephora.sam > `basename -s 1P.fq.gz $library`.Canephora.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Egenioides.sam > `basename -s 1P.fq.gz $library`.Egenioides.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.unplacedContigs.sam > `basename -s 1P.fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz


###################################
#PRJNA609253 will be called "C"#####
#Coffe leaves under differnt temperatures#
###################################

###########
#Opt temperature#
###########

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816486 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816485 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816484 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816483 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816482 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816481 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816480 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816479 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816468 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816467 &


for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./rmdup.all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done
rm *fastq

for library in $( ls *1P.fq.gz )
do
#hisat2 ../genome/canephoraChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Canephora.sam --summary-file `basename -s 1P.fq.gz $library`.Canephora.alignment.report

#hisat2 ../genome/eugenioidesChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Egenioides.sam --summary-file `basename -s 1P.fq.gz $library`.Egenioides.alignment.report

#hisat2 ../genome/unplacedContigs -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.unplacedContigs.sam --summary-file `basename -s 1P.fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Cara_1.0.sam --summary-file `basename -s 1P.fq.gz $library`.Cara_1.0.alignment.report

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Cara_1.0.sam > `basename -s 1P.fq.gz $library`.Cara_1.0.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Canephora.sam > `basename -s 1P.fq.gz $library`.Canephora.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Egenioides.sam > `basename -s 1P.fq.gz $library`.Egenioides.bam

#samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.unplacedContigs.sam > `basename -s 1P.fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz

for file in `ls *bam`;
do
#samtools sort -@ 24 -n  $file >`basename $file`.sorted.bam

  java -jar ~/bin/picard.jar SortSam I=$file O=`basename -s .bam $file`.sorted.bam SO=queryname
done

for file in `ls *sorted.bam`;
do
java -XX:ParallelGCThreads=24 -jar ~/bin/picard.jar MarkDuplicates I=$file O=`basename -s sorted.bam $file`sorted.mkdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=TRUE
done


rm *sorted.bam

#samtools view -F 4 -f 2 -F 100 SRX7816486_1.trimmed_.Cara_1.0.sorted.mkdup.bam | wc -l
#samtools view -bq 1  SRX7816486_1.trimmed_.Cara_1.0.sorted.mkdup.bam | wc -l

#For a PE-reads BAM file, I would like to keep only those pairs that map to a unique location. I know I can filter reads with multiple mappings with samtools ("-bq 1"), but this doesn't take into account paired information. I mean, if mate #1 maps uniquely to locus A, and mate #2 maps both to A and B, then we should be confident that #2 maps to A. Using samtools -bq 1 removes mate #2 in this case but I would like to keep it (because of the uniqueness of mate #1). Any idea how to solve this?

for file in `ls *sorted.mkdup.bam`;
do
#-F 4 remove all unmapped
#-f 2 include all sucessfull paired
#-F 100 remove multmappers (secondary alignment)
samtools view --threads 24  -F 4 -f 2 -F 100 $file | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" {print $2,$1}' | sort -n -k1,1 > `basename -s  sorted.mkdup.bam $file`.hope
done



#Buket one to transfeer
tar -zcvf bucket1.tar.gz *bam


###########
#Wat temperature#
###########

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816470 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816469 &


~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816478 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816477 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816476 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816475 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816474 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816473 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816472 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRX7816471 &


for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./rmdup.all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *P.fq
done
rm *fastq

for library in $( ls *1P.fq.gz )
do
hisat2 ../genome/canephoraChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Canephora.sam --summary-file `basename -s 1P.fq.gz $library`.Canephora.alignment.report

hisat2 ../genome/eugenioidesChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Egenioides.sam --summary-file `basename -s 1P.fq.gz $library`.Egenioides.alignment.report

hisat2 ../genome/unplacedContigs -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.unplacedContigs.sam --summary-file `basename -s 1P.fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Cara_1.0.sam --summary-file `basename -s 1P.fq.gz $library`.Cara_1.0.alignment.report

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Cara_1.0.sam > `basename -s 1P.fq.gz $library`.Cara_1.0.bam

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Canephora.sam > `basename -s 1P.fq.gz $library`.Canephora.bam

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Egenioides.sam > `basename -s 1P.fq.gz $library`.Egenioides.bam

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.unplacedContigs.sam > `basename -s 1P.fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz

for file in `ls *bam`;
do
#samtools sort -@ 24 -n  $file >`basename $file`.sorted.bam

java -jar ~/bin/picard.jar SortSam I=$file O=`basename -s .bam $file`.sorted.bam SO=queryname
done

for file in `ls *sorted.bam`;
do
java -XX:ParallelGCThreads=24 -jar ~/bin/picard.jar MarkDuplicates I=$file O=`basename -s sorted.bam $file`sorted.mkdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics REMOVE_DUPLICATES=TRUE
done


rm *sorted.bam

#samtools view -F 4 -f 2 -F 100 SRX7816486_1.trimmed_.Cara_1.0.sorted.mkdup.bam | wc -l
#samtools view -bq 1  SRX7816486_1.trimmed_.Cara_1.0.sorted.mkdup.bam | wc -l

#For a PE-reads BAM file, I would like to keep only those pairs that map to a unique location. I know I can filter reads with multiple mappings with samtools ("-bq 1"), but this doesn't take into account paired information. I mean, if mate #1 maps uniquely to locus A, and mate #2 maps both to A and B, then we should be confident that #2 maps to A. Using samtools -bq 1 removes mate #2 in this case but I would like to keep it (because of the uniqueness of mate #1). Any idea how to solve this?

for file in `ls *sorted.mkdup.bam`;
do
#-F 4 remove all unmapped
#-f 2 include all sucessfull paired
#-F 100 remove multmappers (secondary alignment)
samtools view --threads 24  -F 4 -f 2 -F 100 $file | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" {print $2,$1}' | sort -n -k1,1 > `basename -s  sorted.mkdup.bam $file`.hope
done

rm *.fq.gz



###################################
#PRJEB32533 will be called "D"#####
#Coffe seeds                      #
###################################
#Remeber, this is single end

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348595 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348589 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348587 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348582 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348581 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348577 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348574 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348573 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX3348571 &


for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 24 $file `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *fq
done
rm *fastq

for library in $( ls *trimmed.fq.gz  )
do
hisat2 ../genome/canephoraChromossomes -U $library -p 24 -S `basename -s .fq.gz $library`.Canephora.sam --summary-file `basename -s .fq.gz $library`.Canephora.alignment.report

hisat2 ../genome/eugenioidesChromossomes -U $library -p 24 -S `basename -s .fq.gz $library`.Egenioides.sam --summary-file `basename -s .fq.gz $library`.Egenioides.alignment.report

hisat2 ../genome/unplacedContigs -U $library -p 24 -S `basename -s .fq.gz $library`.unplacedContigs.sam --summary-file `basename -s .fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -U $library -p 24 -S `basename -s .fq.gz $library`.Cara_1.0.sam --summary-file `basename -s .fq.gz $library`.Cara_1.0.alignment.report

rm *.Cara_1.0.sam

samtools view -bS -@ 24 `basename -s .fq.gz $library`.Canephora.sam > `basename -s .fq.gz $library`.Canephora.bam

samtools view -bS -@ 24 `basename -s .fq.gz $library`.Egenioides.sam > `basename -s .fq.gz $library`.Egenioides.bam

samtools view -bS -@ 24 `basename -s .fq.gz $library`.unplacedContigs.sam > `basename -s .fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz

tar -zcvf bucket2.tar.gz *bam



###################################
#PRJEB32533 will be called "E"#####
#Coffea roots without N sources   #
###################################
#remeber, this is also SINGLE END

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1744633 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1739875 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1739874 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1739873 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1735800 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1735799 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1735798 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1735797 &

~/bin/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ERX1735796 &





for file in `ls *1.fastq`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 24 $file `basename -s fastq $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:2:30:7 SLIDINGWINDOW:4:20 MINLEN:30;

rm *U.fq
pigz -p 24 *fq
done
rm *fastq

for library in $( ls *trimmed.fq.gz  )
do
hisat2 ../genome/canephoraChromossomes -U $library -p 24 -S `basename -s .fq.gz $library`.Canephora.sam --summary-file `basename -s .fq.gz $library`.Canephora.alignment.report

hisat2 ../genome/eugenioidesChromossomes -U $library -p 24 -S `basename -s .fq.gz $library`.Egenioides.sam --summary-file `basename -s .fq.gz $library`.Egenioides.alignment.report

hisat2 ../genome/unplacedContigs -U $library -p 24 -S `basename -s .fq.gz $library`.unplacedContigs.sam --summary-file `basename -s .fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -U $library -p 24 -S `basename -s .fq.gz $library`.Cara_1.0.sam --summary-file `basename -s .fq.gz $library`.Cara_1.0.alignment.report

rm *.Cara_1.0.sam

samtools view -bS -@ 24 `basename -s .fq.gz $library`.Canephora.sam > `basename -s .fq.gz $library`.Canephora.bam

samtools view -bS -@ 24 `basename -s .fq.gz $library`.Egenioides.sam > `basename -s .fq.gz $library`.Egenioides.bam

samtools view -bS -@ 24 `basename -s .fq.gz $library`.unplacedContigs.sam > `basename -s .fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz

tar -zcvf buket3.tar.gz *bam


###################################
#Yudai      will be called "F"#####
#Meristems and vascular tissue   #
###################################

#Vegetative Meristem


for seq in `ls *gz`
do
~/bin/minion search-adapter -i $seq -write-fasta `basename -s gz $seq`.adapter.minion.txt &
done

cat *adapter.minion.txt > all.adp.fasta
#remove duplicated adapters
cat ~/bin/Trimmomatic-0.39/adapters/TruSeq2*PE.fa all.adp.fasta > all.adapter.1.fasta
/home/tcherubino/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin/gt sequniq -o rmdup.all.adp.fasta all.adapter.1.fasta

rm all.adp.fasta
rm *adapter.minion.txt

mv rmdup.all.adp.fasta all.adp.fasta

fastqc *gz -t 24 &

for file in `ls *gz`;
do
echo $file
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -basein $file -baseout `basename -s gz $file`trimmed.fq ILLUMINACLIP:./all.adp.fasta:3:25:6 SLIDINGWINDOW:4:20 MINLEN:30;
#removed :1:TRUE from ILLUMINACLIP
#default 2:30:7
rm *U.fq
pigz -p 24 *P.fq
done

fastqc *P.fq.gz -t 24 &

rm *.fastq.gz

for library in $( ls *1P.fq.gz )
do
hisat2 ../genome/canephoraChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Canephora.sam --summary-file `basename -s 1P.fq.gz $library`.Canephora.alignment.report

hisat2 ../genome/eugenioidesChromossomes -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Egenioides.sam --summary-file `basename -s 1P.fq.gz $library`.Egenioides.alignment.report

hisat2 ../genome/unplacedContigs -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.unplacedContigs.sam --summary-file `basename -s 1P.fq.gz $library`.unplacedContigs.alignment.report

hisat2 ../genome/GCF_003713225.1_Cara_1.0 -1 $library -2 `basename -s 1P.fq.gz $library`2P.fq.gz -p 24 -S `basename -s 1P.fq.gz $library`.Cara_1.0.sam --summary-file `basename -s 1P.fq.gz $library`.Cara_1.0.alignment.report

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Cara_1.0.sam > `basename -s 1P.fq.gz $library`.Cara_1.0.bam

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Canephora.sam > `basename -s 1P.fq.gz $library`.Canephora.bam

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.Egenioides.sam > `basename -s 1P.fq.gz $library`.Egenioides.bam

samtools view -bS -@ 24 `basename -s 1P.fq.gz $library`.unplacedContigs.sam > `basename -s 1P.fq.gz $library`.unplacedContigs.bam

rm *sam

done

rm *.fq.gz

tar -zcvf buket4.tar.gz *bam


#########################
#As I have seen that both sub-genomes are sufficientilly different and the use of both simultaneouly as the C arabica genome increase the number of uniquelly mapped reads. So, although I have backuped all the bam files in seq@ufla.br drive I will only use Carabica for now on!

#first lets remove all unmmaped reds in the bam files, I am doing it in the sala8 desktop


#back up the header

for file in `ls *bam`;
do
    samtools view -HS $file > `basename $file`.header
done

for file in `ls *bam`;
do
#-F 4 remove all unmapped
#-f 2 include all sucessfull paired
#-F 100 remove multmappers (secondary alignment)
samtools view -F 4 -f 2 -F 100 -@ 4 $file > `basename -s sam $file`.onlyMapped.bam

done

#them let sort by queryname to remove PCR duplicates


for file in `ls *_1.0.bam`;
do
java -jar ~/bin/picard.jar SortSam I=$file O=`basename -s .bam $file`.qNameSorted.bam SO=queryname

#samtools sort -@ 24 $file >`basename -s .bam $file`.sorted.bam
done


for file in `ls *.qNameSorted.bam`;
do
java -XX:ParallelGCThreads=8 -jar ~/bin/picard.jar MarkDuplicates I=$file O=`basename -s bam $file`mkdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=`basename -s bam $file`output.metrics REMOVE_DUPLICATES=TRUE
done


for file in `ls *.mkdup.bam`;
do

samtools sort -@ 4 $file > `basename -s .mkdup.bam $file`.sorted.final.bam
done

#remove unmapped (done yet)
#keep only sucessfuly paired
#remove multmappers
#for file in `ls *sorted.final.bam`;
#do
#-F 4 remove all unmapped
#-f 2 include all sucessfull paired
#-F 100 remove multmappers (secondary alignment)
#samtools view --threads 4 -F 4 -f 2 -F 100 $file > `basename -s qNameSorted.sorted.final.bam  $file`.realyFinal.bam
#done

#rm *qNameSorted*

samtools merge -@ 4 C.a.merged.bam *sorted.final.bam

#########

for file in $(ls *Cara_1.0..hope);
do
sed -i "" '/^NC_008535.1/d' $file
sed -i "" '/^NW/d' $file
done
##########################
#prediction of C canephora subgenome
mkdir ccp.subgenome

cd ccp.subgenome

samtools view -@ 24 -b ../C.a.merged.bam NC_039898.1 NC_039900.1 NC_039902.1 NC_039904.1 NC_039907.1 NC_039908.1 NC_039910.1 NC_039913.1 NC_039914.1 NC_039917.1 NC_039918.1 > C.canephoraSubGenome.bam

samtools sort -@ 24 -n C.canephoraSubGenome.bam > name.sorted.bam

/home/tcherubino/miniconda3/bin/filterBam --uniq --paired --pairwiseAlignment --in name.sorted.bam --out Aligned.out.ssf.bam

rm name.sorted.bam

samtools sort -@ 24 Aligned.out.ssf.bam -o Aligned.out.ss.bam

/usr/local/bioinformatic/miniconda2/bin/bam2hints --intronsonly --in=Aligned.out.ss.bam --out=introns.gff
#Most RNA-seq data is by nature unstranded, i.e., it is unclear from which strand an aligned read stems. We try to guess the correct strand by using genomic splice site information (and discard all reads that do not have appropriate splice site information)

#Please check the headers of the genome FASTA file. If the headers are long and contain whitespaces, some RNA-Seq alignment tools will truncate sequence names in the BAM file. This leads to an error with BRAKER. Solution: shorten/simplify FASTA headers in the genome file before running the RNA-Seq alignment and BRAKER.

#I WILL SINPLIFI AFTER
sed "s/ .*//g" canephoraChromossomes.fa > c.c.fasta

#perl /usr/local/bioinformatic/Augustus/docs/tutorial2018/BRAKER_v2.0.4+/filterIntronsFindStrand.pl ../../genome/c.c.fasta introns.gff --score > introns.f.gff 2> error.gff

/home/tcherubino/miniconda3/bin/filterIntronsFindStrand.pl ../../genome/c.c.fasta introns.gff --score > introns.f.gff 2> error.gff


rm error.gff

~/bin/gmes_linux_64/gmes_petap.pl --verbose --sequence=../../genome/c.c.fasta --ET=introns.f.gff --cores=24 #--et_score=0 --debug

#Next, we will filter the GeneMark-ET predictions for those that have support in all introns by RNA-seq data (single-exon genes will be selected randomly in an appropriate proportion):


/home/tcherubino/miniconda3/bin/filterGenemark.pl genemark.gtf introns.f.gff


ln -sf genemark.f.good.gtf bonafide.gtf


#AUGUSTUS requires not only information about the coding sequence, but also about the non-coding sequence.We provide information about non-coding sequence in the flanking region of genes in the GenBank flatfile. The flanking region chosen should not be too short. On the other hand, setting it to maximum length may increase the runtime of training a lot. As a rule of thumb, you may use half of the mRNA size. Compute this value, e.g., with computeFlankingRegion.pl:


#perl /home/tcherubino/miniconda3/pkgs/augustus-3.3.3-pl526h0faeac2_5/bin/computeFlankingRegion.pl bonafide.gtf

/home/tcherubino/miniconda3/bin/computeFlankingRegion.pl bonafide.gtf


#Total length gene length (including introns): 18927873. Number of genes: 41019. Average Length: 461.441600234038


#The flanking_DNA value is: 230 (the Minimum of 10 000 and 230)

#Convert all GeneMark gene structures to GenBank flatfile format. Whenever two neighboring genes are separated by fewer than 2*flanking_DNA bases, the end of the flanking DNA of both is set to the midpoint between the two genes. Therefore, we do not use the filtered GeneMark gene structures here, because by using all of them we exclude potential coding regions from the flanking region of AUGUSTUS training genes, which will increase specificity of AUGUSTUS:

#perl /home/tcherubino/miniconda3/pkgs/augustus-3.3.3-pl526h0faeac2_5/bin/gff2gbSmallDNA.pl  genemark.gtf ../../genome/c.c.fasta 232 tmp.gb

#/home/tcherubino/miniconda3/bin/gff2gbSmallDNA.pl genemark.gtf ../../genome/c.c.fasta 232 tmp.gb
#/usr/local/bioinformatic/Augustus/scripts/gff2gbSmallDNA.pl genemark.gtf ../../genome/c.c.fasta 232 tmp.gb


/home/tcherubino/miniconda3/bin/gff2gbSmallDNA.pl genemark.gtf ../../genome/c.c.fasta 230 tmp.gb


#I had to do some black magic to make this work....
perl filterGenesIn_mRNAname.pl bonafide.gtf tmp.gb > bonafide.gb

ls -lh

#The file bonafide.gb can serve as an input to the protocol for eliminating redundant training genes in Support Protocol 2.

grep gene bonafide.gb | sed  's/gene="//' | sort -u > traingenes.lst

#then some nano processing


grep -w -f traingenes.lst -F bonafide.gtf > bonafide.f.gtf

#If all genes from bonafide.gtf are present in bonafide.gb, you can skip the two steps above and simply create a softlink from bonafide.gtf to bonafide.f.gtf
#pretty mutch the same the same go on!

#####################
#SUPPORT PROTOCOL 2#
###################

ln -sf bonafide.gtf bonafide.f.gtf

perl /home/tcherubino/miniconda3/bin/gtf2aa.pl ../../genome/c.c.fasta bonafide.f.gtf prot.aa

#Next, BLAST all training gene amino acid sequences against themselves and output only those protein sequences that are less than 80% redundant with any other sequence in the set:

/home/tcherubino/miniconda3/bin/aa2nonred.pl --cores=24 prot.aa prot.nr.aa

#In file prot.nr.aa, we will now find a nonredundant subset of genes. Delete the leading > characters of their names and store them in a file nonred.lst:
grep ">" prot.nr.aa | perl -pe 's/>//' > nonred.lst


cat bonafide.gb | perl -ne ' if ($_ =~ m/LOCUS\s+(\S+)\s/) { $txLocus = $1;
} elsif ($_ =~ m/\/gene=\"(\S+)\"/) { $txInGb3{$1} = $txLocus
}
if(eof()) { foreach (keys %txInGb3) { print "$_\t$txInGb3{$_}\n";
} }' > loci.lst

grep -f nonred.lst loci.lst | cut -f2 > nonred.loci.lst

#Filter bonafide.gb to keep only those loci that are contained in nonred.loci.lst:

/home/tcherubino/miniconda3/bin/filterGenesIn.pl nonred.loci.lst bonafide.gb > bonafide.f.gb


#bonafide.gb:18785
#bonafide.f.gb:15770
#84%


mv bonafide.f.gb bonafide.gb

# etraining

/home/tcherubino/miniconda3/bin/new_species.pl --species=subCanephora #--AUGUSTUS_CONFIG_PATH=/home/lfmp/bin/Augustus/config/
#/home/tcherubino/miniconda3/config/

/home/tcherubino/miniconda3/bin/etraining  --species=subCanephora bonafide.gb &> bonafide.out


#Count how many times the output file bonafide.out contains the expression Variable stopCodonExcludedFromCDS set right that is actually part of an error message

grep -c "Variable stopCodonExcludedFromCDS set right" bonafide.out

# Comparar se esse número é maior que o número de genes treinado

grep -c LOCUS bonafide.gb
#15770
#otimo, NEHNUM!


#3000 genes for test!

/home/tcherubino/miniconda3/bin/randomSplit.pl bonafide.gb 3000

mv bonafide.gb.test test.gb
mv bonafide.gb.train train.gb

/home/tcherubino/miniconda3/bin/etraining --species=subCanephora train.gb &> etrain.out

#3000 genes for test!

/home/tcherubino/miniconda3/bin/augustus --species=subCanephora  test.gb > test.out


#Frequency of stop codons:
#tag: 4645 (0.294)
#taa: 4570 (0.29)
#tga: 6555 (0.416)



#BASIC PROTOCOL 3  - PREDICTION USING EXTRINSIC EVIDENCE

#This section describes the structural genome annotation process when experimental evidence is used to identify (parts of) gene structures, to uncover alternative splicing, or to overall improve annotation quality

#Each hint can be given a group name, by specifying group=groupname;or grp=groupname; in the last column for the hint in the GFF file. In the case of alignments of longer sequences like IsoSeq (Alternate Protocol 6), proteins (Alternate Protocol 7), or ESTs, this should be used to group together all the hints coming from the alignment of the same sequence to the genome. For

#for the integration of transcriptome data into gene prediction. Coverage information of RNA-seq data will not only cover coding regions, but also UTRs. The absence of UTR parameters during gene prediction with such data may lead to false positive coding sequence prediction in longer UTRs because coverage exceeds the coding sequence and stretches into UTRs.

#Generating Hints From RNA-Seq Data RNA-seq

#For generation of intron hints from spliced read information, please follow the first six steps in Basic Protocol 1. The resulting file introns.gff contains hints from spliced RNA-seq reads.

#Begin with the file Aligned.out.ss.bam from step 4 or step 5 in Basic Protocol 1. Convert this file to wiggle format:


/home/tcherubino/miniconda3/bin/bam2wig Aligned.out.ss.bam > rnaseq.wig

cat rnaseq.wig | /home/tcherubino/miniconda3/bin/wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --UCSC=unstranded.track --radius=4.5 --pri=4 --strand="."  > rnaseq.gff

#When using this type of hints, remember to call AUGUSTUS with the option --UTR=on

#ALTERNATE PROTOCOL 7

#GENERATING HINTS FROM PROTEIN DATA

#Suitable hints for AUGUSTUS can be generated by various protein-to-genome alignment tools. In general, there are tools that produce alignments with a high specificity at the expense of sensitivity (e.g., GenomeThreader; Gremme, 2013) and more sensitive but less specific tools (e.g., Exonerate; Slater & Birney, 2005). The optimal choice of alignment tools depends on the number of proteins, on the degree of relatedness of proteins and genome, and on the available computational resources (e.g., Exonerate is slow in comparison to other aligners, and will require data parallelization).


#Generally, we recommend using an aligner with a high specificity rather than an aligner
#with a high sensitivity, in particular if you are aligning more than one proteome to the
#target genome.

#THIS WAS DONE previously! SO JUST GO TO THE RESPECTIVE PROTEIN FOLDER
/home/tcherubino/bin/BRAKER/scripts/startAlign.pl  --CPU=24 --genome=../genome/c.c.fasta --prot=../all.transcriptomes/0.8.finalProteins.fasta --prg=gth

/home/tcherubino/bin/BRAKER/scripts/align2hints.pl --in=gth.concat.aln --out=prot.hints --prg=gth

cp prot.hints ~/doc/RNAseqBasedPrediction/ccp.subgenome/

#ALTERNATE PROTOCOL 8
#GENERATING HINTS FROM ESTs AND MEDIUM SIZE RNA-SEQ READS ESTs

#STs and medium size RNA-seq reads (>ß450 nt) are suitable for generating intron, exonpart (sometimes abbreviated as “ep”), and exon hints.

cd ../EST

cp /home/tcherubino/doc/all.transcriptomes/0.80.colapsed.all.transcriptomes.fasta.transdecoder.cds ./

#Align ESTs/cDNAs/medium sized RNA-seq reads to genome using BLAT:
/home/tcherubino/bin/blat -noHead -minIdentity=92 ../../genome/c.c.fasta 0.80.colapsed.all.transcriptomes.fasta.transdecoder.cds cdna.psl

#Filter alignments to obtain those that are potentially most useful for gene prediction:
/home/tcherubino/bin/pslCDnaFilter -minId=0.9 -localNearBest=0.005 -ignoreNs -bestOverlap cdna.psl cdna.f.psl

#Sort filtered alignments according to start position and target sequence name

cat cdna.f.psl | sort -n -k 16,16 | sort -s -k 14,14 > cdna.fs.psl

#Convert psl-file to hints: bash

/home/tcherubino/miniconda3/bin/blat2hints.pl --in=cdna.fs.psl --out=cdna.hints 20

#BASIC PROTOCOL 4
#RUNNING AUGUSTUS WITH HINTS

#Concatenate all hints. If you have created more than one hints file, e.g., hints.proteins.gff from Alternate Protocol 7 and hints.rnaseq.gff from Basic Protocol 3, then join them in a single file:

mkdir finalPrediction
cd finalPrediction


cp ../align_gth/prot.hints ./
cp ../EST/cdna.hints ./

cp ../../RNAseqBasedPrediction/ccp.subgenome/rnaseq.gff ./


cp ../../RNAseqBasedPrediction/ccp.subgenome/introns.f.gff ./

#and change the priorities

sed -i 's/pri=4/pri=3/g' rnaseq.gff
sed -i 's/pri=4/pri=5/g' prot.hints

cat prot.hints cdna.hints rnaseq.gff introns.f.gff > hints.gff

rm prot.hints cdna.hints rnaseq.gff introns.f.gff


#Hints from different sources can still be distinguished by augustus via the src attribute in the last column, e.g., with src=P and src=E for hints from proteins and RNA-seq, respectively.

cp /home/tcherubino/miniconda3/config/extrinsic/extrinsic.M.RM.E.W.P.cfg ./

/home/tcherubino/miniconda3/bin/augustus --species=subCanephora --UTR=on --extrinsicCfgFile=extrinsic.M.RM.E.W.P.cfg  --allow_hinted_splicesites=atac ../../genome/c.c.fasta --codingseq=on --outfile=predicted.C.canephora.gff --progress=true --genemodel=complete --hintsfile=hints.gff

/home/tcherubino/miniconda3/bin/getAnnoFasta.pl --seqfile ../../genome/c.c.fasta predicted.C.canephora.gff

sed 's/>/>SubC.c_/g' predicted.C.canephora.aa > subCanephoraProteins.aa

busco -f -i subCanephoraProteins.aa -l eudicots_odb10 -o buscoProtResult -m prot -c 24 --long

awk '{if($2 == "Complete" || $2 == "Duplicated") {print $0}}' full_table.tsv > complete.C.c.busco.tab

#remove unecessary data from GFF.

sed '/Delete/d' predicted.C.canephora.gff > Temp

sed -i"" '/Error/d' Temp

sed -i"" '/b2h/d' Temp

mv Temp processed.C.canephora.gff

sed -i"" 's/"g/"SubC.c_/g' processed.C.canephora.gff


###########################
#Busco NCBI annotation ####
###########################

cd /home/tcherubino/doc/busco.95.similarity.NCBI.proteins



busco -f -i GCF_003713225.1_Cara_1.0_protein.faa -l eudicots_odb10 -o buscoProtResult -m prot -c 24 --long


##########################
########Some stats#######
#######################
R

chormossomes <- read.table("chrm.table.tab")

#aca.opt.1
#aca.opt.1 <- list()

#aca.opt.1[[1]] <- read.table("SRX7816467_1.trimmed_.Cara_1.0..hope")

#aca.opt.1[[2]] <- read.table("SRX7816467_1.trimmed_.Canephora..hope")

#aca.opt.1[[3]] <- read.table("SRX7816467_1.trimmed_.Egenioides..hope")

#aca.opt.1[[4]] <- 5329566


#arabica

#(sum(aca.opt.1[[1]]$V2))/ aca.opt.1[[4]]

#Canephora

#(sum(aca.opt.1[[2]]$V2))/ aca.opt.1[[4]]

#Eugenioides

#(sum(aca.opt.1[[3]]$V2))/ aca.opt.1[[4]]

#aca.opt.2
aca.opt.2 <- list()

aca.opt.2[[1]] <- read.table("SRX7816468_1.trimmed_.Cara_1.0..hope")

aca.opt.2[[2]] <- read.table("SRX7816468_1.trimmed_.Canephora..hope")

aca.opt.2[[3]] <- read.table("SRX7816468_1.trimmed_.Egenioides..hope")

t <-read.table("SRX7816468_1.trimmed_.unplacedContigs..hope")

aca.opt.2[[4]] <- 8679307


#arabica

(sum(aca.opt.2[[1]]$V2))/ aca.opt.2[[4]]

#Canephora

(sum(aca.opt.2[[2]]$V2))/ aca.opt.2[[4]]

#Eugenioides

(sum(aca.opt.2[[3]]$V2))/ aca.opt.2[[4]]


#aca.opt.3
aca.opt.3 <- list()

aca.opt.3[[1]] <- read.table("SRX7816479_1.trimmed_.Cara_1.0..hope")

aca.opt.3[[2]] <- read.table("SRX7816479_1.trimmed_.Canephora..hope")

aca.opt.3[[3]] <- read.table("SRX7816479_1.trimmed_.Egenioides..hope")

aca.opt.3[[4]] <- 11443394


#arabica

(sum(aca.opt.3[[1]]$V2))/ aca.opt.3[[4]]

#Canephora

(sum(aca.opt.3[[2]]$V2))/ aca.opt.3[[4]]

#Eugenioides

(sum(aca.opt.3[[3]]$V2))/ aca.opt.3[[4]]


#aca.opt.4
aca.opt.4 <- list()

aca.opt.4[[1]] <- read.table("SRX7816480_1.trimmed_.Cara_1.0..hope")

aca.opt.4[[2]] <- read.table("SRX7816480_1.trimmed_.Canephora..hope")

aca.opt.4[[3]] <- read.table("SRX7816480_1.trimmed_.Egenioides..hope")

aca.opt.4[[4]] <- 5983576


#arabica

(sum(aca.opt.4[[1]]$V2))/ aca.opt.4[[4]])

#Canephora

(sum(aca.opt.4[[2]]$V2))/ aca.opt.4[[4]]

#Eugenioides

(sum(aca.opt.4[[3]]$V2))/ aca.opt.4[[4]]


#aca.opt.5
aca.opt.5 <- list()

aca.opt.5[[1]] <- read.table("SRX7816481_1.trimmed_.Cara_1.0..hope")

aca.opt.5[[2]] <- read.table("SRX7816481_1.trimmed_.Canephora..hope")

aca.opt.5[[3]] <- read.table("SRX7816481_1.trimmed_.Egenioides..hope")

aca.opt.5[[4]] <- 8795369


#arabica

(sum(aca.opt.5[[1]]$V2))/ aca.opt.5[[4]]

#Canephora

(sum(aca.opt.5[[2]]$V2))/ aca.opt.5[[4]]

#Eugenioides

(sum(aca.opt.5[[3]]$V2))/ aca.opt.5[[4]]



#cat.opt.1
cat.opt.1 <- list()

cat.opt.1[[1]] <- read.table("SRX7816482_1.trimmed_.Cara_1.0..hope")

cat.opt.1[[2]] <- read.table("SRX7816482_1.trimmed_.Canephora..hope")

cat.opt.1[[3]] <- read.table("SRX7816482_1.trimmed_.Egenioides..hope")

cat.opt.1[[4]] <- 6274944


#arabica

(sum(cat.opt.1[[1]]$V2))/ cat.opt.1[[4]]

#Canephora

(sum(cat.opt.1[[2]]$V2))/ cat.opt.1[[4]]

#Eugenioides

(sum(cat.opt.1[[3]]$V2))/ cat.opt.1[[4]]


#cat.opt.2
cat.opt.2 <- list()

cat.opt.2[[1]] <- read.table("SRX7816483_1.trimmed_.Cara_1.0..hope")

cat.opt.2[[2]] <- read.table("SRX7816483_1.trimmed_.Canephora..hope")

cat.opt.2[[3]] <- read.table("SRX7816483_1.trimmed_.Egenioides..hope")

cat.opt.2[[4]] <- 5735496


#arabica

(sum(cat.opt.2[[1]]$V2))/ cat.opt.2[[4]]

#Canephora

(sum(cat.opt.2[[2]]$V2))/ cat.opt.2[[4]]

#Eugenioides

(sum(cat.opt.2[[3]]$V2))/ cat.opt.2[[4]]


#cat.opt.3
cat.opt.3 <- list()

cat.opt.3[[1]] <- read.table("SRX7816484_1.trimmed_.Cara_1.0..hope")

cat.opt.3[[2]] <- read.table("SRX7816484_1.trimmed_.Canephora..hope")

cat.opt.3[[3]] <- read.table("SRX7816484_1.trimmed_.Egenioides..hope")

cat.opt.3[[4]] <- 8159130


#arabica

(sum(cat.opt.3[[1]]$V2))/ cat.opt.3[[4]]

#Canephora

(sum(cat.opt.3[[2]]$V2))/ cat.opt.3[[4]]

#Eugenioides

(sum(cat.opt.3[[3]]$V2))/ cat.opt.3[[4]]
########################

#cat.opt.4
cat.opt.4 <- list()

cat.opt.4[[1]] <- read.table("SRX7816485_1.trimmed_.Cara_1.0..hope")

cat.opt.4[[2]] <- read.table("SRX7816485_1.trimmed_.Canephora..hope")

cat.opt.4[[3]] <- read.table("SRX7816485_1.trimmed_.Egenioides..hope")

cat.opt.4[[4]] <- 9165147


#arabica

(sum(cat.opt.4[[1]]$V2))/ cat.opt.4[[4]]

#Canephora

(sum(cat.opt.4[[2]]$V2))/ cat.opt.4[[4]]

#Eugenioides

(sum(cat.opt.4[[3]]$V2))/ cat.opt.4[[4]]



#cat.opt.5
cat.opt.5 <- list()

cat.opt.5[[1]] <- read.table("SRX7816486_1.trimmed_.Cara_1.0..hope")

cat.opt.5[[2]] <- read.table("SRX7816486_1.trimmed_.Canephora..hope")

cat.opt.5[[3]] <- read.table("SRX7816486_1.trimmed_.Egenioides..hope")

cat.opt.5[[4]] <- 9792147


#arabica

(sum(cat.opt.5[[1]]$V2))/ cat.opt.5[[4]]

#Canephora

(sum(cat.opt.5[[2]]$V2))/ cat.opt.5[[4]]

#Eugenioides

(sum(cat.opt.5[[3]]$V2))/ cat.opt.5[[4]]

#aca.wat.1
aca.wat.1 <- list()

aca.wat.1[[1]] <- read.table("SRX7816469_1.trimmed_.Cara_1.0..hope")

aca.wat.1[[2]] <- read.table("SRX7816469_1.trimmed_.Canephora..hope")

aca.wat.1[[3]] <- read.table("SRX7816469_1.trimmed_.Egenioides..hope")

aca.wat.1[[4]] <- 9660156


#arabica

(sum(aca.wat.1[[1]]$V2))/ aca.wat.1[[4]]

#Canephora

(sum(aca.wat.1[[2]]$V2))/ aca.wat.1[[4]]

#Eugenioides

(sum(aca.wat.1[[3]]$V2))/ aca.wat.1[[4]]

#aca.wat.2
aca.wat.2 <- list()

aca.wat.2[[1]] <- read.table("SRX7816470_1.trimmed_.Cara_1.0..hope")

aca.wat.2[[2]] <- read.table("SRX7816470_1.trimmed_.Canephora..hope")

aca.wat.2[[3]] <- read.table("SRX7816470_1.trimmed_.Egenioides..hope")

aca.wat.2[[4]] <- 8407959


#arabica

(sum(aca.wat.2[[1]]$V2))/ aca.wat.2[[4]]

#Canephora

(sum(aca.wat.2[[2]]$V2))/ aca.wat.2[[4]]

#Eugenioides

(sum(aca.wat.2[[3]]$V2))/ aca.wat.2[[4]]



#aca.wat.3
aca.wat.3 <- list()

aca.wat.3[[1]] <- read.table("SRX7816472_1.trimmed_.Cara_1.0..hope")

aca.wat.3[[2]] <- read.table("SRX7816472_1.trimmed_.Canephora..hope")

aca.wat.3[[3]] <- read.table("SRX7816472_1.trimmed_.Egenioides..hope")

aca.wat.3[[4]] <- 9779528


#arabica

(sum(aca.wat.3[[1]]$V2))/ aca.wat.3[[4]]

#Canephora

(sum(aca.wat.3[[2]]$V2))/ aca.wat.3[[4]]

#Eugenioides

(sum(aca.wat.3[[3]]$V2))/ aca.wat.3[[4]]


#aca.wat.4
aca.wat.4 <- list()

aca.wat.4[[1]] <- read.table("SRX7816473_1.trimmed_.Cara_1.0..hope")

aca.wat.4[[2]] <- read.table("SRX7816473_1.trimmed_.Canephora..hope")

aca.wat.4[[3]] <- read.table("SRX7816473_1.trimmed_.Egenioides..hope")

aca.wat.4[[4]] <- 10272079


#arabica

(sum(aca.wat.4[[1]]$V2))/ aca.wat.4[[4]]

#Canephora

(sum(aca.wat.4[[2]]$V2))/ aca.wat.4[[4]]

#Eugenioides

(sum(aca.wat.4[[3]]$V2))/ aca.wat.4[[4]]



#cat.wat.1
cat.wat.1 <- list()

cat.wat.1[[1]] <- read.table("SRX7816474_1.trimmed_.Cara_1.0..hope")

cat.wat.1[[2]] <- read.table("SRX7816474_1.trimmed_.Canephora..hope")

cat.wat.1[[3]] <- read.table("SRX7816474_1.trimmed_.Egenioides..hope")

cat.wat.1[[4]] <- 6712081


#arabica

(sum(cat.wat.1[[1]]$V2))/ cat.wat.1[[4]]

#Canephora

(sum(cat.wat.1[[2]]$V2))/ cat.wat.1[[4]]

#Eugenioides

(sum(cat.wat.1[[3]]$V2))/ cat.wat.1[[4]]

#cat.wat.2
#cat.wat.2 <- list()

#cat.wat.2[[1]] <- read.table("SRX7816475_1.trimmed_.Cara_1.0..hope")

#cat.wat.2[[2]] <- read.table("SRX7816475_1.trimmed_.Canephora..hope")

#cat.wat.2[[3]] <- read.table("SRX7816475_1.trimmed_.Egenioides..hope")

#cat.wat.2[[4]] <- 9665306


#arabica

#(sum(cat.wat.2[[1]]$V2))/ cat.wat.2[[4]]

#Canephora

#(sum(cat.wat.2[[2]]$V2))/ cat.wat.2[[4]]

#Eugenioides

#(sum(cat.wat.2[[3]]$V2))/ cat.wat.2[[4]]

#cat.wat.3
cat.wat.3 <- list()

cat.wat.3[[1]] <- read.table("SRX7816476_1.trimmed_.Cara_1.0..hope")

cat.wat.3[[2]] <- read.table("SRX7816476_1.trimmed_.Canephora..hope")

cat.wat.3[[3]] <- read.table("SRX7816476_1.trimmed_.Egenioides..hope")

cat.wat.3[[4]] <- 8747562


#arabica

(sum(cat.wat.3[[1]]$V2))/ cat.wat.3[[4]]

#Canephora

(sum(cat.wat.3[[2]]$V2))/ cat.wat.3[[4]]

#Eugenioides

(sum(cat.wat.3[[3]]$V2))/ cat.wat.3[[4]]


#cat.wat.4
cat.wat.4 <- list()

cat.wat.4[[1]] <- read.table("SRX7816477_1.trimmed_.Cara_1.0..hope")

cat.wat.4[[2]] <- read.table("SRX7816477_1.trimmed_.Canephora..hope")

cat.wat.4[[3]] <- read.table("SRX7816477_1.trimmed_.Egenioides..hope")

cat.wat.4[[4]] <- 9965254


#arabica

(sum(cat.wat.4[[1]]$V2))/ cat.wat.4[[4]]

#Canephora

(sum(cat.wat.4[[2]]$V2))/ cat.wat.4[[4]]

#Eugenioides

(sum(cat.wat.4[[3]]$V2))/ cat.wat.4[[4]]

#cat.wat.5
cat.wat.5 <- list()

cat.wat.5[[1]] <- read.table("SRX7816478_1.trimmed_.Cara_1.0..hope")

cat.wat.5[[2]] <- read.table("SRX7816478_1.trimmed_.Canephora..hope")

cat.wat.5[[3]] <- read.table("SRX7816478_1.trimmed_.Egenioides..hope")

cat.wat.5[[4]] <- 12646761


#arabica

(sum(cat.wat.5[[1]]$V2))/ cat.wat.5[[4]]

#Canephora

(sum(cat.wat.5[[2]]$V2))/ cat.wat.5[[4]]

#Eugenioides

(sum(cat.wat.5[[3]]$V2))/ cat.wat.5[[4]]

########
###old##
########


for file in $( ls *Canephora.alignment.report);
do
grep "aligned concordantly exactly 1 time" $file >> Canephora.one.time
grep "overall alignment rate" $file >> Canephora.overall
done

#Cara_1.0.alignment.report

for file in $( ls *Cara_1.0.alignment.report);
do
grep "aligned concordantly exactly 1 time" $file >> Cara.one.time
grep "overall alignment rate" $file >> Cara.overall
done

#Egenioides.alignment.report

for file in $( ls *Egenioides.alignment.report);
do
grep "aligned concordantly exactly 1 time" $file >> Egenioides.one.time
grep "overall alignment rate" $file >> Egenioides.overall
done

#unplacedContigs.alignment.report


for file in $( ls *unplacedContigs.alignment.report);
do
grep "aligned concordantly exactly 1 time" $file >> unplacedContigs.one.time
grep "overall alignment rate" $file >> unplacedContigs.overall
done


sed -i" " 's/\% overall alignment rate//g' Canephora.overall
sed -i" " 's/\% overall alignment rate//g' Cara.overall
sed -i" " 's/\% overall alignment rate//g' Egenioides.overall
sed -i" " 's/\% overall alignment rate//g'unplacedContigs.overall


paste Cara.overall Canephora.overall Egenioides.overall unplacedContigs.overall > all.overall

cat temp all.overall > all.o

#oneTime

awk '{print $2}' unplacedContigs.one.time | sed 's/\%//g' > unplacedContigs.ot
awk '{print $2}' Cara.one.time | sed 's/\%//g' > Cara.ot
awk '{print $2}' Egenioides.one.time | sed 's/\%//g' > Egenioides.ot
awk '{print $2}' Canephora.one.time | sed  's/\%//g' > Canephora.ot

#had to remove ( and ) with nano


paste Cara.ot Canephora.ot Egenioides.ot unplacedContigs.ot > all.one.time


C.arabicaOneTime        C.canephoraOneTime      C.eugenioidesOneTime    unplacedContigsOneTime
cat temp2 all.one.time > all.ot

R

one.time <- read.table("all.one.time",header=T)
overall <- read.table("all.overall",header=T)

one.timeMeans <- colMeans(one.time)
overallMeans <-  colMeans(overall)

overallMeans <- overallMeans - one.timeMeans

testeMatrix <- cbind(one.timeMeans, overallMeans)

rownames(testeMatrix) <- c("C. arabica","sub-C.canephora","sub-C.eugenioides", "unplaced")

#cols <- c("#2c76ab", "#3791d1", "#442a9c", "#5d39d4", "#2a9456","#3bcc75","#52db3a" ,"#40b02f")

oneTimeCIara <- as.matrix(confint(lm(one.time$Cara.ot ~ 1)))
oneTimeCIcan <- as.matrix(confint(lm(one.time$Canephora.ot ~ 1)))
oneTimeCIeug <- as.matrix(confint(lm(one.time$Egenioides.ot ~ 1)))
oneTimeCIunp <- as.matrix(confint(lm(one.time$unplacedContigs.ot ~ 1)))


OverAllCIara <- as.matrix(confint(lm(overall$C.arabicaOneTime ~ 1)))
OverAllCIcan <- as.matrix(confint(lm(overall$C.canephoraOneTime ~ 1)))
OverAllCIeug <- as.matrix(confint(lm(overall$C.eugenioidesOneTime ~ 1)))
OverAllCIunp <- as.matrix(confint(lm(overall$unplacedContigsOneTime ~ 1)))



pdf("overall.alignmet.pdf",h=4.5,w=8.5)
barplot <- barplot(t(testeMatrix), ylim=c(0,100),cex.names=1, col = c("#522c28","#c24a3a"),ylab="% of mapped reads",main="Overall aligment to genome",space=0.5)

arrows(barplot, c(oneTimeCIara[,1], oneTimeCIcan[,1], oneTimeCIeug[,1], oneTimeCIunp[,1]), barplot, y1=c(oneTimeCIara[,2], oneTimeCIara[,2], oneTimeCIeug[,2], oneTimeCIunp[,2]),code = 3, angle=90, lwd=2,col= "#7bc437")



arrows(barplot, c(OverAllCIara[,1], OverAllCIcan[,1], OverAllCIeug[,1], OverAllCIunp[,1]), barplot, y1=c(OverAllCIara[,2], OverAllCIcan[,2], OverAllCIeug[,2], OverAllCIunp[,2]),code = 3, angle=90, lwd=2,col = "#7bc437" )



legend("topright",c("uniquely mapped", "mult mapped"), col=c("#522c28","#c24a3a"),pch=15)
dev.off()



limMax = 10

nice.matrix1 <- matrix(nrow=17,ncol=22)

nice.matrix2 <- matrix(nrow=17,ncol=11)

nice.matrix3<- matrix(nrow=17,ncol=11)
i=1

for (whateevah in ls()){
print(whateevah)
print(class(get(whateevah)))
    if(is.list(get(whateevah))) {
        print("ok")
        pdf(paste0(whateevah,".pdf"),w=7,h=7)
        par(mfrow=c(3,1))

whateevah <- "aca.opt.2"

        temp1 <- get(whateevah)[[1]]$V2/get(whateevah)[[4]]*100

        barplot <- barplot(get(whateevah)[[1]]$V2/get(whateevah)[[4]]*100,ylim=c(0,limMax),col =c("#c24a3a","#522c28"))

        title(main=paste(whateevah,"aligned to complete C. arabica genome"),ylab ="% of uniquely mapped reads",xlab=)

        axis(side= 1,at=barplot,labels=chormossomes$V1[match(get(whateevah)[[1]]$V1, chormossomes$V2)],las=2,xpd=T )

        temp2 <- get(whateevah)[[2]]$V2/get(whateevah)[[4]]*100

        barplot <- barplot(get(whateevah)[[2]]$V2/get(whateevah)[[4]]*100,ylim=c(0,limMax),col =c("#c24a3a"))

        title(main=paste(whateevah,"aligned to C. canephora sub-genome"),ylab ="% of uniquely mapped reads")


        axis(side= 1,at=barplot,labels=chormossomes$V1[match(get(whateevah)[[2]]$V1, chormossomes$V2)],las=2,xpd=T)


        temp3 <- get(whateevah)[[3]]$V2/get(whateevah)[[4]]*100

        barplot <- barplot(get(whateevah)[[3]]$V2/get(whateevah)[[4]]*100,ylim=c(0,limMax),col =c("#522c28"))

        title(main=paste(whateevah,"aligned to C. eugenioides sub-genome"),ylab ="% of uniquely mapped reads")


        axis(side= 1,at=barplot,labels=chormossomes$V1[match(get(whateevah)[[3]]$V1, chormossomes$V2)],las=2,xpd=T )

        dev.off()

        nice.matrix1[i,] <- temp1

        nice.matrix2[i,] <- temp2

        nice.matrix3[i,] <- temp3

        i<- i+1

        }

}


nice.matrix1.ci <- matrix(nrow=2,ncol=22)

for(i in 1:22){

nice.matrix1.ci[,i] <-  confint(lm(nice.matrix1[,i] ~ 1))

}

nice.matrix2.ci <- matrix(nrow=2,ncol=11)

for(i in 1:11){

nice.matrix2.ci[,i] <-  confint(lm(nice.matrix2[,i] ~ 1))

}

nice.matrix3.ci <- matrix(nrow=2,ncol=11)

for(i in 1:11){

nice.matrix3.ci[,i] <-  confint(lm(nice.matrix3[,i] ~ 1))

}


####Plot all chromossomes with CI
pdf("FINALFigure2.pdf",w=7,h=7)


par(mfrow=c(3,1))

 barplot <- barplot(colMeans(nice.matrix1),ylim=c(0,limMax),col =c("#c24a3a","#522c28"))

title(main="aligned to complete C. arabica genome",ylab ="% of uniquely mapped reads",xlab=)


arrows(barplot, nice.matrix1.ci[1,], barplot, y1=nice.matrix1.ci[2,],code = 3, angle=90, lwd=1.5,length=0.05,col = "#7bc437" )

axis(side= 1,at=barplot,labels=c("Chr_1c","Chr_1e","Chr_2c","Chr_2e","Chr_3c","Chr_3e","Chr_4c","Chr_4e","Chr_5e","Chr_5c","Chr_6c","Chr_6e","Chr_7c","Chr_7e","Chr_8e","Chr_8c","Chr_9c","Chr_9e","Chr_10e","Chr_10c","Chr_11c","Chr_11e") ,las=2,xpd=T )


barplot <- barplot(colMeans(nice.matrix2),ylim=c(0,limMax),col =c("#c24a3a","#522c28"))

title(main="aligned to C. canephora sub-genome",ylab ="% of uniquely mapped reads")

arrows(barplot, nice.matrix2.ci[1,], barplot, y1=nice.matrix2.ci[2,],code = 3, angle=90, lwd=1.5,length=0.1,col = "#7bc437" )

axis(side= 1,at=barplot,labels=c("Chr_1c","Chr_2c","Chr_3c","Chr_4c","Chr_5c","Chr_6c","Chr_7c","Chr_8c","Chr_9c","Chr_10c","Chr_11c") ,las=2,xpd=T )


barplot <- barplot(colMeans(nice.matrix3),ylim=c(0,limMax),col =c("#c24a3a","#522c28"))

title("aligned to C. eugenioides sub-genome",ylab ="% of uniquely mapped reads")

arrows(barplot, nice.matrix3.ci[1,], barplot, y1=nice.matrix3.ci[2,],code = 3, angle=90, lwd=1.5,length=0.1,col = "#7bc437" )

axis(side= 1,at=barplot,labels=c("Chr_1e","Chr_2e","Chr_3e","Chr_4e","Chr_5e","Chr_6e","Chr_7e","Chr_8e","Chr_9e","Chr_10e","Chr_11e") ,las=2,xpd=T )
dev.off()


aca.opt.2
aca.opt.3
aca.opt.4
aca.opt.5
aca.wat.1
aca.wat.2
aca.wat.3
aca.wat.4
cat.opt.1
cat.opt.2
cat.opt.3
cat.opt.4
cat.opt.5
cat.wat.1
cat.wat.3
cat.wat.4
cat.wat.5

ChromossomeUsage <- read.table("ChromossomeUsage.txt", header= T)

ChromossomeUsage$MultComp <- paste0(ChromossomeUsage$chromossomeNumber,".",ChromossomeUsage$chromossomeHomolog)

model <- lm(mapped ~ 0 + ChromossomeUsage$MultComp + cultivar + temperature ,data = ChromossomeUsage)



aov(model)

t <- TukeyHSD(aov(model))

model2<- lm(mapped ~ 0 + chromossomeHomolog+ cultivar + temperature  ,data = ChromossomeUsage)

aov(model2)


TukeyHSD(aov(model2))
#library("lme4")
#library("lmerTest")
#model2 <- lmer(mapped ~ 0 + cultivar + temperature + chromossomeNumber * chromossomeHomolog + (1|library),data = ChromossomeUsage)

library(multcomp)

contr<- as.matrix(read.delim("contrast.table.tab",row.names=1))

resp.contr <- glht(model, linfct = contr)

summary(resp.contr)

##########################
#length Stats##############
##########################


Arabica	 <- read.delim("arabicaLengths.tab")
Canephora	 <- read.delim("canephoraLengths.tab")
Eugenioides <- read.delim("ceugenioidesLengths.tab")


#BUSCO GO
#C.canephora

busco -f -i subCanephoraProteins.aa -l eudicots_odb10 -o C.c.buscoProtResult -m prot -c 24 --long --config  ~/bin/busco-4.1.4/config/config.ini


grep -w -f exclusive.BUSCO.C.canephora.tab complete.C.c.busco.tab | awk '{print $3}' > C.c.exclusiveIDS.txt

#retrieve their Descriptions

grep -w -f C.c.exclusiveIDS.txt C.canephora.GenDescriptions.tab> C.c.exclusiveIDS.Descriptions.tab

#retrieve their GO terms

grep -w -f C.c.exclusiveIDS.txt GO.onePerRow.subcanephora.tab | awk '{if( $2 != "") print $0}' > C.c.exclusiveIDS.GO.tab

awk '{if($2 == "Fragmented") {print $0}}' full_table.tsv > Fragmented.C.c.busco.tab

cp Fragmented.C.c.busco.tab ~/temp

grep -w -f Frag.BUSCO.C.canephora.tab Fragmented.C.c.busco.tab | awk '{print $3}' > FRAG.C.c.exclusiveIDS.txt

#retrieve their Descriptions

grep -w -f FRAG.C.c.exclusiveIDS.txt C.canephora.GenDescriptions.tab> FRAG.C.c.exclusiveIDS.Descriptions.tab

#retrieve their GO terms

grep -w -f  FRAG.C.c.exclusiveIDS.txt GO.onePerRow.subcanephora.tab | awk '{if( $2 != "") print $0}' > FRAG.C.c.exclusiveIDS.GO.tab

################
#C.eugenioides##
################

busco -f -i subEugenioidesProteins.aa -l eudicots_odb10 -o C.e.buscoProtResult -m prot -c 24 --long --config  ~/bin/busco-4.1.4/config/config.ini


grep -w -f exclusive.BUSCO.C.eugenioides.tab complete.C.e.busco.tab | awk '{print $3}' > C.e.exclusiveIDS.txt

#retrieve their Descriptions

grep -w -f C.e.exclusiveIDS.txt C.eugeniodes.GenDescriptions.tab> C.e.exclusiveIDS.Descriptions.tab

grep -w -f C.e.exclusiveIDS.txt GO.onePerRow.subEugenoides.tab | awk '{if( $2 != "") print $0}' > C.e.exclusiveIDS.GO.tab

awk '{if($2 == "Fragmented") {print $0}}' full_table.tsv > Fragmented.C.e.busco.tab

cp Fragmented.C.e.busco.tab ~/temp



grep -w -f Fragmented.BUSCO.C.eugenioides.tab Fragmented.C.e.busco.tab | awk '{print $3}' > FRAG.C.e.exclusiveIDS.txt

#retrieve their Descriptions

grep -w -fFRAG.C.e.exclusiveIDS.txt C.eugeniodes.GenDescriptions.tab > FRAG.C.e.exclusiveIDS.Descriptions.tab

#retrieve their GO terms

grep -w -f  FRAG.C.e.exclusiveIDS.txt GO.onePerRow.subEugenoides.tab | awk '{if( $2 != "") print $0}' > FRAG.C.e.exclusiveIDS.GO.tab

#######################################################################

awk -F'\t' '{print $1"\t"$2}' direct_go_count_bp_subcanephoraproteins.txt  > Filter.direct_go_count_bp_subcanephoraproteins.txt

awk -F'\t' '{print $1"\t"$2}' direct_go_count_cc_subcanephoraproteins.txt  > Filter.direct_go_count_cc_subcanephoraproteins.txt

awk -F'\t' '{print $1"\t"$2}' direct_go_count_mf_subcanephoraproteins.txt  > Filter.direct_go_count_mf_subcanephoraproteins.txt

#Eugenioides
awk -F'\t' '{print $1"\t"$2}' direct_go_count_bp_subeugenioidesproteins.txt  > Filter.direct_go_count_bp_subeugenioidesproteins.txt

awk -F'\t' '{print $1"\t"$2}' direct_go_count_cc_subeugenioidesproteins.txt  > Filter.direct_go_count_cc_subeugenioidesproteins.txt

awk -F'\t' '{print $1"\t"$2}' direct_go_count_mf_subeugenioidesproteins.txt  > Filter.direct_go_count_mf_subeugenioidesproteins.txt

###################
R

subCanephoraBP <- read.delim("Filter.direct_go_count_bp_subcanephoraproteins.txt", header=T)
subCanephoraCC <- read.delim("Filter.direct_go_count_cc_subcanephoraproteins.txt", header=T)
subCanephoraMF <- read.delim("Filter.direct_go_count_mf_subcanephoraproteins.txt", header=T)

subEugenioidesBP <- read.delim("Filter.direct_go_count_bp_subeugenioidesproteins.txt", header=T)
subEugenioidesCC <- read.delim("Filter.direct_go_count_cc_subeugenioidesproteins.txt", header=T)
subEugenioidesMF <- read.delim("Filter.direct_go_count_mf_subeugenioidesproteins.txt", header=T)

#BP
setdiff(subEugenioidesBP$GO,subCanephoraBP$GO)

BParea1 <- length(subEugenioidesBP$GO)
BParea2 <- length(subCanephoraBP$GO)

BPcross <- length(intersect(subEugenioidesBP$GO,subCanephoraBP$GO))

length(setdiff(subEugenioidesBP$GO,subCanephoraBP$GO))

length(setdiff(subCanephoraBP$GO,subEugenioidesBP$GO))

#########
#CC

setdiff(subEugenioidesCC$GO,subCanephoraCC$GO)

union(subEugenioidesCC$GO,subCanephoraCC$GO)

length(setdiff(subEugenioidesCC$GO,subCanephoraCC$GO))

length(setdiff(subCanephoraCC$GO,subEugenioidesCC$GO))

CCarea1 <- length(subEugenioidesCC$GO)
CCarea2 <- length(subCanephoraCC$GO)

CCcross <- length(intersect(subEugenioidesCC$GO,subCanephoraCC$GO))

#########
#MF

setdiff(subEugenioidesMF$GO,subCanephoraMF$GO)

union(subEugenioidesMF$GO,subCanephoraMF$GO)

length(setdiff(subEugenioidesMF$GO,subCanephoraMF$GO))

length(setdiff(subCanephoraMF$GO,subEugenioidesMF$GO))

MFarea1 <- length(subEugenioidesMF$GO)
MFarea2 <- length(subCanephoraMF$GO)

MFcross <- length(intersect(subEugenioidesMF$GO,subCanephoraMF$GO))

############################################
###Drawn venn
###########################################

library(VennDiagram)
#Biological Process

draw.pairwise.venn(area1 = BParea1,
area2 = BParea2,
cross.area = BPcross,
fill = c("#BF7F58", "#698595"), scaled = TRUE,
cex = rep (0,3),
ext.text=F,
)

#Celular Component

draw.pairwise.venn(area1 = CCarea1,
area2 = CCarea2,
cross.area = CCcross,
fill = c("#BF7F58", "#698595"), scaled = TRUE,
cex = rep (0,3),
ext.text=F,
)
#Molecular Function

draw.pairwise.venn(area1 = MFarea1,
area2 = MFarea2,
cross.area = MFcross,
fill = c("#BF7F58", "#698595"), scaled = TRUE,
cex = rep (1,3),
ext.text=F,
)
#####################################
#GO quantity per sub genome
w=40/2.54
h=w/1.618
#merge the individual dataframes into a single one
#BP

bothBP <- merge(subEugenioidesBP,subCanephoraBP, by="GO", all.x=T, all.y=T, sort=F)

colnames(bothBP) <- c("GO", "C.eu.numberSeqs","C.ca.numberSeqs")
#Filter the TOP 10 and add the others as a 11th line

#othersBP <- colSums(bothBP[c(seq(-1,-10,-1)),c(2,3)],na.rm=T)

processedBothBP <- as.matrix(head(bothBP,10)[,c(2,3)])

#processedBothBP <- rbind(processedBothBP,othersBP)

rownames(processedBothBP) <- as.character(head(bothBP,10)$GO)

rownames(processedBothBP)[3] <- "regulation of transcription"
rownames(processedBothBP)[6] <- "metabolism of carbohydrates"


processedBothBP <- apply(processedBothBP, 2, rev)

#library(plotrix)

#Direct GO Count BP
pdf("Direct.BP.GO.Count.pdf",h=h,w=w)
par(mar=c(4,30,2,2))
bp <- barplot(t(processedBothBP),horiz=T,beside=T,xaxt='n', yaxt='n', xlim=c(0,4000),col=c("#BF7F58", "#698595"),space=c(0,0.5))

axis(side=1, at=seq(0,19000,1000),labels=F)
mtext(side=1,at=seq(0,19000,1000),text=seq(0,19,1), cex = 2.7, line = 1.6)

mtext(side=1,at=-370,text="1000 x ", cex = 2.7, line = 1.6)

axis(side=2, at=bp[1,]+0.5,labels=F)
mtext(side=2,at=bp[1,]+0.5,text=rownames(processedBothBP), cex = 2.7, line = .6, las=2)


dev.off()

system("open Direct.BP.GO.Count.pdf")

#########################################

#GO quantity per sub genome
w=40/2.54
h=w/1.618
#merge the individual dataframes into a single one
#CC

bothCC <- merge(subEugenioidesCC,subCanephoraCC, by="GO", all.x=T, all.y=T, sort=F)

colnames(bothCC) <- c("GO", "C.eu.numberSeqs","C.ca.numberSeqs")
#Filter the TOP 10 and add the others as a 11th line

#othersCC <- colSums(bothCC[c(seq(-1,-10,-1)),c(2,3)],na.rm=T)

processedBothCC <- as.matrix(head(bothCC,10)[,c(2,3)])

#processedBothCC <- rbind(processedBothCC,othersCC)

rownames(processedBothCC) <- as.character(head(bothCC,10)$GO)

rownames(processedBothCC)[1] <- "membrane component"
#rownames(processedBothCC)[6] <- "metabolism of carbohydrates"


processedBothCC <- apply(processedBothCC, 2, rev)

#library(plotrix)

#Direct GO Count CC
pdf("Direct.CC.GO.Count.pdf",h=h,w=w)
par(mar=c(4,30,2,2))
CC <- barplot(t(processedBothCC),horiz=T,beside=T,xaxt='n', yaxt='n', xlim=c(0,4000),col=c("#BF7F58", "#698595"),space=c(0,0.5))

axis(side=1, at=seq(0,19000,1000),labels=F)
mtext(side=1,at=seq(0,19000,1000),text=seq(0,19,1), cex = 2.7, line = 1.6)

mtext(side=1,at=-370,text="1000 x ", cex = 2.7, line = 1.6)

axis(side=2, at=CC[1,]+0.5,labels=F)
mtext(side=2,at=CC[1,]+0.5,text=rownames(processedBothCC), cex = 2.7, line = .6, las=2)


dev.off()

system("open Direct.CC.GO.Count.pdf")

#GO quantity per sub genome
w=40/2.54
h=w/1.618
#merge the individual dataframes into a single one
#MF

bothMF <- merge(subEugenioidesMF,subCanephoraMF, by="GO", all.x=T, all.y=T, sort=F)

colnames(bothMF) <- c("GO", "C.eu.numberSeqs","C.ca.numberSeqs")
#Filter the TOP 10 and add the others as a 11th line

#othersMF <- colSums(bothMF[c(seq(-1,-10,-1)),c(2,3)],na.rm=T)

processedBothMF <- as.matrix(head(bothMF,10)[,c(2,3)])

#processedBothMF <- rbind(processedBothMF,othersMF)

rownames(processedBothMF) <- as.character(head(bothMF,10)$GO)

rownames(processedBothMF)[5] <- "hybrid ribonuclease activity"
#rownames(processedBothMF)[6] <- "metabolism of carbohydrates"


processedBothMF <- apply(processedBothMF, 2, rev)

#library(plotrix)

#Direct GO Count MF
pdf("Direct.MF.GO.Count.pdf",h=h,w=w)
par(mar=c(4,30,2,2))
MF <- barplot(t(processedBothMF),horiz=T,beside=T,xaxt='n', yaxt='n', xlim=c(0,8000),col=c("#BF7F58", "#698595"),space=c(0,0.5))

axis(side=1, at=seq(0,19000,1000),labels=F)
mtext(side=1,at=seq(0,19000,1000),text=seq(0,19,1), cex = 2.7, line = 1.6)

mtext(side=1,at=-770,text="1000 x ", cex = 2.7, line = 1.6)

axis(side=2, at=MF[1,]+0.5,labels=F)
mtext(side=2,at=MF[1,]+0.5,text=rownames(processedBothMF), cex = 2.7, line = .6, las=2)


dev.off()

system("open Direct.MF.GO.Count.pdf")


######################################
#GO BP top terms
#SubCanephora
grep 'DNA integration' subCanephoraAnnot.tab | awk -F'\t' '{print $2}'  > C.c.DNAintegration.tab
#GO CC top term
grep 'integral component of membrane' subCanephoraAnnot.tab | awk -F'\t' '{print $2}'  >  C.c.integralComponentMembrane.tab
#GO MF to term
grep 'nucleic acid binding' subCanephoraAnnot.tab | awk -F'\t' '{print $2}'  >  C.c.nucleicAcidBinding.tab

cp C.c.nucleicAcidBinding.tab /Users/johnjoyce/Dropbox/myPapers/C.arabica.Gen.Prediction/misc/

cp C.c.integralComponentMembrane.tab /Users/johnjoyce/Dropbox/myPapers/C.arabica.Gen.Prediction/misc/

cp C.c.DNAintegration.tab /Users/johnjoyce/Dropbox/myPapers/C.arabica.Gen.Prediction/misc/

##################################

#GO BP top terms
#SubEugenioides
grep 'DNA integration' subEugenoidesAnnot.tab | awk -F'\t' '{print $2}'  > C.e.DNAintegration.tab
#GO CC top term
grep 'integral component of membrane' subEugenoidesAnnot.tab | awk -F'\t' '{print $2}'  >  C.e.integralComponentMembrane.tab
#GO MF to term
grep 'nucleic acid binding' subEugenoidesAnnot.tab | awk -F'\t' '{print $2}'  >  C.e.nucleicAcidBinding.tab

cp C.e.nucleicAcidBinding.tab /Users/johnjoyce/Dropbox/myPapers/C.arabica.Gen.Prediction/misc/

cp C.e.integralComponentMembrane.tab /Users/johnjoyce/Dropbox/myPapers/C.arabica.Gen.Prediction/misc/

cp C.e.DNAintegration.tab /Users/johnjoyce/Dropbox/myPapers/C.arabica.Gen.Prediction/misc/

#########
#Back to R

C.c.nucleicAcidBinding <- read.delim("C.c.nucleicAcidBinding.tab",header=F)

tableC.c.nucleicAcidBinding <- sort(table(C.c.nucleicAcidBinding))

C.c.integralComponentMembrane <- read.delim("C.c.integralComponentMembrane.tab",header=F)

tableC.c.integralComponentMembrane <- sort(table(C.c.integralComponentMembrane))

C.c.DNAintegration <- read.delim("C.c.DNAintegration.tab",header=F)

tableC.c.DNAintegration <- sort(table(C.c.DNAintegration))

########################

C.e.nucleicAcidBinding <- read.delim("C.e.nucleicAcidBinding.tab",header=F)

tableC.e.nucleicAcidBinding <- sort(table(C.e.nucleicAcidBinding))

C.e.integralComponentMembrane <- read.delim("C.e.integralComponentMembrane.tab",header=F)

tableC.e.integralComponentMembrane <- sort(table(C.e.integralComponentMembrane))

C.e.DNAintegration <- read.delim("C.e.DNAintegration.tab",header=F)

tableC.e.DNAintegration <- sort(table(C.e.DNAintegration))

#######################
#col=c("#BF7F58", "#698595")

ColorPieC.c <- colorRampPalette(c("#698595", "white"))(7)
ColorPieC.e <- colorRampPalette(c("#BF7F58", "white"))(7)

#nucleicAcidBinding

#C.c
otherstableC.c.nucleicAcidBinding <- sum(tableC.c.nucleicAcidBinding[c(seq(1,length(tableC.c.nucleicAcidBinding)-5,1))])
names(otherstableC.c.nucleicAcidBinding) <- "others"
processedC.c.nucleicAcidBinding<- rev(tail(tableC.c.nucleicAcidBinding,5))

processedC.c.nucleicAcidBinding <- c(processedC.c.nucleicAcidBinding,otherstableC.c.nucleicAcidBinding)

pie(processedC.c.nucleicAcidBinding ,labels='', col =ColorPieC.c)

#C.e"

otherstableC.e.nucleicAcidBinding <- sum(tableC.e.nucleicAcidBinding[c(seq(1,length(tableC.e.nucleicAcidBinding)-5,1))])
names(otherstableC.e.nucleicAcidBinding) <- "others"
processedC.e.nucleicAcidBinding<- rev(tail(tableC.e.nucleicAcidBinding,5))

processedC.e.nucleicAcidBinding <- c(processedC.e.nucleicAcidBinding,otherstableC.e.nucleicAcidBinding)

pie(processedC.e.nucleicAcidBinding ,labels='', col =ColorPieC.e)

##########################
#DNA integration

#C.c"
otherstableC.c.DNAintegration <- sum(tableC.c.DNAintegration[c(seq(1,length(tableC.c.DNAintegration)-5,1))])
names(otherstableC.c.DNAintegration) <- "others"
processedC.c.DNAintegration<- rev(tail(tableC.c.DNAintegration,5))

processedC.c.DNAintegration <- c(processedC.c.DNAintegration,otherstableC.c.DNAintegration)

pie(processedC.c.DNAintegration ,labels='', col =ColorPieC.c)

#C.e"

otherstableC.e.DNAintegration <- sum(tableC.e.DNAintegration[c(seq(1,length(tableC.e.DNAintegration)-5,1))])
names(otherstableC.e.DNAintegration) <- "others"
processedC.e.DNAintegration<- rev(tail(tableC.e.DNAintegration,5))

processedC.e.DNAintegration <- c(processedC.e.DNAintegration,otherstableC.e.DNAintegration)

pie(processedC.e.DNAintegration ,labels='', col =ColorPieC.e)

##########################
#integralComponentMembrane

#C.c"
otherstableC.c.integralComponentMembrane <- sum(tableC.c.integralComponentMembrane[c(seq(1,length(tableC.c.integralComponentMembrane)-5,1))])
names(otherstableC.c.integralComponentMembrane) <- "others"
processedC.c.integralComponentMembrane<- rev(tail(tableC.c.integralComponentMembrane,5))

processedC.c.integralComponentMembrane <- c(processedC.c.integralComponentMembrane,otherstableC.c.integralComponentMembrane)

pie(processedC.c.integralComponentMembrane ,labels='', col =ColorPieC.c)

#C.e"

otherstableC.e.integralComponentMembrane <- sum(tableC.e.integralComponentMembrane[c(seq(1,length(tableC.e.integralComponentMembrane)-5,1))])
names(otherstableC.e.integralComponentMembrane) <- "others"
processedC.e.integralComponentMembrane<- rev(tail(tableC.e.integralComponentMembrane,5))

processedC.e.integralComponentMembrane <- c(processedC.e.integralComponentMembrane,otherstableC.e.integralComponentMembrane)

pie(processedC.e.integralComponentMembrane ,labels='', col =ColorPieC.e)
