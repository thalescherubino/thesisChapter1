#!/bin/sh

#  SubArabica.predition.sh
#
#
#  Created by THCRibeiro on 23/01/2021.
#
mkdir subEugenioides
cd subEugenioides

/home/tcherubino/miniconda3/bin/new_species.pl --species=subEugenioides #--AUGUSTUS_CONFIG_PATH=/home/lfmp/bin/Augustus/config/


sed 's/\s.*$//' ../genome/eugenioidesChromossomes.fa > ../genome/c.e.fa


export PATH=$PATH:/home/tcherubino/bin/gth-1.7.1-Linux_x86_64-64bit/bin

/home/tcherubino/bin/BRAKER/scripts/startAlign.pl  --CPU=24 --genome=../genome/c.e.fa --prot=../all.transcriptomes/0.8.finalProteins.fasta --prg=gth


#Convert align_gth/gth.concat.aln to GTF format:
/home/tcherubino/miniconda3/bin/gth2gtf.pl gth.concat.aln bonafide.gtf


#Compute a flanking region length for training gene structures in GenBank flatfile format

/home/tcherubino/miniconda3/bin/computeFlankingRegion.pl bonafide.gtf

#Total length gene length (including introns): 13121831. Number of genes: 6550. Average Length: 2003.33297709924

#The flanking_DNA value is: 1001 (the Minimum of 10 000 and 1001)

/home/tcherubino/miniconda3/bin/gff2gbSmallDNA.pl bonafide.gtf ../../genome/c.e.fa 1001 bonafide.gb

/home/tcherubino/miniconda3/bin/randomSplit.pl bonafide.gb 1200

mv bonafide.gb.test test.gb
mv bonafide.gb.train train.gb

/home/tcherubino/miniconda3/bin/etraining --species=subEugenioides train.gb &> etrain.out


#nano /home/tcherubino/miniconda3/config/species/subEugenioides/subEugenioides_parameters.cfg

#Frequency of stop codons:
#tag:  858 (0.23)
#taa: 1316 (0.353)
#tga: 1551 (0.416)


/home/tcherubino/miniconda3/bin/augustus --species=subEugenioides test.gb > test.out


grep -c LOCUS bonafide.gb
#6203
grep -c "Variable stopCodonExcludedFromCDS set right" etrain.out
#955

################################################
#GENERATING TRAINING GENE STRUCTURES FROM ESTs##
################################################
#Alternate protocol 2

#I now believe that I should not use only the CDS because this will remove the UTRs and this will fuck with the UTR prediction.... So some tinkering will be done to copy the complete EST using the IDs from complete CDSs

#filter only complete sequences

grep complete 0.80.colapsed.all.transcriptomes.fasta.transdecoder.cds | awk '{print $1}' > complente.cds.txt


awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < 0.80.colapsed.all.transcriptomes.fasta > temp.fasta


cat complente.cds.txt |  awk -F '.p' '{print $1}' > new.complente.cds.txt

grep -A 1 -F -f new.complente.cds.txt temp.fasta > 0.8.finalEST.fasta

scp 0.8.finalEST.fasta lfmp@177.105.2.190:/home/lfmp
scp c.e.fa lfmp@177.105.2.190:/home/lfmp


mkdir subEugenioides
cd subEugenioides
cp ~/0.8.finalEST.fasta ./
mv ~/c.e.fa ./

/home/lfmp/bin/PASApipeline.v2.4.1/bin/seqclean 0.8.finalEST.fasta


cp /home/lfmp/bin/PASApipeline.v2.4.1/pasa_conf/pasa.alignAssembly.Template.txt  ./alignAssembly.config

#Edit the file alignAssembly.
#
# database settings
#DATABASE=ESTsubEugenioides_pasa

sudo /home/lfmp/bin/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g c.e.fa -t 0.8.finalEST.fasta.clean --ALIGNERS blat --CPU 8

#Call open reading frames (ORFs) in PASA assemblies as follows:

~/bin/PASApipeline.v2.4.1/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta  ESTsubEugenioides_pasa.assemblies.fasta --pasa_transcripts_gff3  ESTsubEugenioides_pasa.pasa_assemblies.gff3



#Create a list of complete ORFs:
grep complete  ESTsubEugenioides_pasa.assemblies.fasta.transdecoder.cds | sed 's/>//' | awk '{print $1}' > complete.orfs

grep -F -f complete.orfs  ESTsubEugenioides_pasa.assemblies.fasta.transdecoder.genome.gff3 | grep -P "(\tCDS\t|\texon\t)" | perl -pe 's/cds\.//; s/\.exon\d+//; ' >  trainingSetComplete.gff3

mv trainingSetComplete.gff3 trainingSetComplete.temp.gff3

cat trainingSetComplete.temp.gff3 | perl -pe 's/\t\S*(asmbl_ \d+).*/\t$1/' | sort -n -k 4 | sort -s -k 9 | sort -s -k 1,1 > trainingSetComplete.gff3

cp trainingSetComplete.gff3 bonafide.gtf


#Follow steps 3 and 4 of Alternate Protocol 1.

/home/lfmp/bin/Augustus/scripts/computeFlankingRegion.pl bonafide.gtf


#Total length gene length (including introns): 53964810. Number of genes: 20088. Average Length: 2686.42025089606

#The flanking_DNA value is: 1343 (the Minimum of 10 000 and 1343)

/home/lfmp/bin/Augustus/scripts/gff2gbSmallDNA.pl bonafide.gtf c.e.fa 1343 bonafide.gb

#Now lets transfeer the bonafide.gb to Fiocruz!

scp bonafide.gb tcherubino@bagre.cpqrr.fiocruz.br:/home/tcherubino/doc/subCanephora/EST

/home/tcherubino/miniconda3/bin/randomSplit.pl bonafide.gb 1900


mv bonafide.gb.test test.gb
mv bonafide.gb.train train.gb

/home/tcherubino/miniconda3/bin/etraining --species=subEugenioides train.gb &> etrain.out


/home/tcherubino/miniconda3/bin/augustus --species=subEugenioides test.gb > test.out


#SUPPORT PROTOCOL 5 - UTR TRAINING
#UTRparameters are not required for gene prediction withAUGUSTUS.When integrating expression data, however, the accuracy ofAUGUSTUS may benefit fromUTRparameters (reduction of false positive CDS feature prediction in long UTRs).

#I will try in the folder used to predict based on EST

mv  test.out no.UTR.out


mv bonafide.gb.test test.gb
mv bonafide.gb.train train.gb

/home/tcherubino/miniconda3/bin/optimize_augustus.pl  --cpus=24 --species=subEugenioides --kfold=24 train.gb --trainOnlyUtr=1 --AUGUSTUS_CONFIG_PATH=/home/tcherubino/miniconda3/config/ --metapars=/home/tcherubino/miniconda3/config/species/subEugenioides/subEugenioides_metapars.utr.cfg | tee optimizeReport.UTR.out

/home/tcherubino/miniconda3/bin/augustus --species=subEugenioides  test.gb --UTR=on --print_utr=on > test.utr.out


#optimize other parameters than UTR

/home/tcherubino/miniconda3/bin/optimize_augustus.pl  --cpus=24 --species=subEugenioides --kfold=24 train.gb --trainOnlyUtr=0 --AUGUSTUS_CONFIG_PATH=/home/tcherubino/miniconda3/config/ | tee optimizeReport.Final.out


cd /home/tcherubino/doc/RNAseqBasedPrediction/ceu.subgenome

samtools view -@ 24 -b ../C.a.merged.bam NC_039899.1 NC_039901.1 NC_039903.1 NC_039905.1 NC_039906.1 NC_039909.1 NC_039911.1 NC_039912.1 NC_039915.1 NC_039916.1 NC_039919.1 > C.eugenioidesSubGenome.bam

samtools sort -@ 24 -n C.eugenioidesSubGenome.bam > name.sorted.bam

/home/tcherubino/miniconda3/bin/filterBam --uniq --paired --pairwiseAlignment --in name.sorted.bam --out Aligned.out.ssf.bam

rm name.sorted.bam

samtools sort -@ 24 Aligned.out.ssf.bam -o Aligned.out.ss.bam

rm Aligned.out.ssf.bam

/usr/local/bioinformatic/miniconda2/bin/bam2hints --intronsonly --in=Aligned.out.ss.bam --out=introns.gff

#Most RNA-seq data is by nature unstranded, i.e., it is unclear from which strand an aligned read stems. We try to guess the correct strand by using genomic splice site information (and discard all reads that do not have appropriate splice site information)

#Please check the headers of the genome FASTA file. If the headers are long and contain whitespaces, some RNA-Seq alignment tools will truncate sequence names in the BAM file. This leads to an error with BRAKER. Solution: shorten/simplify FASTA headers in the genome file before running the RNA-Seq alignment and BRAKER.

/home/tcherubino/miniconda3/bin/filterIntronsFindStrand.pl ../../genome/c.e.fa introns.gff --score > introns.f.gff 2> error.gff

rm error.gff

~/bin/gmes_linux_64/gmes_petap.pl --verbose --sequence=../../genome/c.e.fa --ET=introns.f.gff --cores=24 #--et_score=0 --debug

#Next, we will filter the GeneMark-ET predictions for those that have support in all introns by RNA-seq data (single-exon genes will be selected randomly in an appropriate proportion):


/home/tcherubino/miniconda3/bin/filterGenemark.pl genemark.gtf introns.f.gff

ln -sf genemark.f.good.gtf bonafide.gtf

#AUGUSTUS requires not only information about the coding sequence, but also about the non-coding sequence.We provide information about non-coding sequence in the flanking region of genes in the GenBank flatfile. The flanking region chosen should not be too short. On the other hand, setting it to maximum length may increase the runtime of training a lot. As a rule of thumb, you may use half of the mRNA size. Compute this value, e.g., with computeFlankingRegion.pl:

/home/tcherubino/miniconda3/bin/computeFlankingRegion.pl bonafide.gtf


#Total length gene length (including introns): 20569872. Number of genes: 44601. Average Length: 461.197551624403

#The flanking_DNA value is: 230 (the Minimum of 10 000 and 230)
#Convert all GeneMark gene structures to GenBank flatfile format. Whenever two neighboring genes are separated by fewer than 2*flanking_DNA bases, the end of the flanking DNA of both is set to the midpoint between the two genes. Therefore, we do not use the filtered GeneMark gene structures here, because by using all of them we exclude potential coding regions from the flanking region of AUGUSTUS training genes, which will increase specificity of AUGUSTUS:

/home/tcherubino/miniconda3/bin/gff2gbSmallDNA.pl genemark.gtf ../../genome/c.e.fa 230 tmp.gb

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

perl /home/tcherubino/miniconda3/bin/gtf2aa.pl ../../genome/c.e.fa bonafide.f.gtf prot.aa

#Next, BLAST all training gene amino acid sequences against themselves and output only those protein sequences that are less than 80% redundant with any other sequence in the set:


/home/tcherubino/miniconda3/bin/aa2nonred.pl --cores=24 prot.aa prot.nr.aa

grep ">" prot.nr.aa | perl -pe 's/>//' > nonred.lst


cat bonafide.gb | perl -ne ' if ($_ =~ m/LOCUS\s+(\S+)\s/) { $txLocus = $1;
} elsif ($_ =~ m/\/gene=\"(\S+)\"/) { $txInGb3{$1} = $txLocus
}
if(eof()) { foreach (keys %txInGb3) { print "$_\t$txInGb3{$_}\n";
} }' > loci.lst

grep -f nonred.lst loci.lst | cut -f2 > nonred.loci.lst
#Filter bonafide.gb to keep only those loci that are contained in nonred.loci.lst:

/home/tcherubino/miniconda3/bin/filterGenesIn.pl nonred.loci.lst bonafide.gb > bonafide.f.gb


#bonafide.gb 876586
#bonafide.f.gb 775208
# 88%

mv bonafide.f.gb bonafide.gb

# etraining

/home/tcherubino/miniconda3/bin/etraining  --species=subEugenioides bonafide.gb &> bonafide.out

#Count how many times the output file bonafide.out contains the expression Variable stopCodonExcludedFromCDS set right that is actually part of an error message
grep -c "Variable stopCodonExcludedFromCDS set right" bonafide.out

#0!

grep -c LOCUS bonafide.gb
#17136
#otimo, NEHNUM!

#3000 genes for test

/home/tcherubino/miniconda3/bin/randomSplit.pl bonafide.gb 3000


mv bonafide.gb.test test.gb
mv bonafide.gb.train train.gb

/home/tcherubino/miniconda3/bin/etraining --species=subEugenioides train.gb &> etrain.out


/home/tcherubino/miniconda3/bin/augustus --species=subEugenioides  test.gb > test.out


/home/tcherubino/miniconda3/bin/augustus --species=subEugenioides  test.gb --UTR=on > test.utr.out

#tag: 4243 (0.3)
#taa: 4117 (0.291)
#tga: 5776 (0.409)


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

#THIS WAS DONE BEFORE! SO JUST GO TO THE RESPECTIVE PROTEIN FOLDER

cp prot.hints ~/doc/RNAseqBasedPrediction/ceu.subgenome/


#ALTERNATE PROTOCOL 8
#GENERATING HINTS FROM ESTs AND MEDIUM SIZE RNA-SEQ READS ESTs

#STs and medium size RNA-seq reads (>ß450 nt) are suitable for generating intron, exonpart (sometimes abbreviated as “ep”), and exon hints.

cd ../EST

cp /home/tcherubino/doc/all.transcriptomes/0.80.colapsed.all.transcriptomes.fasta.transdecoder.cds ./


#Align ESTs/cDNAs/medium sized RNA-seq reads to genome using BLAT:
/home/tcherubino/bin/blat -noHead -minIdentity=92 ../../genome/c.e.fa 0.80.colapsed.all.transcriptomes.fasta.transdecoder.cds cdna.psl


#Filter alignments to obtain those that are potentially most useful for gene prediction:
/home/tcherubino/bin/pslCDnaFilter -minId=0.9 -localNearBest=0.005 -ignoreNs -bestOverlap cdna.psl cdna.f.psl


#Sort filtered alignments according to start position and target sequence name

cat cdna.f.psl | sort -n -k 16,16 | sort -s -k 14,14 > cdna.fs.psl

#Convert psl-file to hints: bash

/home/tcherubino/miniconda3/bin/blat2hints.pl --in=cdna.fs.psl --out=cdna.hints --minintronlen=35 --trunkSS

cd ../../subEugenioides/

mkdir finalPrediction
cd finalPrediction


cp ../align_gth/prot.hints ./
cp ../EST/cdna.hints ./
cp ../../RNAseqBasedPrediction/ceu.subgenome/rnaseq.gff ./

cp ../../RNAseqBasedPrediction/ceu.subgenome/introns.f.gff ./

sed -i 's/pri=4/pri=3/g' rnaseq.gff
sed -i 's/pri=4/pri=5/g' prot.hints

cat prot.hints cdna.hints rnaseq.gff introns.f.gff > hints.gff


#Hints from different sources can still be distinguished by augustus via the src attribute in the last column, e.g., with src=P and src=E for hints from proteins and RNA-seq, respectively.

cp /home/tcherubino/miniconda3/config/extrinsic/extrinsic.M.RM.E.W.P.cfg ./

/home/tcherubino/miniconda3/bin/augustus --species=subEugenioides --UTR=on --extrinsicCfgFile=extrinsic.M.RM.E.W.P.cfg  --allow_hinted_splicesites=atac ../../genome/c.e.fa --codingseq=on --outfile=predicted.C.eugenioides.gff --progress=true --genemodel=complete --hintsfile=hints.gff

/home/tcherubino/miniconda3/bin/getAnnoFasta.pl --seqfile ../../genome/c.e.fa predicted.C.eugenioides.gff

sed 's/>/>SubC.e_/g' predicted.C.eugenioides.aa > subEugenioidesProteins.aa

busco -f -i subEugenioidesProteins.aa -l eudicots_odb10 -o buscoProtResult -m prot -c 24 --long

#remove unecessary data from GFF.

sed '/Delete/d' predicted.C.eugenioides.gff > Temp

sed -i"" '/Error/d' Temp

sed -i"" '/b2h/d' Temp

mv Temp processed.C.eugenioides.gff

sed -i"" 's/"g/"SubC.e_/g' processed.C.eugenioides.gff
