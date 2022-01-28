
#!/bin/bash
set -e # this makes the whole script exit on any error.
#fill these variables with the arguments given to this script
userDir = alzhang_bmeg22
sample=$1
fq1=$2
fq2=$3
logDir=MyLogDirectory # this is where all the files to keep track of progress will go.
mkdir -p MyLogDirectory # make the directory where log files will go, if it doesn't exist already
echo running pipeline for $sample
if [ ! -e $logDir/$sample.fastqc.done ] #run this code only if $logDir/$sample.fastqc.done is missing
then
    	echo Performing fastqc of sample $sample with the following fastqs:
        ls /projects/bmeg/A2/$fq1 /projects/bmeg/A2/$fq2


        #enter commands to run fastqc here
		fastqc /projects/bmeg/A2/$fq1 --outdir=/home/$userDir/assignment2
		fastqc /projects/bmeg/A2/$fq2 --outdir=/home/$userDir/assignment2

	

        touch $logDir/$sample.fastqc.done #create the file that we were looking for at the beginning of this if statement so that this same code is not run next time
else # $logDir/$sample.fastqc.done was not missing
        echo Already performed fastqc of $sample
fi

if [ ! -e $logDir/$sample.align.done ]
then
		echo Running Alignment for $sample now, please wait around 30 minutes.

		bowtie2 -x /projects/bmeg/indexes/hg38/hg38_bowtie2_index -1 /projects/bmeg/A2/$fq1 /projects/bmeg/A2/$fq2 -S /home/$userDir/assignment2/$sample.align.sam
        
        touch $logDir/$sample.align.done #create the file that we were looking for at the beginning of this if statement so that this same code is not run next time
else # $logDir/$sample.fastqc.done was not missing
        echo Already performed alignment of $sample

fi
	echo Converting the sam file to bam file now
	samtools view -S -b -h $sample.align.sam > $sample.align.bam

	echo Sorting the bam file. And removing the sam file
	sambamba sort $sample.align.bam
	rm $sample.align.sam

	echo Now, filtering the sorted bam file 

	sambamba view -F "[XS] == null and not unmapped and not duplicate" -h --format=bam  $sample.align.sorted.bam  > $sample.align.sorted.filtered.bam

	echo The sorted and filtered file is now saved as $sample.align.sorted.filtered.bam

	echo The filtered and sorted file has the following reads:

	sambamba view -c $sample.align.sorted.filtered.bam 





