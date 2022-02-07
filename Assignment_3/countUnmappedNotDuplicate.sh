#!/bin/bash
set -e
bam1=$1
echo Read length is:
echo $bam1 | cut -d "_" -f 7 
echo And the count is: 
sambamba view -F "[XS] == null and not unmapped and not duplicate" /projects/bmeg/A3/$bam1 -c
