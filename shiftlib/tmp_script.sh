#!/bin/bash

/home/tools/bedtools/bedtools-v2.16.2/bin/intersectBed -a $1 -b $2 -wao | \
  awk '$23 != 0 {print $4"\t"$5"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$20"\t"$21"\t"$23}' | \
  /home/tools/bedtools/bedtools-v2.16.2/bin/intersectBed -a stdin -b $3 -wao | \
  awk '($22 != 0) && ($4 == $19) {print $5"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$9"\t"$11"\t"$12"\t"$13"\t"$14"\t"$20"\t"$21"\t"$15"\t"$22"\t"$10}' | \
  sort | uniq
