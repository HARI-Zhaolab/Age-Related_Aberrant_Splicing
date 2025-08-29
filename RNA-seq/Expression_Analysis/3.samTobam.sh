cat sample.txt | parallel -j 6 "samtools sort -o ./4.samTobam/{}.bam ./3.align/{}.sam > ./4.samTobam/{}_sort.log 2>&1"
