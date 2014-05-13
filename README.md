paired_end_demux
================

This python script is used to take input FASTQ files that contain indexed samples and split the samples by index.

A single-base mismatch is allowed for indexes by default, but can be disabled. If a set of reads matches the index, the ID for the read is appropriately adjusted to reflect what the index sequence was (whether an exact match or not). Reads are also filtered for N-base content, i.e. reads with high N composition are dropped. Reads without an index are printed to a separate file. Summary statistics are printed to standard out at the end of a run, including total reads, reads with index, counts by index, and amount of total sequence.

##Requirements
- Input in FASTQ format

- Input with Phred quality scale

- Paired-end reads

- Indexed

 -Appended to read 1 *or*
	
 -Third FASTQ with index read 

##Usage
	# index appended to read 1
    python paired_end_demux r1.fq r2.fq ATCACGA,CGATGTA,TTAGGCA 100 output > log

	# index as separate FASTQ
	python paired_end_demux r1.fq r2.fq ATCACGA,CGATGTA,TTAGGCA 100 output --indexRead ndx.fq > log
	
	# index appended to read 1, dropping reads with > 50% N bases
	python paired_end_demux r1.fq r2.fq ATCACGA,CGATGTA,TTAGGCA 100 output --nPercentCutoff 0.50 > log
	
	# index in separate FASTQ, disallow mismatches to index, drop reads with > 75% N content
	python paired_end_demux r1.fq r2.fq ATCACGA,CGATGTA,TTAGGCA 100 output --indexRead ndx.fq --nPercentCutoff 0.75 --noIndexMismatch > log
