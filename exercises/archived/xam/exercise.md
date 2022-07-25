# Exercise 5: XAM.jl
For this exercise, you will be working with SAM and BAM files.

### Exercise 5.1: Look at the SAM file
SAM files are usually huge compared to BAM files while storing the same information, and therefore don't see much use.

They are usually useful because they are human-readable, so that's what we'll use it for today.

The file `subset.bam` was created by mapping the _M. abscessus_ reads from the FASTX exercise to the _M. abscessus_ reference genome.
Because BAM files are big, the provided BAM file only include some of the reads from the FASTQ file.

The file `head.sam` contains the first few lines from the SAM representation of `subset.bam`.

Open the SAM file and look at it. At the top are the header lines beginning with @.
The lines beginning with `@SQ` contain all the reference sequences that were mapped against - in practise, the identifiers and sequence lengths of the FASTA file used as a reference.

* Q 5.1.1: There is only one SQ-line in the file. What is its identifier and sequence length?

Look at the first line that represents an alignment.
Given the following BAM flags:

       1: Is paired
       2: Is proper pair
       4: Is unmapped
       8: Is next read in group unmapped
      16: Is reversed
      32: Is next read in group reversed
      64: Is first read in group
     128: Is last read in group
     256: Is secondary alignment
     512: Is quality controlled failure
    1024: Is PCR or optical duplicate
    2048: Is supplementary alignment

* Q 5.1.2: Which flags does the first alignment have?
* Q 5.1.3: What is the CIGAR string of the first alignment?

In the cigar string, `H` means "hard clip" - that is, part of the read is not part of the alignment.
`M` means match, meaning it matches the reference.

* Q 5.1.4: Look at the leftmost position of the first alignment. What is it?
* Q 5.1.5: Given the position of the first read, and the fact that the reference is the reference genome of _M. abscessus_, why do you think the alignment is clipped? Think of the biology of bacterial genomes.

### Exercise 5.2: Read the BAM file.
Open the `subset.bam` and store the records an a variable.

* Q 5.2.1: How many records are in the file?

Some alignments are marked "unmapped" or "next read in group unmapped" "supplementary". Filter away these reads.

* Q 5.2.2: How many reads were filtered away?

Get the largest position of any read in the file (that is, the largest "leftmost mapping position").

* Q 5.2.3: Using this information, the number of reads, and the mean read length, what is the depth in this file? How well does this match with the depth estimated in the FASTX exercise?

Notice that SAM/BAM files, include following information:
    * Read name (also known as template name)
    * Read sequence
    * Read quality

This is also the information stored in FASTQ files.
Hence, we can create FASTQ files from BAM files.

* Q 5.2.4: Create a FASTQ file from the filtered reads.