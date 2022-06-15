# Exercise 3: FASTX
In this exercise, you'll learn how to read and write FASTA and FASTQ files.
In the last part, you'll also see a concrete usecase for kmers, by creating a simple kmer spectrum from some FASTQ data.

You will notice the files are compressed with gzip (they have the .gz extension).
This is common with FASTA and FASTQ files, because their size can make storage expensive, and because they compress well due to them being plaintext files.
To read the files, you need to first import `CodecZlib` and use `GzipDecompressorStream`.
Here is an example

```julia
using CodecZlib
reader = FASTA.Reader(GzipDecompressorStream(open(path)))
# [ do stuff ... ]
close(reader) # always close after use
```

### Exercise 3.1: FASTA files
* Q 3.1.1: What is the identifier and description of the fifth-from-last record in `ns.fna`?

* Q 3.1.2: How many records in `ns.fna` have both ID and description?

* Q 3.1.3: How many records have a sequence length longer than the average in `ns.fna`?

* Q 3.1.4: Load in `ns.fna` into a `Vector{Tuple{String, LongDNA{4}}}` with the ID (identifier, not description) and sequence of each FASTA record.

* Q 3.1.5: Write a new file `ns.revcomp.fna` which contains the same record in the same order as `ns.fna`, but with all the sequence reverse-complemented, and where all records have had their description removed

### Exercise 3.2: FASTQ files
The file `M_abscessus.fq.gz` contains DNA reads originating from the bacterium _M. abscessus_, a mycobacterium capable of causing persistent infections in humans.

The file is large enough that you might want to consider performance of the code you write, but still small enough that you don't need to do any serious optimization.
It is possible to make code for these exercises very fast.

* Q 3.2.1: What is the mean and median length of reads in this file?

* Q 3.2.2: How many total basepairs are in the file?

* Q 3.2.3: How big a fraction of the reads contain only A, C, G or T in their sequences?

### Exercise 3.3: Make a kmer spectrum
For this exercise, let's consider 21-mers (of type `DNAKMer{21}`).

* Q 3.3.1: How many total 21-mers in the _M. abscessus_ file begin with `AAA`? Remember to count both the 21-mers and their reverse complement.

A "canonical" kmer is the lexographically smallest of a kmer and its reverse-complement. In all the following exercises, consider only canonical kmers.
_Hint: See the function `canonical`_

* Q 3.3.2: How many distinct 21-mers are there in the file? How many total kmers?

Make kmer spectrum like I showed on the slide and identify where the peak is.
You can use your favorite Julian plotting package to plot the spectrum if you want, but it's not necessary, you can just look at the raw data.

* Q 3.3.3: At what X-position is the rightmost peak? In other words, what the most common coverage for any given kmer in the genome?

* Q 3.3.4: You will notice there are two peaks. The rightmost peak is at around the mode (most common) depth of the genome.
The leftmost peak represents kmers that are only observed a few times.
What do you think could cause the leftmost peak?

In the following calculations, disregard all kmers in the leftmost peak.

Let's pretend you're about to assemble the genome using an old school assembler.
For this you will need to calculate the depth `D`.

First, calculate `P`, the mean number of times a kmer has been observed using the data in the right peak.
Be careful! This is not the same as the mode of the peak!
Rather, it's the total number of kmers divided by the number of distinct kmers.
Again, we disregard any kmers in the leftmost peak.

`D` can now be calculated with:

D = PT/Tk
where
- D = depth
- P = Mean number of times a kmer has been observed
- T = Total number of base pairs in the FASTQ file
- Tk = Total number of kmers (disregarding the leftmost peak)

* Q: 3.3.5: What is the depth `D`?

With the depth, you can estimate the genome size `S` of _M. abscessus_ using S = T/D.

* Q 3.3.6: What is the genome size?

Use a search engine to look up the genome size of the bactierum (it's Mycobacterium abscessus subsp. bolletii BDT)
Verify this estimate is within a few percent of the reported genome size of the bacterium.
