# Exercise 2: BioSequences.jl
This time, you'll be doing some basic processing of biological sequences using BioSequences.jl.

### Exercise 2.1: Load in the data
The file `sequences.json` is a JSON file containing the human chromosome Y, and a vector of `(name, sequence)` pairs, where `sequence` is a HA segment from Influenza viruses, and `name` is the corresponding identifier.

Search the internet for the JSON3 package documentation, and figure out how to load in the data in `sequences.json` into Julia.

* Q 2.1.1: What Julian data type are the DNA sequences stored as in the JSON file?

Find a way to convert them to `LongSequence{DNAAlphabet{4}}`.

* Q 2.1.2: How many HA segments are there in the file?

* Q 2.1.3: What is the sum of the lengths of the HA segments?

### Exercise 2.2: Human chromosome Y
Look at the human chromosome Y, named `:chromy` in the serialization file.

* Q 2.2.1: How many total times do each of the following symbols appear in the provided human chromosome Y?: W, Y, C, N.
Try to solve it without converting to Strings!

Consider all substrings of length 10 that does not contain `N`.

* Q 2.2.2: How many of these substrings are palindromic, i.e. they are the same as their reverse complement?

The _alu_ element is a parasitic transposable genetic element in supraprimates, including humans. 
It is present about 1 million times in the human genome, comprising 10% of human DNA.
_alu_ sequences include a promoter box with the sequence "TGGCTCACGCC".

* Q 2.2.3: How many times does this exact sequence occur in the provided reference human chromosome Y?
Skip bases called as "N" and remember to count reverse-complemented sequences.

### Exercise 2.3: Influenza HA segments
* Q 2.3.1: How many of these sequences could be encoded by `DNAAlphabet{2}`?

* Q 2.3.2: Exactly one sequence does not start with `AGC`. Which?

### Exercise 2.4: The genetic code
RNA can be translated to proteins in a process called _translation_. BioSequences.jl has methods for this.

* Q 2.4.1: What does this RNA sequence translate to using the _standard genetic code_? `AUGCGCGGAGAUAUAUUA`

Consider this following DNA sequence:

ATGATCTCTATGGTATGGAATATGTGTAACGTTCGTGCAGATGCACCTATGGCGTGGCAGAAATTGTTTC
AAGACCCTGCAACGTCCAACATGGAAGGAATTGTAGACCTTCACGCAGATATTTGCTTTTTTTTGATTGT
AATTTTGGTTTTGGTTTTGTGGCTCGGAGCTCGAATTCTTGTAAGTTTTCATCACAACTTGCAACCTGTA
CCAGAGCGTTTTAACCACCACACTTCTTTGGAATTGGTTTGGGCAATTTTGCCTTCTGTTATTGTAACAT
TGATTGCATTGCCTTCCTAGTCTCTTGTGTATACTTATGATGATTAGGTAAGTAAACCCGCATTGACCGT
GAAGGTTACTGGGCGTCAGTGGTATTGGAGCTACGCTATGAATGAGCATGTACAAATGAATTTGAGTCAA
CAGGCGAAAGATTTGTTGCTTCAGTCTTCA

This is a gene (coding sequence, excluding STOP codons) from a mitochrondrial genome.

* Q 2.4.2: Is the sequence from a vertebrate mitochondrion or not?
_Hint: Look at the different genetic codes. How can you tell if a certain code is used for this gene?_

The next question pertains to the influenza HA sequence named "AY650270" in your dataset. The coding part of the sequence is `22:1701`.
When doing analysis of selection pressure (dN/dS ratio) in DNA sequences, one needs to count the number of possible mutations that cause amino acid substitutions versus the number of possible silent mutations.

* Q 2.4.3: How many point mutations (that is, a mutation that changes a single base pair) is possible in this sequence which will NOT alter the protein it is translated to?

* Q 2.4.4: How many different RNA sequences can code for this small protein using the standard genetic code: `RTHLMPYDQDISYIMDHKGP`?
