# Exercise 2: Kmers
This exercise pertains to kmers.

From a theoretical point of view, a kmer is any biosequence with a length `k`.

In bioinformatics, kmers are important, because if the value of `k` is fixed and small, kmers can be stored in memory as integers which is far more efficient than `LongSequence`, which stores the data using an array.

From a practical point of view, a kmer is a small, immutable `BioSequence`.
Besides being immutable, they should be very similar to the `LongSequence` you got familiar with in the last exercise.

### Exercise 2.1: Basic exercises using kmers
RNA codons can be though of as RNA 3-mers.
Indeed, `RNACodon` is an alias in the package Kmers for a length 3 kmer of the alphabet `RNAAlphabet{2}`.

* Q 2.1.1: Create a `Vector` containing all 64 different RNA codons. What is its type?

* Q 2.1.2: Filter away all codons C whose reverse complement is before C when put in alphabetical order (i.e non-canonical codons). How many are left?

* Q 2.1.3: Create a `Dict` that maps codons to their corresponding amino acid in the standard genetic code.

### Exercise 2.2: Efficient computation using kmers
If kmers are restricted to non-ambiguous bases, that is A, C, G, T (DNA only) or U (RNA only), it can be stored using only 2 bits per base, increasing efficiency.

For this reason, it is common to ignore all kmers containing ambiguous bases like `N` or `W` when extracting kmers from a sequence.

In the Kmers package, the function `EveryKmer{T}(s::BioSequence)` where `T` is a `Kmer` type will automatically skip ambiguous bases if `T` uses a 2-bit alphabet.

In the biosequence exercise, you counted the number of palindromic subsequences of length 10 does not contain `N` in the human chromosome Y.

* Q 2.2.2: Repeat the exercise, this time skipping all ambiguous bases, and using kmers.

You can time functions with the `@timed` macro.
For example, `@timed 1+1` will time how long it takes to run `1+1`.

Put your solution from the biosequences exercise in a function, and the equivalent solution you just did using kmers in another function, and time both using `@timed`.

* Q 2.2.2: How fast are each of the solutions?

* Q 2.2.3: Also repeat the biosequence exercise where you searched for the alu sequence. Time the two solutions and report the timings.

We will return to kmers in the FASTX exercise.
