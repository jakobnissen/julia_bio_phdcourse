# Exercise 1: PopGen.jl
Welcome to the first exercise!
This time, you'll be doing some basic population genetics with PopGen.jl

PopGen's functions work on data of the abstract type `PopObj`.
Currently, the only concrete subtype of `PopObj` is `PopData`.
Hence, you will be working a lot with `PopData` objects in this exercise.

### Exercise 1.1: Load in the data
The file `lobsters.vcf.gz` in the data directory is a gzipped VCF-file containing genetic markers of some wild American lobsters (_Homarus americanus_).

Use the `vcf` function to load in the file as a `PopData` object.

* Q 1.1.1: Inspect the type by typing it into the terminal. How many loci (SNP markers) and samples are there in the file?

First, get a feeling of how the object is stored in memory.

* Q 1.1.2: What fields do the `PopData` have? What are their types?

You will note when you loaded the data in there was a warning that:

> population info must be added

Look at the sample names of the lobsters, e.g. `BON_01`.
This name means "Lobster 1 sampled at location BON".
Let's consider each sampling location a population of lobsters.

* Q 1.1.3: Given the string `"BON_01"` composed of a prefix, and underscore, and a suffix, how do you extract the prefix? 

The method `populations!(::PopData, ::Vector{String}, ::Vector{String})` can be used to add population information to a `PopData` object.
Add population information to your lobster dataset.

Extract the sample names from the `.meta` dataframe of the lobsters `PopData`.
Then create a corresponding vector of population names from the sample names vector,

* Q 1.1.4: How many distinct populations are there?

Use your two vectors with the `populations!` function to set the population names of your `PopData`.

Now let's do some data exploration to see what kind of data we have in a typical `PopData` object.
The `.loci` field of the `PopData` contains the locus information for the lobsters.

Focus only on that same locus "un_1002665". You can do this by figuring out how to subset your dataframe, or by extracting the relevant columns of the dataframe and filtering them.

* Q 1.1.5: What distinct genotypes does the loci "un_1002665" have?

One of the genotypes is called `(1, 4)`, meaning one allele of type 1 and one allele of type 4, i.e. a lobster heterozygous at that loci.

Group the subset of the dataframe by the different populations using the `groupby` function of the package `DataFrames`.

* Q 1.1.6: Does the genotype (1, 4) of this locus exist in every population sampled?

In general, you can use the `groupby` function on column X to answer questions about the different values of X - for example, if you wanted to know which sample was most homozygous, you could group by sample and compute homozygosity for each group.

Let's clean up the data a bit.
When dealing with data, missing values are a constant source of annoyance.
Julia has a quite principled way of dealing with missing values, using the value `missing`.
This value has special behaviour:

```julia
julia> missing == 1
missing

julia> missing == missing
missing
```

The reason for this behaviour is that if a value is missing, it's unknown whether its true value is equal to 1 or not, so operations like `==` can neither return true, nor false.
I personally hate this behaviour with a fiery passion, but then again, I'm not a data scientist.

* Q 1.1.7: Find out how to check if some object `x` is `missing`. How?

Time to clean the data a little.
For example, can you throw any useless data away?
If there are any loci without genetic diversity, i.e. where there is only 1 genotype (plus perhaps a `missing` genotype), then that locus gives no information and can be removed.

For each genotype, count the number of distinct genotypes.
Remember to not consider `missing` when you count!

First, try to simply count _all_ the distinct genotypes in the entire dataset across all loci.
Yeah, this doesn't really make biological sense (because two genotypes at two different loci both called `(1, 4)` refers to different alleles), but do it just to get a hang of it.

* Q 1.1.8: How many distinct values are there in the "genotype" column of the `.loci` dataframe, not counting `missing`?

Figure out how to remove `missing` from a vector of values such as `[5, 2, missing, 19]`.

Now do the same counting, but for each locus individually.

* Q 1.1.9: How many loci has only 1 distinct genotype (plus possibly `missing`)?

Figure out how to remove all those loci from the `PopData`.

That missing data in the `.loci` dataframe sure is annoying.
It would be nice if we could remove all loci and samples where any genotype was missing.
Of course, you could risk having to remove too much data, so you need to figure out how much data would be removed if you removed all samples/loci with at least one missing genotype.

To count the number of samples with at least one missing genotype, 

* Q 1.1.10: How many samples have at least one genotype as `missing`?
* Q 1.1.11: How many different _genotypes_ are missing for at least one sample?
* Q 1.1.12: Make a judgement call: Does it make sense to remove samples and/or loci with at least one missing value?

### Exercise 1.2: Allele frequencies
The function `allele_freqtable` gives the allele frequencies.
Get the allele frequencies of all loci across all samples - that is - just of the entire `PopData` object.
Explore the returned object.

* Q 1.2.1. What type is the returned value? Explain its layout a little.

From just looking at the value in the REPL, it looks like most loci has only 2 different alleles.

* Q 1.2.2: Do ALL loci have 2 different alleles only?

The next few exercises take quite a bit of memory (PopGen is not very optimized...), so you will need to restrict your data to not risk running out of memory.

Create a new `PopData` object which only contains data from lobsters of the populations "LOB" and "SEA".
* Q 1.2.3: How many lobsters and loci do you have in this subset?

Let's do a crude, homemade population analysis.
The function `pairwise_identical` will create a table with the fraction of loci that have identical genotypes for each sample-sample pair.
Run it on your subset of lobsters from only the two aforementioned populations.
Like before, read the documentation of `pairwises_identical` and look at the output type and figure out what information it contains.
The `n` column is not that well documented: I think it's the number of loci used for the given row. 

Compute three mean values of pairwise identity:
    - Lobster pairs where both are from the LOB population
    - Pairs where both are from the SEA population
    - Pairs where they are from different populations

There is no built-in function to do this in PopGen, as far as I'm aware.
You can use the following sub-steps:
- How to you loop over the rows in a `DataFrame`?
- When iterating over the rows, given a row, how can you tell which of the three categories (LOB/LOB, SEA/SEA, LOB/SEA) the row represents, in order to put the `.identical` value into a corresponding vector?
- Given a vector of floats in the interval `[0.0, 1.0]`, how do you compute the mean value of it?

* Q 1.2.4: What is the mean value of pairwise identity for each of the three groups LOB/LOB, SEA/SEA and LOB/SEA?

* Q 1.2.5: Based on this simple measure, are the populations similar? That is, is the between- and within-population identities similar?

### Exercise 1.3: Check for siblingship
Let's do a more serious job of analysing these populations.
First though, we have to get rid of closely related lobsters from our dataset.
If we accidentally included a bunch of siblings in a population, our sample will be much more homogenous than the real population actually is.

We can check for relatedness by checking whether some individuals have more alleles in common than we would expect by random chance.
There are many statistical methods of doing this.
I don't understand the difference between them (I'm not a population geneticist), so let's just pick the Lynch-Li method.

This analysis is slow, so let's restrict out dataset even further to only the "LOB" population and a random selection of 250 loci.

To avoid including loci not present in the "LOB" population, first subset your dataframe to only the LOB population. After that, you can draw the random loci.

* Q 1.3.1: How do you sample 250 values from a vector without replacement?

Subset your `PopData` further to only include 50 loci chosen at random.
You can now calculate relatedness of your data.
Use the `relatedness` function on your sub-population with the Lynch-Li method, without using bootstrap (i.e. 1 iteration only, as is the default).

Inspect the result.

* Q 1.3.2: What is the type and layout of the result?

Have a look at the relatedness column (`LynchLi`).
To plot it, you need to make a copy of the column without missing values

* Q 1.3.3: How do you create a copy of a vector with values removed?

Use the package `UnicodePlots` to plot a histogram of the relatedness, without the missing values.

One of the problems of relatedness analysis is that they're not well-calibrated.
That is, while we expect a half-sibling to have relatedness of 0.25, due to the random nature of inheritance of any one locus, any observed half-sibling pair might have more or less than 25% of their alleles in common relative to a background distribution.

* Q 1.3.4: Simply looking at the distribution of relatedness, do you think it's likely there's a problem with too many siblings in your data?

To be sure, let's simulate some siblings and do a comparison between the observed data and simulated data.

The package PopGenSims.jl has a function `simulate_sibship` for this use-case.
Use it to simulate at least 500 unrelated, 500 half-sibling and 500 sibling lobster pairs from your sub-population.

Inspect the resulting `PopData` object.
Note the names of the populations of the simulated lobsters.
These represent unrelated, full siblings and half-siblings.

Use `exclude` to create three subgroups from your simulated data:
A) The simulated unrelated population
B) The simulated half-siblings
C) The simulated full-siblings

Then, for each group, calculate relatedness using the same method as you did for your original non-simulated population, and plot a histogram of it.
Hint: To show a histogram in a loop, call `display` on the histogram.

Compare the distribution of relatedness across the three groups with your data.

* Q 1.3.5: Based on these distributions of relatedness, does your original lobster population have a problem with too many siblings?

### Exercise 1.4: Fst and HWE measures
In an earlier question, you checked pairwise identity between two populations.
Now, let's try to do a better job of figuring out whether the populations are diverged by estimating the index of fixation Fst.

Use `pairwise_fst` on your whole dataset with all the lobsters (not the subset)!
Let's consider populations diverged if the Fst is above 0.125.

* Q 1.4.1: What Fst values do you get?
Compare it to your home-made population analysis in the earlier question comparing the LOB and SEA populations.
Does the Fst analysis broadly agree?

It could be that these lobsters just mate completely randomly with each other, thereby erasing any kind of population signal.
Maybe there are is no population structure at all in this dataset.
A more sensitive test of non-random mating is a test for Hardy-Weinberg equilibrium (HWE).

Search the PopGen docs to find out how to do a chi-square test for HWE for each locus.
If mating was random, you would expect the P value of the test for each locus to be randomly distributed in [0,1].

* Q 1.4.2: Is mating random in the total population, or are the P-values clearly skewed towards 0?
