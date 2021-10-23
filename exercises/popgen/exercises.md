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

The method `populations!(::PopData, ::Vector{String}, ::Vector{String})` can be used to add population information to a `PopData` object.
Add population information to your lobster dataset.

* Q 1.1.3: How many distinct populations are there?

Now let's do some data wrangling.
The `.loci` field of the `PopData` contains the locus information for the lobsters.

Focus only on that same locus "un_1002665", by subsetting your dataframe.

* Q 1.1.4: What distinct genotypes does the loci "un_1002665" have?

One of the genotypes is called `(1, 4)`, meaning one allele of type 1 and one allele of type 4, i.e. a lobster heterozygous at that loci.

Group the subset of the dataframe by the different populations using the `groupby` function of the package `DataFrames`.

* Q 1.1.5: Does the genotype (1, 4) of this locus exist in every population sampled?

Let's clean up the data a bit.
For example, can we throw any data away?

Check if there are any loci without genetic diversity, i.e. where there is only 1 genotype (plus perhaps a `missing` genotype).

Remove any such loci.

* Q 1.1.6: How many loci were removed?

That missing data in the `.loci` dataframe sure is annoying.
It would be nice if we could remove loci or samples where any genotype was missing.
Of course, you could risk having to remove too much data, so you need to figure out how much data would be removed if you removed all samples/loci with at least one missing genotype.

Here, you can use `groupby again`

* Q 1.1.7: How many samples/loci would we need to remove?

### Exercise 1.2: Allele frequencies
The function `allele_freqtable` gives the allele frequencies.
Get the allele frequencies of all loci across all samples - that is - just of the entire `PopData` object.

* Q 1.2.1: Which loci has the most evenly distributed allele frequency (i.e. closest to 50-50), and which the most lopsided (closest to 0-100)?

Let's do a crude, homemade population analysis.
The function `pairwise_identical` will create a table with the fraction of loci that have identical genotypes for each sample-sample pair (it takes some time to compute).
For each of the N populations, calculate the mean identity between samples of those populations. Put it in an NxN matrix.

* Q 1.2.2: Based on this simple measure, are the populations similar? That is, is the between- and within-population identities similar?

### Exercise 1.3: Check for siblingship
Let's do a more serious job of analysing these populations.
First though, we have to get rid of closely related lobsters from our dataset.
If we accidentally included a bunch of siblings in a population, our sample will be much more homogenous than the real population actually is.

We can check for relatedness by checking whether some individuals have more alleles in common than we would expect by random chance.
There are many statistical methods of doing this.
I don't understand the difference between them, so let's just pick the Lynch-Li method.

This analysis is slow, so let's restrict out dataset:
Create a `PopObj` with only the population `TRI`, and a random selection of 250 loci. Call it `subpop`

You can now calculate relatedness of your data.
Use the `relatedness` function on your sub-population with the Lynch-Li method and 100 bootstrap iterations.

Inspect the dataframe.
Siblings have relatedness of 0.5, and half-siblings of 0.25.

Q 1.3.1: What is the distribution of mean relatedness (`LynchLi_mean`)? Could there be a problems with siblings in the data?

One of the problems of relatedness analysis is that they're not well-calibrated.
That is, while we expect a half-sibling to have relatedness of 0.25, due to the random nature of inheritance, any observed half-sibling pair might have more or less than 25% of their alleles in common relative to a background distribution.

Let's simulate some siblings and do a comparison between the observed data and simulated data.

Use `simulate_sibship` to simulate at least 500 unrelated, 500 half-sibling and 500 sibling lobsters from your sub-population.
Then compute relatedness for each population using the same method (Lynch-Li).

Compare the distribution of `LynchLi` values between:
* Your observed subpopulation
* The simulated unrelated population
* The simulated half-siblings
* The simulated full-siblings

* Q 1.3.2: Based on these distributions, are there still a problem with siblingship in your subpopulation?

### Exercise 1.4: Fst and HWE measures
It looks like we're OK w.r.t siblingship.

Compared to your home-made test in Q 1.2.2, the Fst measure is a more rigorous test for within- vs between population differences.

Do a `pairwise_fst` test on your whole data (not the subset)!

* Q 1.4.1: What values do you get? Compare it to your home-made population analysis. Do they agree?

It could be that these lobsters just mate completely randomly with each other, thereby erasing any kind of population signal.
Maybe there are is no population structure at all in this dataset.
A more sensitive test of non-random mating is a test for Hardy-Weinberg equilibrium (HWE).

Use PopGen.jl to do a chi-square test for HWE for each locus across all samples. Then do HWE tests by population.

If mating was random, you would expect the P value to be randomly distributed in [0,1].

* Q 1.4.2: Is mating random in the total population, or are the P-values clearly skewed towards 0?

* Q 1.4.3: How does it look when doing the analysis by population? Can you conclude anything from this?
