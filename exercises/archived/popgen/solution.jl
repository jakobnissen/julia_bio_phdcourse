using PopGen

# * Q 1.1.1: Inspect the type by typing it into the terminal.'
# How many loci (SNP markers) and samples are there in the file?

# First, load in the file. Simply running `vcf("lobsters.vcf.gz")` tells me to also
# use the GZip package and GeneticVaritation, so I do that:

using GeneticVariation, GZip
lobsters = vcf("/home/jakob/Downloads/data/lobsters.vcf.gz")

# 8411 loci and 586 samples
# the display of the lobsters value shows it.

# * Q 1.1.2: What fields do the `PopData` have? What are their types?
dump(PopData)
# :meta field and :loci field, both of type DataFrame

# * Q 1.1.3: Given the string `"BON_01"` composed of a prefix, and underscore, and a suffix,'
# how do you extract the prefix?

# 
first(split("BON_01", '_'))

# * Q 1.1.4: How many distinct populations are there?
sample_names = lobsters.meta.name
population_names = [first(split(i, '_')) for i in lobsters.meta.name]
println(length(unique(population_names)))

# If you try doing this:
# populations!(lobsters, sample_names, population_names)
# It will fail with a MethodError, because `population_name` is a vector of SubString, not a vector of String.
# Here, the PopGen package should be more liberal in what it accepts, but we, as users, can't change that,
# so let's just convert our population_names to a vector of String
populations!(lobsters, sample_names, map(String, population_names))

# * Q 1.1.5: What distinct genotypes does the loci "un_1002665" have?
subset_df = lobsters.loci[lobsters.loci.locus .== "un_1002665", :]
println(unique(subset_df.genotype))

# Three genotypes plus a missing genotype.

# * Q 1.1.6: Does the genotype (1, 4) of this locus exist in every population sampled?
using DataFrames
by_pop = groupby(lobsters.loci, :population)
all(by_pop) do group
    (1, 4) in group.genotype
end

# Yes, looks like it.

# * Q 1.1.7: Find out how to check if some object `x` is `missing`.
x = missing
ismissing(x) # or x === missing, they are equivalent.


# * Q 1.1.8: How many distinct values are there in the "genotype" column of the `.loci` dataframe,
# not counting `missing`?
genotypes = unique(lobsters.loci.genotype)

# There are 10 genotypes + missing.
# To remove the missing value programatically, one could do
# `setdiff(genotypes, [missing])`
# or one could use `skipmissing`, or something fancy using `coalesce`.

# * Q 1.1.9: How many loci has only 1 distinct genotype (plus possibly `missing`)?
n_nondiverse_loci = 0
for grouped_by_locus in groupby(lobsters.loci, :locus)
    genotypes = setdiff(unique(grouped_by_locus.genotype), [missing])
    if length(genotypes) < 2 # can also be empty!
        n_nondiverse_loci += 1
    end
end
println(n_nondiverse_loci)

# There are no such loci, so I don't have to figure out how to remove them.

# * Q 1.1.10: How many samples have at least one genotype as `missing`?
# * Q 1.1.11: How many different _genotypes_ are missing for at least one sample?
for column in (:locus, :name)
    grouped = groupby(lobsters.loci, column)
    println("Column $column has $(length(grouped)) values total")
    n_has_one_missing = count(grouped) do group
        any(ismissing, group.genotype)
    end
    println("...and $n_has_one_missing has at least one missing value")
end

# * Q 1.1.12: No, we would remove everything except one locus.

# * Q 1.2.1. What type is the returned value? Explain its layout a little.
freqs = allele_freqtable(lobsters)

# It's a DataFrame with 16822 rows one for each locus-allele pair, and the following columns:
# :locus, :allele, :count, :frequency
# Noticably, there are 8411 loci and 2*8411 locus-allele pairs, meaning there are only two
# alleles per locus.

# * Q 1.2.2: Do ALL loci have 2 different alleles only?
all(groupby(freqs, :locus)) do locus
    nrow(locus) == 2
end

# * Q 1.2.3: Create a new `PopData` object which only contains data from lobsters
# of the populations "LOB" and "SEA".

# Use the function `exclude`.
all_pops = map(String, population_names)
pops_to_exclude = setdiff(unique(all_pops), ["LOB", "SEA"])
sublob = exclude(lobsters, population=pops_to_exclude)

# 69 lobsters, 8411 loci

# * Q 1.2.4: What is the mean value of pairwise identity for each of the three groups
# LOB/LOB, SEA/SEA and LOB/SEA?

pairwises = pairwise_identical(sublob)
identities = Dict()
for row in eachrow(pairwises)
    pop1 = first(split(row.sample_1, '_'))
    pop2 = first(split(row.sample_2, '_'))
    # sort the two pops to avoid having both (LOB, SEA) and (SEA, LOB) as distinct keys
    key = minmax(pop1, pop2)
    if haskey(identities, key)
        push!(identities[key], row.identical)
    else
        identities[key] = [row.identical]
    end
end

for (key, id_vector) in identities
    println(key, " ", sum(id_vector) / length(id_vector))
end

# * Q 1.2.5: Based on this simple measure, are the populations similar?
# That is, is the between- and within-population identities similar?

# Yes, they are very similar!

# * Q 1.3.1: How do you sample 250 values from a vector without replacement?

# Several different ways, but one is:
using Random
sample_no_replace(v::Vector, n::Integer) = shuffle(v)[1:n]

# * Q 1.3.1: Get a vector of the unique loci in your dataset, then draw 250 random values using `rand`.
all_pops = map(String, population_names)
pops_to_exclude = setdiff(unique(all_pops), ["LOB"])
sublob = exclude(lobsters, population=pops_to_exclude)
loci_to_exclude = setdiff(sublob.loci.locus, sample_no_replace(unique(sublob.loci.locus), 250))
sublob = exclude(sublob, locus=loci_to_exclude)
rel = relatedness(sublob, method = LynchLi)

# * Q 1.3.2: What is the type and layout of the result?

# It's a DataFrame with rows having the two sample names, n loci compared between the the samples,
# and a relatedness estimate for all lobster pairs.

# * Q 1.3.3: How do you create a copy of a vector with values removed?
remove_missing(v::Vector) = collect(skipmissing(v))

# * Q 1.3.4: Simply looking at the distribution of relatedness,
# do you think it's likely there's a problem with siblings in your data?

using UnicodePlots
histogram(remove_missing(rel.LynchLi))

# Probably not - relatedness is distributed around approx -0.2, which is actually less
# than expected from random mating!

# * Q 1.3.5: Based on these distributions of relatedness, does your original lobster population
# have a problem with too many siblings?
sims = simulate_sibship(sublob; fullsib=500, halfsib=500, unrelated=500)
for group in ("fullsib", "halfsib", "unrelated")
    pop = exclude(sims, population=setdiff(unique(sims.meta.population), [group]))
    rel_ = relatedness(pop, method = LynchLi)
    display(histogram(remove_missing(rel_.LynchLi)))
end

# No. In fact, the population is LESS related than simulated unrelated lobsters!

# * Q 1.4.1: What Fst values do you get?
# Compare it to your home-made population analysis in the earlier question
# comparing the LOB and SEA populations.
# Does the Fst analysis broadly agree?
fst = pairwise_fst(lobsters)

# Looking at it, all the values seem to be around 0, so not diverged.
# Let's find the LOB/SEA Fst to compare the two
lob_ind, sea_ind = [findfirst(isequal(i), names(fst.results)) for i in ("LOB", "SEA")]
fst.results[max(lob_ind, sea_ind), min(lob_ind, sea_ind)]

# Yes, it's near 0, which agreed with the earlier estimate that the two
# populations were nearly identical.

# Q 1.4.2
hwe_all = hwe_test(lobsters)
histogram(remove_missing(hwe_all.P))

# It's clearly skewed towards 0, by a lot
