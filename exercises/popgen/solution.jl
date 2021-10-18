using PopGen, GeneticVariation, GZip, DataFrames, UnicodePlots, PopGenSims

lobsters = vcf("lobsters.vcf.gz")

# Q 1.1.1: 8411 loci and 586 samples
# the display of the lobsters value shows it.

samplenames = collect(samples(lobsters))
popnames = map(samplenames) do samplename
    String(first(split(samplename, '_')))
end
populations!(lobsters, samplenames, popnames)

# Q 1.1.2:
n_pops = length(unique(popnames))
println(n_pops)

# Q 1.1.3:
# How many distinct genotypes does the loci "un_1002665" have?
subdf = lobsters.loci[lobsters.loci.locus .== "un_1002665", :]
unique(subdf.genotype)
# one is missing, so 3

# Q 1.1.4:
mask = map(subdf.genotype) do genotype
    genotype !== missing && genotype == (1, 4)
end
println(length(unique(subdf[mask, :].population)) == n_pops)

# Q 1.1.5:
groups = groupby(lobsters.loci, :locus)
filtered_groups = filter(groups) do group
    genotypes = Set(group.genotype)
    delete!(genotypes, missing)
    length(genotypes) < 2
end
println(length(filtered_groups))

# Q 1.1.6
for sym in (:locus, :name)
    groups = groupby(lobsters.loci, sym)
    n_groups = length(groups)
    has_missing = count(groups) do group
        any(ismissing, group.genotype)
    end
    println(sym, " ", has_missing, " ", n_groups)
end

# Q 1.2.1:
freq = allele_freqtable(lobsters)
byfreq = sort!(collect(groupby(freq, :locus)), by=i -> minimum(i.frequency))
for f in (first, last)
    println(f(byfreq).locus[1])
end

# Q 1.2.2
pairwises = pairwise_identical(lobsters)
ids = Dict{Tuple{String, String}, Vector{Float64}}()
for r in eachrow(pairwises)
    pop1 = first(split(r.sample_1, '_'))
    pop2 = first(split(r.sample_2, '_'))
    push!(get!(Vector{Float64}, ids, (pop1, pop2)), r.identical)
end
indexof = Dict(x => i for (i, x) in enumerate(populations(lobsters)))
heatmap = fill(1.0, (length(indexof), length(indexof)))
for ((pop1, pop2), v) in ids
    mean_id = sum(v) / length(v)
    index1, index2 = indexof[pop1], indexof[pop2]
    heatmap[index1, index2] = mean_id
    heatmap[index2, index1] = mean_id
end

# Yes, they are very similar

# Q 1.3.1:
select_loci = rand(unique(lobsters.loci.locus), 250)
subpop = keep(keep(lobsters, locus=select_loci), population="TRI")

rel = relatedness(subpop, method = LynchLi, iterations = 100).LynchLi
histogram(collect(skipmissing(rel.LynchLi)))

# It's distributed around -1, 0.4. Clearly this means there 
# are huge error bars for this estimate. The ones at 0.4
# could be half- or full siblings, hard to know.

# Q 1.3.2:
sims = simulate_sibship(subpop; fullsib=500, halfsib=500, unrelated=500)
simsrel = relatedness(sims, method=LynchLi)

for kind in ["fullsib", "halfsib", "unrelated"]
    mask = map(s -> occursin(kind, s), simsrel.sample_1)
    v = collect(skipmissing(simsrel.LynchLi[mask]))
    println(kind)
    display(histogram(v))
end

# The observed population has even lower siblingship than the
# randomly simulated individuals. That's, uh, unexpected,
# but at least it looks like there are no problem with
# siblingship in general. There still could be individuals which are
# siblings, but it's hard to tell.

# Q 1.4.1
fst = pairwise_fst(lobsters)
println(maximum(Matrix(fst.results)))

# Yeah, both approaches show the populations are very similar

# Q 1.4.2
hwe_all = hwe_test(lobsters)
histogram(collect(skipmissing(hwe_all.P)))

# It's clearly skewed towards 0, by a lot

# Q 1.4.3
hwe_pop = hwe_test(lobsters; by="population")
histogram(collect(skipmissing(hwe_pop.P)))

# It's much less skewed! But, it might just be higher P-values
# because there are fewer samples.
# I would need to correct for sample size to say anything here...