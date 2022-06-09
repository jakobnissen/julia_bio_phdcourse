using JSON3, Kmers, BioSequences

# Q 2.1.1
bases = (RNA_A, RNA_C, RNA_G, RNA_U)
codons = RNACodon[]
for i in bases, j in bases, k in bases
    push!(codons, RNACodon((i, j, k)))
end
println(typeof(codons))

# Q 2.1.2
filtered = filter(codons) do codon
    codon ≤ reverse_complement(codon)
end
println(length(filtered))

# Q 2.1.3
Dict(
    c => translate(c)
    for c in codons
)

# Q 2.2.1
chrom = LongDNA{4}(open(JSON3.read, "../../data/sequences.json")[:chromy])

function count_palindromic_1(seq)
    sum(EveryKmer{DNAKmer{10}}(seq)) do (n, kmer)
        kmer == reverse_complement(kmer)
    end
end

println(count_palindromic_1(chrom))

# Q 2.2.2

# Copy-pasted from the biosequences solution with slight
# modifications
function count_palindromic_2(seq)
    n_palindromes = 0
    for i in 1:length(seq)-9
        subseq = seq[i:i+9]
        DNA_N in subseq && continue
        n_palindromes += subseq == reverse_complement(subseq)
    end
    return n_palindromes
end

# Run both to not time compilation
count_palindromic_1(chrom)
count_palindromic_2(chrom)

println(
    "Kmer solution:         ", (@timed count_palindromic_1(chrom)).time, '\n',
    "LongSequence solution: ", (@timed count_palindromic_2(chrom)).time
)

# Q 2.2.3
function count_alu_1(seq)
    n_alu = 0
    #alu = mer"TGGCTCACGCC"
    alu = DNAKmer{11}(dna"TGGCTCACGCC")
    for (i, kmer) in EveryKmer{typeof(alu)}(seq)
        if kmer == alu || reverse_complement(kmer) == alu
            n_alu += 1
        end
    end
    return n_alu
end

# Copy-pasted from the biosequences solution with slight
# modifications
function count_alu_2(seq)
    alu = dna"TGGCTCACGCC"
    n_alu = 0
    len_alu = length(alu)
    for seq in (seq, reverse_complement(seq))
        for i in 1:length(seq)-len_alu+1
            n_alu += seq[i:i+len_alu-1] == alu
        end
    end
    return n_alu
end

# Run both to not time compilation
count_alu_1(chrom)
count_alu_2(chrom)

println(
    "Kmer solution:         ", (@timed count_alu_1(chrom)).time, '\n',
    "LongSequence solution: ", (@timed count_alu_2(chrom)).time
)
