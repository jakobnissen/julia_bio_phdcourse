using JSON3, BioSequences

# Q 1.1.1
data = open(JSON3.read, "sequences.json")
println(typeof(data[:chromy]))
println(typeof(data[:has][1][2]))

# The DNA sequences are String

data = Dict(
    :has => [(n, LongDNA{4}(s)) for (n,s) in data[:has]],
    :chromy => LongDNA{4}(data[:chromy])
)

# Q 1.1.2
println(length(data[:has]))

# Q 1.1.3
println(sum(length(seq) for (name, seq) in data[:has]))

# Q 1.2.1
for symbol in (DNA_W, DNA_Y, DNA_C, DNA_N)
    println(count(isequal(symbol), data[:chromy]))
end

# Q 1.2.2

# This solution takes a few seconds to run. It can be solved more efficiently
# using kmers, but whatever.
n_palindromes = 0
seq = data[:chromy]
for i in 1:length(seq)-9
    global n_palindromes
    subseq = seq[i:i+9]
    DNA_N in subseq && continue
    n_palindromes += subseq == reverse_complement(subseq)
end

println(n_palindromes)

# Q 1.2.3
alu = dna"TGGCTCACGCC"
n_alu = 0
len_alu = length(alu)
for seq in (data[:chromy], reverse_complement(data[:chromy]))
    for i in 1:length(seq)-len_alu+1
        global n_alu
        n_alu += seq[i:i+len_alu-1] == alu
    end
end

println(n_alu)

# Q 1.3.1
count(data[:has]) do (name, seq)
    iszero(count(isambiguous, seq))
end |> println

# Q 1.3.2
filter(data[:has]) do (name, seq)
    seq[1:3] != dna"AGC"
end |> only |> first |> println

# Q 1.4.1
println(translate(rna"AUGCGCGGAGAUAUAUUA"))

# Q 1.4.2
seq = dna"ATGATCTCTATGGTATGGAATATGTGTAACGTTCGTGCAGATGCACCTATGGCGTGGCAGAAATTGTTTC
AAGACCCTGCAACGTCCAACATGGAAGGAATTGTAGACCTTCACGCAGATATTTGCTTTTTTTTGATTGT
AATTTTGGTTTTGGTTTTGTGGCTCGGAGCTCGAATTCTTGTAAGTTTTCATCACAACTTGCAACCTGTA
CCAGAGCGTTTTAACCACCACACTTCTTTGGAATTGGTTTGGGCAATTTTGCCTTCTGTTATTGTAACAT
TGATTGCATTGCCTTCCTAGTCTCTTGTGTATACTTATGATGATTAGGTAAGTAAACCCGCATTGACCGT
GAAGGTTACTGGGCGTCAGTGGTATTGGAGCTACGCTATGAATGAGCATGTACAAATGAATTTGAGTCAA
CAGGCGAAAGATTTGTTGCTTCAGTCTTCA"

findfirst(
    isequal(AA_Term),
    translate(seq, code=BioSequences.vertebrate_mitochondrial_genetic_code)
)
# Since there is a stop at position 100, this cannot be from a
# vertebrate mitochondrion.
# In fact, it is from Scenedesmus obliquus, a single-celled algae.

# Q 1.4.3
seq = filter(data[:has]) do (name, seq)
    name == "AY650270"
end |> only |> last
coding = seq[22:1701]

protein = translate(coding)
n_mutations = 0
for i in eachindex(protein)
    codon = coding[3i-2:3i]
    for j in 1:3
        current_base = codon[j]
        for base in (DNA_A, DNA_C, DNA_G, DNA_T)
            global n_mutations
            base == current_base && continue
            cp = copy(codon)
            cp[j] = base
            n_mutations += only(translate(cp)) == protein[i]
        end
    end
end
println(n_mutations)

# Q 1.4.4
seq = aa"RTHLMPYDQDISYIMDHKGP"
n_codons = Dict{AminoAcid, Int}()
bases = (RNA_A, RNA_C, RNA_G, RNA_U)
for i in bases, j in bases, k in bases
    codon = LongRNA{2}((i, j, k))
    aa = only(translate(codon))
    n_codons[aa] = get(n_codons, aa, 0) + 1
end
println(prod(i -> n_codons[i], seq))
