using FASTX, BioSequences, Kmers, CodecZlib

# Q 3.1.1
reader = FASTA.Reader(
    GzipDecompressorStream(open("../../data/ns.fna.gz"))
)
ns_records = collect(reader)
close(reader)

println(
    FASTA.identifier(ns_records[end-4]), '\n',
    FASTA.description(ns_records[end-4])
)

# Q 3.1.2
count(ns_records) do record
    FASTA.hasidentifier(record) &&
    FASTA.hasdescription(record)
end |> println

# Q 3.1.3
mean_len = sum(length(FASTA.sequence(rec)) for rec in ns_records) / length(ns_records)
println(count(length(FASTA.sequence(rec)) > mean_len for rec in  ns_records))

# Q 3.1.4
seqpairs = map(ns_records) do record
    id = let
        i = FASTA.identifier(record)
        isnothing(i) ? "" : i
    end
    (id, FASTA.sequence(LongDNA{4}, record))
end

# Q 3.1.5
open(FASTA.Writer, "/tmp/ns.revcomp.fna") do writer
    for (id, seq) in seqpairs
        write(writer, FASTA.Record(id, reverse_complement(seq)))
    end
end

# Q 3.2.1
# Low-effort solution: Load it all into memory
reader = FASTQ.Reader(
    GzipDecompressorStream(open("../../data/M_abscessus.fq.gz"))
)
m_abs_records = collect(reader)
close(reader)

lens = map(m_abs_records) do record
    length(record.sequence)
end

println("Mean length: ", sum(lens) / length(lens))
println("Median length: ", sort(lens)[div(length(lens), 2)])

# Q 3.2.2:
println(sum(lens))

# Q 3.2.3:
n_unambiguous = count(m_abs_records) do record
    iszero(count(isambiguous, FASTA.sequence(LongDNA{4}, record)))
end
println(n_unambiguous / length(lens))

# Q 3.3.1
# This is the same as counting the occurrence of the 3-mer AAA
# in the first N-18 bases of each read
function count_aaa_kmers(record)
    n = 0
    seq = FASTA.sequence(LongDNA{4}, record)
    for (_, kmer) in EveryKmer{Kmer{DNAAlphabet{4}, 21}}(seq)
        n += kmer[1:3] == mer"AAA"d
        n += reverse_complement(kmer)[1:3] == mer"AAA"d
    end
    return n
end

n_aaa_kmers = 0
for record in m_abs_records
    global n_aaa_kmers += count_aaa_kmers(record)
end
println(n_aaa_kmers)

# Q 3.3.2
function update_counts!(counts, record)
    seq = FASTA.sequence(LongDNA{4}, record)
    for (i, kmer) in EveryCanonicalKmer{DNAKmer{21}}(seq)
        counts[kmer] = get(counts, kmer, 0) + 1
    end
end

counts = Dict{DNAKmer{21, 1}, Int}()
for record in m_abs_records
    update_counts!(counts, record)
end

println(length(counts))
println(sum(values(counts)))

# Q 3.3.3:
spectrum = zeros(Int, maximum(values(counts)))
for v in values(counts)
    spectrum[v] += 1
end

# I haven't added a plotting package to the Project,
# because it would bloat it up.
# But I would do it like this:

# using UnicodePlots; lineplot(spectrum[10:100])
# looks like peak is at 31
# and the left peak is at 1:10

# Q 3.3.4:
# DNA sequencers are not perfectly accurate.
# Sequencing errors introduce spurious kmers,
# that only show up a few times.

# Q 3.3.5:
T = sum(lens)
Tk = sum((i*j) for (i,j) in enumerate(spectrum) if i > 10)
distinct_kmers = sum(spectrum[11:end])
P = Tk / distinct_kmers
D = (P*T)/Tk
println(D)

# Q 3.3.6:
S = T/D
println(S)

# M abscessus is around 5.08-5.09 Mbp, within 0.2% of the result
