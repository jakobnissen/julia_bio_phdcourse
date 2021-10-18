using FASTX, BioSequences, CodecZlib

reader = FASTA.Reader(GzipDecompressorStream(open("ns.fna.gz")))
ns_records = collect(reader)
close(reader)

# Q 3.1.1
println(FASTA.identifier(ns_records[end-4]))
println(FASTA.description(ns_records[end-4]))

# Q 3.1.2
count(ns_records) do record
    !isnothing(FASTA.identifier(record)) &&
    !isnothing(FASTA.description(record))
end |> println

# Q 3.1.3
mean_len = sum(i -> length(FASTA.sequence(i)), ns_records) / length(ns_records)
println(count(i -> length(FASTA.sequence(i)) > mean_len, ns_records))

# Q 3.1.4
seqpairs = map(ns_records) do record
    id = let
        i = FASTA.identifier(record)
        isnothing(i) ? "" : i
    end
    (id, FASTA.sequence(LongDNASeq, record))
end

# Q 3.1.5
open(FASTA.Writer, "/tmp/ns.revcomp.fna") do writer
    for (id, seq) in seqpairs
        write(writer, FASTA.Record(id, reverse_complement(seq)))
    end
end

# Q 3.2.1
# Low-effort solution: Load it into memory
reader = FASTQ.Reader(GzipDecompressorStream(open("M_abscessus.fq.gz")))
m_abs_records = collect(reader)
close(reader)

lens = map(m_abs_records) do record
    length(record.sequence)
end

println(sum(lens) / length(lens))
println(sort(lens)[div(length(lens), 2)])

# Q 3.2.2:
println(sum(lens))

# Q 3.2.3:
n_unambiguous = count(m_abs_records) do record
    iszero(count(isambiguous, FASTA.sequence(LongDNASeq, record)))
end
println(n_unambiguous / length(lens))

# Q 3.3.1
# This is the same as counting the occurrence of the 3-mer AAA
# in the first N-18 bases of each read
n_aaa_kmers = 0
for record in m_abs_records
    seq = FASTA.sequence(LongDNASeq, record)
    for kmerres in each(DNAMer{3}, seq[1:end-21+3])
        n_aaa_kmers += kmerres.fw == mer"AAA" || kmerres.bw == mer"AAA"
    end
end
println(n_aaa_kmers)

# Q 3.3.2
counts = Dict{DNAMer{21}, Int}()
for record in m_abs_records
    seq = FASTA.sequence(LongDNASeq, record)
    for kmerres in each(DNAMer{21}, seq)
        canon = canonical(kmerres)
        counts[canon] = get(counts, canon, 0) + 1
    end
end
println(length(counts))
println(sum(values(counts)))

# Q 3.3.3:
spectrum = zeros(Int, maximum(values(counts)))
for v in values(counts)
    spectrum[v] += 1
end
# using UnicodePlots; lineplot(spectrum[10:100])
# looks like peak is at 31
# and the left peak is at 1:10

# Q 3.3.4:
# DNA sequencers are not perfectly accurate.
# Sequencing errors introduce spurious kmers.

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