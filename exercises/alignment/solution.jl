using BioAlignments, BioSequences, FASTX, CodecZlib

# Q 4.1.1
goose = aa"MEKIVLLLAIVSLVKSDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKTHNGKLCDL
NGVKPLILRDCSVAGWLLGNPMCDEFINVPEWSYIVEKASPANDLCYPGDFNDYEELKHL
LSRTNHFEKIQIIPKSSWSNHDASSGVSSACPYHGRSSFFRNVVWLIKKNSAYPTIKRSY
NNTNQEDLLVLWGIHHPNDAAEQTKLYQNPTTYISVGTSTLNQRLVPEIATRPKVNGQSG
RMEFFWTILKPNDAINFESNGNFIAPEYAYKIVKKGDSAIMKSELEYGNCNTKCQTPMGA
INSSMPFHNIHPLTIGECPKYVKSNRLVLATGLRNTPQRERRRKKRGLFGAIAGFIEGGW
QGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGVTNKVNSIIDKMNTQFEAVGREFNNLE
RRIENLNKQMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKVRLQLRDNAKELG
NGCFEFYHKCDNECMESVKNGTYDYPQYSEEARLNREEISGVKLESMGTYQILSIYSTVA
SSLALAIMVAGLSLWMCSNGSLQCRICI"

brevig = aa"MEARLLVLLCAFAATNADTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCK
LKGIAPLQLGKCNIAGWLLGNPECDLLLTASSWSYIVETSNSENGTCYPGDFIDYEELRE
QLSSVSSFEKFEIFPKTSSWPNHETTKGVTAACSYAGASSFYRNLLWLTKKGSSYPKLSK
SYVNNKGKEVLVLWGVHHPPTGTDQQSLYQNADAYVSVGSSKYNRRFTPEIAARPKVRDQ
AGRMNYYWTLLEPGDTITFEATGNLIAPWYAFALNRGSGSGIITSDAPVHDCNTKCQTPH
GAINSSLPFQNIHPVTIGECPKYVRSTKLRMATGLRNIPSIQSRGLFGAIAGFIEGGWTG
MIDGWYGYHHQNEQGSGYAADQKSTQNAIDGITNKVNSVIEKMNTQFTAVGKEFNNLERR
IENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVRNLYEKVKSQLKNNAKEIGNG
CFEFYHKCDDACMESVRNGTYDYPKYSEESKLNREEIDGVKLESMGVYQILAIYSTVASS
LVLLVSLGAISFWMCSNGSLQCRICI"

function custom_align(a::AminoAcidSeq, b::AminoAcidSeq)
    pairalign(
        GlobalAlignment(),
        a, b,
        AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-2)
    )
end

aln = custom_align(goose, brevig).aln
println(count_insertions(aln))
println(count_deletions(aln))
println(count_insertions(aln) + count_deletions(aln))

println(count_matches(aln))
println(count_mismatches(aln))

# For amino acid ID, it can be calculated in a variety of ways depending on
# whether you include gaps. I'll simply take the number of matches
# divided by the length of the larger sequence

id = count_matches(aln) / max(length(goose), length(brevig))
println(id)

# Q 4.1.2:
println(cigar(aln.a.aln))

# Q 4.1.3:
# There is a clear gap around the sequence RRRK, which is missing from
# the A/Brevik Mission/1/1918 virus

# Challenge version - one possible version:
# the version you come up with depends on your prior knowledge of
# the problem. Just identifying a stretch of missing amino acids from one
# of the two seqeucnes could work!
function has_polybasic_cleavage_site(aaseq::AminoAcidSeq)
    aln = custom_align(aaseq, aa"LATGLRNSPLREKRRKRGLFGAIAGFIEGGW").aln
    site = aa""

    # Polybasic cleavage site here goes from PL motive to GLF.
    for refpos in 9:20
        (seqpos, operation) = ref2seq(aln, refpos)
        operation ∈ (OP_SEQ_MATCH, OP_SEQ_MISMATCH) && push!(site, aaseq[seqpos])
    end

    # Check if we got mandatory GLF motif in our site,
    # and only return true if there are >3 basic sites as per recommendations
    # of the animal influenza network OFFLU
    return site[end-2:end] == aa"GLF" && count(∈((AA_R, AA_K)), site) > 3
end

# Q 4.2.1
hiv_records = open(collect, FASTA.Reader, "hiv.faa")
seqs = map(i -> ungap!(FASTA.sequence(i)), hiv_records)
aln = custom_align(seqs...).aln
seqa, seqb = LongAminoAcidSeq(first.(aln)), LongAminoAcidSeq(last.(aln))
open(FASTA.Writer, "/tmp/hiv.aln.faa") do writer
    bytes = sum(zip(hiv_records, (seqa, seqb))) do (r, s)
        write(writer, FASTA.Record(FASTA.identifier(r), s))
    end
    println(bytes)
end

# Q 4.3.1
reader = FASTA.Reader(GzipDecompressorStream(open("../../data/ns.fna.gz")))
ns_records = collect(reader)
close(reader)

ref = ns_records[findfirst(i -> FASTA.identifier(i) == "KM070129", ns_records)]
ref = FASTA.sequence(ref)

# Choose any alignment model
function custom_align(a::NucleotideSeq, b::NucleotideSeq)
    pairalign(
        GlobalAlignment(), a, b,
        AffineGapScoreModel(EDNAFULL, gap_open=-12, gap_extend=-2)
    )
end

sum(ns_records) do record
    seq = FASTA.sequence(record)
    aln = custom_align(seq, ref).aln
    (seqpos, op) = ref2seq(aln, 151)
    op ∈ (OP_SEQ_MATCH, OP_SEQ_MISMATCH) || return false
    seq[seqpos] == DNA_C
end |>  println

# It's not something to worry about. The fact that it's near-universal
# already, it just means that mutation is fixed in the population of
# influenza and now just part of how influenza works.

# Q 4.3.2:
dists = map(ns_records) do record
    seq = FASTA.sequence(record)
    aln = custom_align(seq, ref)
    (FASTA.identifier(record), aln.value)
end
println(first(sort!(dists, by=last)[end-1]))