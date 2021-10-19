using XAM, FASTX

# The following solutions were obtained by simply readin the 
# SAM file:

# Q 5.1.1:
# Identifier: NZ_CP034181.1
# Sequence length: 5067231

# Q 5.1.2:
# 2209 = 1 + 32 + 128 + 2048
# So:
#   Read is mapped
#   Next read in read group reversed
#   Read is last in group
#   Alignment is supplementary

# Q 5.1.3:
# 65H35M

# Q 5.1.4:
# 1 - the very first position in the genome

# Q 5.1.5:
# Bacterial genomes are circular. Hence, the read really ought to map
# to both the beginning and the end of our reference.
# However, it seems the aligner used (BWA) cannot map to circular references
# by default.

# Q 5.2.1:
records = open(collect, BAM.Reader, "../../data/subset.bam")
println(length(records))

# Q 5.2.2:
filter!(records) do record
    BAM.ismapped(record) && BAM.isnextmapped(record) && BAM.isprimary(record)
end

println(100000 - length(records))

# Q 5.2.3:
depth = sum(BAM.seqlength, records) / BAM.position(last(records))
println(depth)

# Yes, this is quite similar to the estimated depth from the FASTX
# exercise.

# Q 5.2.4:
open(FASTQ.Writer, "/tmp/delme.fq") do writer
    for record in records
        write(writer, FASTQ.Record(
            BAM.tempname(record),
            BAM.sequence(record),
            BAM.quality(record)
        ))
    end
end
