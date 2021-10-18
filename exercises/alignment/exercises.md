# Exercise 4: BioAlignments.jl
This exercise is an introduction to pairwise alignment - multiple sequence alignment (MSA)'s smaller brother.
Unfortunately, we will not be going into MSA's in this course.

Besides the file in this directory, we will also re-use the files from earlier exercises.

### Exercise 4.1: Getting a taste of alignment
The _goose Guangdong virus_ was a strain of highly pathogenic avian influenza isolated in Guangdong, China in 1996.
The virus is considered an ancestor to the Z genotype of avian H5N1 influenza virus, which has since spread globally in birds.
It is capable of infecting humans and have a case fatality rate of approximately 50%.
There is considerable worry that a descendant of this virus will cause a future pandemic which could be catastrophic.
This is the hemagglutinin sequence of the original goose Guangdong virus:

>A/Goose/Guangdong/1/1996
MEKIVLLLAIVSLVKSDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKTHNGKLCDL
NGVKPLILRDCSVAGWLLGNPMCDEFINVPEWSYIVEKASPANDLCYPGDFNDYEELKHL
LSRTNHFEKIQIIPKSSWSNHDASSGVSSACPYHGRSSFFRNVVWLIKKNSAYPTIKRSY
NNTNQEDLLVLWGIHHPNDAAEQTKLYQNPTTYISVGTSTLNQRLVPEIATRPKVNGQSG
RMEFFWTILKPNDAINFESNGNFIAPEYAYKIVKKGDSAIMKSELEYGNCNTKCQTPMGA
INSSMPFHNIHPLTIGECPKYVKSNRLVLATGLRNTPQRERRRKKRGLFGAIAGFIEGGW
QGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGVTNKVNSIIDKMNTQFEAVGREFNNLE
RRIENLNKQMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKVRLQLRDNAKELG
NGCFEFYHKCDNECMESVKNGTYDYPQYSEEARLNREEISGVKLESMGTYQILSIYSTVA
SSLALAIMVAGLSLWMCSNGSLQCRICI

The virus isolate _A/Brevig Mission/1/1918_ was extracted from a human corpse preserved in the arctic permafrost.
The corpses were from the small village called Brevig Mission in USA, where the 1918 influenza pandemic (called the Spanish Flu) killed 72 out of 80 residents of the village, as well as 25 million other people across the world.
The sequence below is the hemagglutinin protein of the virus responsible for the
pandemic of 1918:

>A/Brevig Mission/1/1918
MEARLLVLLCAFAATNADTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCK
LKGIAPLQLGKCNIAGWLLGNPECDLLLTASSWSYIVETSNSENGTCYPGDFIDYEELRE
QLSSVSSFEKFEIFPKTSSWPNHETTKGVTAACSYAGASSFYRNLLWLTKKGSSYPKLSK
SYVNNKGKEVLVLWGVHHPPTGTDQQSLYQNADAYVSVGSSKYNRRFTPEIAARPKVRDQ
AGRMNYYWTLLEPGDTITFEATGNLIAPWYAFALNRGSGSGIITSDAPVHDCNTKCQTPH
GAINSSLPFQNIHPVTIGECPKYVRSTKLRMATGLRNIPSIQSRGLFGAIAGFIEGGWTG
MIDGWYGYHHQNEQGSGYAADQKSTQNAIDGITNKVNSVIEKMNTQFTAVGKEFNNLERR
IENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVRNLYEKVKSQLKNNAKEIGNG
CFEFYHKCDDACMESVRNGTYDYPKYSEESKLNREEIDGVKLESMGVYQILAIYSTVASS
LVLLVSLGAISFWMCSNGSLQCRICI

Use BioAlignments.jl to align the two sequences using Needleman-Wunch.
For parameters, choose an affine
gap score model using the BLOSUM62 substitution matrix, a gap opening score of -10
and a gap extension score of -2.

* Q 4.1.1: How many gaps are there in the alignment?
How many gaps in the goose Guangdong sequence, and how many in the Spanish Flu sequence?
How many matches and mismatches?
What is the amino acid identity % between the two proteins?

* Q 4.1.2: Find a way to create a CIGAR string from the alignment.

The defining characteristic of highly pathogenic avian influenzas (and a trait of some other viruses, such as SARS-CoV2) is their _polybasic cleavage site_, which enhances virulence, and in influenza allows the infection to spread systemically into for example the brain of the patient.

The goose Guangdong virus have such a polybasic cleavage site, but it is missing in human influenzas, including the Spanish Flu.

* Q 4.1.3: Looking at the alignment, identify the polybasic cleavage site.
Bonus challenge: See if you can do it programatically!

### Exercise 4.2: FASTX
The file `hiv.faa` contains two sequences: The envelop protein from human immunodeficiency virus-1 (HIV-1) and from the related simian immunodeficiency virus (SIV).
I half-hearted tried to align the two sequences by hand, but it's pretty hard, so I gave up and now the file is a mess.
Now it's your job to clean it up and align it properly.

Load in the sequences using `FASTX.jl`, remove the gaps in the sequences, then re-align them using Needleman-Wunch, then write the resulting alignment in FASTA format to a file `hiv.aln.faa`.

* Q 4.2.1: What is the precise file size of `hiv.aln.fna`?

### Exercise 4.3: Mutations
Go back to the `ns.fna.gz` file from the FASTX exercise.

Previous researchers have reported (PubMed ID: 18032512) that a T -> C mutation in a certain position of the NS1 gene of influenza is critical for virus to evade the innate host immune system.

This position is equivalent of position 151 in the sequence named "KM070129" in `ns.fna.gz`.
Use alignment to figure out which positions are equivalent to position 150 for each sequence.

* Q 4.3.1: How many of the representative set of sequences in `ns.fna.gz` has this mutation (i.e. `C` at that position)?
Given that the `ns.fna.gz` file contains a broad set of representative influenza sequences, is the spread of this mutation something to worry about?

You want to find which sequence in `ns.fna.gz` is most closely related to KM070129.
To do this properly, you would need to use actual phylogenetic software, but for this exercise, just use the alignment score.

* Q 4.3.2: Which sequence(s) in `ns.fna.gz` is closest to KM070129?
