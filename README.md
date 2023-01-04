# hope
Assessment of homopolymer performance of ONT sequencing

## Installation

hope is available as a precompiled binary in a release on this repo. Alternatively, you can clone this repo and build from source.

### Build from source

To build from source, you will need rust installed (https://www.rust-lang.org/tools/install)

clone this repo and execute the following command within the repo directory:

    cargo build --release
    
The built binary will then be present at the path `./target/release/hope`

## Usage

Information about required inputs and optional settings can be found by running `hope -h`.

    USAGE:
        hope [OPTIONS] --input-homos <INPUT_HOMOS> --assembly <ASSEMBLY> --bam <BAM> --outprefix <OUTPREFIX>

    OPTIONS:
        -a, --assembly <ASSEMBLY>          the input assembly file
        -b, --bam <BAM>                    the input bam file
        -c, --context                      include sequence context in outfile?
        -h, --help                         Print help information
        -i, --input-homos <INPUT_HOMOS>    file with homopolymer locations and bases
        -o, --outprefix <OUTPREFIX>        the outprefix
        -V, --version                      Print version information

Option details are described in the following section

### Options

#### BAM

BAM file must be located in the same directory as a corresponding index file (.bai)

As hope is aimed at identifying sequencing errors that result in incorrect homopolymer lengths in reads, it is important to map your reads using settings that tolerate indels. In our tests using `minimap2` we found that the options `-A 2 -B 10` result in a good read mapping for this analysis.

#### Assembly

Fasta format. Not gzipped

#### input-homos

Tab-delimited file describing the location of homopolymers to be assessed. Contigs must correspond to sequence IDs in the provided assembly and BAM file. 

Expected columns: Contig ID, start, stop, base, length

N.B., start and stop are 1-base coordinates

#### Context (Optional)

If set, the aligned sequence from 30 bases upstream to 30 bases downstream of the homopolymer will be included in the output file. In cases when a homopolymer is fewer than 30 bases from the end of the mapped portion of the read, all mapped bases will be returned.

#### Outprefix

Output file will be written to a path constructed by adding "out.txt" to whatever you provide.

## Output file

The columns in the output file are: homopolymer_length, homopolymer_base, difference, read_context, assembly_context, homo_start, read_ID

Note that if `-c` is not used, read_context and assembly_context will not be present.

Each entry in the output file corresponds to the information in one read at one homopolymer position. Data are sorted by reads. The data for each read is sorted by homopolymer start position.

The difference column indicates how `hope` scored the sequencing of the homopolymer in each read. negative numbers indicate a deletion of the stated number of bases, while positive numbers indicate insertions. 0 indicates that the homopolymer was correctly sequenced.

In addition to numerical scores, the score may also be reported as either "skip" or "?". "Skip" indicates that no flanking sequence was avaialble on one side of the homopolymer so no score could be assigned. "?" indicates that something more complex than a simple homopolymer error was found. In the below example output, for example, an insertion of CAG is seen in the homopolymer in the read. As different bases were inserted, this is not considered by hope to be a simple extension of the homopolymer.


### Example output

The following shows a subset of output from analysis of real data.

Example output without `-c`

    homopolymer_length      homopolymer_base        difference      homo_start      read_ID
    5       T       0       26      6325e650-5c95-4988-97fc-a7daa945193c
    5       G       ?       990     6325e650-5c95-4988-97fc-a7daa945193c
    5       T       -2      1514    6325e650-5c95-4988-97fc-a7daa945193c
    5       A       -1      1684    6325e650-5c95-4988-97fc-a7daa945193c
   
Example output with `-c`

    homopolymer_length      homopolymer_base        difference      read_context    assembly_context        homo_start     read_ID
    5       T       0       ACGAAACAAACAACGTGAAACGTCAATTTTTTATTTTAGATGCTGA-ACAAGCTAAC---ATT     ACGAAACAAACAACGTGAAACGTCAA-TTTTTATTTTAGATGCT-AGACAAACTAACTTTATT 26      6325e650-5c95-4988-97fc-a7daa945193c
    5       G       ?       GCCGCAAGGCTGAAACTCAAAGG----GCCGGCAGGGGCCCG-----GCGGGTGGAGCATGTGGTTTAA       GCCGCAAGGCTGAAACTCAAAGGAATTGACG---GGGGCCCGCACAAGC-GGTGGAGCATGTGGTTTAA   990     6325e650-5c95-4988-97fc-a7daa945193c
    5       T       -2      ---TT-TAACACCCGAAGTCGGTGGGGTAACC--TTTGGAACCAACCACCT-AGGTGGGACAGATGA         AGTTTGTAACACCCGAAGTCGGTGGGGTAACCTTTTTGGAGCCAGCCGCCTAAGGTGGGACAGATGA     1514    6325e650-5c95-4988-97fc-a7daa945193c
    5       A       -1      TTCG-TTTCGTTTAGTTTTGAGAGTTCAAT-AAAAGTATTGACTCTTAAATGAGGATATGATATA           TTCGTTTTCGTTTAGTTTTGAGAGTTCAATAAAAAGTATTGACTCTTAAATGAGGATATGATATA       1684    6325e650-5c95-4988-97fc-a7daa945193c




    
