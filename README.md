## basecounts

`basecounts` is a simple tool to extract nucleotide counts from BAM files for the user specific loci. This tool can also be used to generate VAF table for known somatic varaints without the proper variant calling. 

### Installation

Clone the repository and run `make`. Requires [htslib](https://github.com/samtools/htslib). 

```
git clone https://github.com/PoisonAlien/basecounts
cd basecounts
make
```

```
$ basecounts
-------------------------------------------------------------------------
basecounts: A tool to extract nucleotide counts
            of genomic loci from the BAM file. Version: 0.1.0
USAGE:
    basecounts [OPTIONS] <loci> <bam>
    e.g; basecounts target_loci.tsv Tumor.bam
OPTIONS:
    -f  Indexed fasta file. If provided, extracts and adds reference base to the ouput tsv
    -q  Min. mapping quality threshold. Default 10 [Read filter]
    -F  Exclude reads with FLAGS >= INT. Default 3852 [Read filter]
    -o  Output file basename. Default parses from basename of BAM file
ARGS:
    <loci>   loci file (chr\tpos)
    <bam>    Indexed BAM file
OUPUT:
    <output>.tsv    TSV file with nucletide counts for all variants
-------------------------------------------------------------------------

```

### Usage

`basecounts` takes a 2 column tsv file with loci and the sorted/indexed BAM file as positional arguments.


First two columns of loci file should be chromsome and position.
```
$ head loci.tsv
1	6304615
1	10407830
1	46269731
1	47133825
1	50884996
1	84662414
1	152328587
1	216019313
1	226252155
```

```
#basic usage
$ basecounts loci.tsv tumor.bam 

#providing fasta file will add the reference base at the loci
$ basecounts -f GRCh37.fa loci.tsv tumor.bam 

#Output is a TSV file with readcounts suporting all 4 nucleotide bases and INDELS
$ head tumor.tsv 
loci    fa_ref  A       T       G       C       Ins     Del
1:6304615       A       70      0       0       19      0       0
1:10407830      C       0       0       39      50      0       0
1:46269731      A       50      3       1       21      0       0
1:47133825      T       3       37      14      4       0       0
1:50884996      A       51      0       0       15      0       0
1:84662414      G       0       2       81      0       0       0
1:152328587     A       85      3       0       0       0       0
1:216019313     T       0       53      17      0       0       0
1:226252155     G       0       41      52      0       0       0
```
