## basecounts

`basecounts` is a simple tool to extract nucleotide counts from BAM files for user specific loci. This tool can also be used to generate VAF table for known somatic varaints without the proper variant calling. 

### Usage

```rust
>basecounts
basecounts 0.1.0
A tool to generate nucleotide counts and VAF table from user specific loci or gentotypes

USAGE:
    basecounts [FLAGS] [OPTIONS] <loci> <bam>...

FLAGS:
    -h, --help       Prints help information
    -r               Assumes input loci are in 4 column format. i.e, Chr\tStart\tRef\tAlt
    -a               Prints allelic fractions (vaf) instead of counts
    -V, --version    Prints version information
    -v               Show progress bar

OPTIONS:
    -f <fasta>              Indexed fasta file. If provided, extracts and adds reference base to the ouput
    -p <placeholder>        Optional column id from loci file which should be included in placeholder slot of output

ARGS:
    <loci>      A 2 column tsv file with chr\tpos or a 4 column tsv with Chr\tStart\tRef\tAlt
    <bam>...    One or more BAM files

```

## Get basecounts at specific loci

### Input
Input can be a 2 column tsv file with `chr\tpos` or 
```
10	100995296
10	102506042
10	104182701
10	104678772
10	134579307
```

a 4 column tsv with `Chr\tStart\tRef\tAlt`
```
10	100995296	A	C
10	102506042	C	T
10	104182701	T	A
10	104678772	C	T
10	134579307	T	C
10	16990523	G	T
```
### Output

Below command generates nucleotide frequency at the target sites in the format `A|T|G|C|Ins|Del` 
```bash
$basecounts vars.tsv Normal.bam Tumor.bam | head
chr	pos	refbase	spaceholder	Normal	Tumor
10	100995295	-	-	104|0|0|0|0|0	106|0|0|6|0|0
10	102506041	-	-	1|0|0|48|0|0	0|12|0|23|0|0
10	104182700	-	-	0|95|1|0|0|0	8|62|1|0|0|0
10	104678771	-	-	1|0|0|82|0|0	0|38|1|65|0|0
10	134579306	-	-	0|70|0|0|0|0	0|40|0|6|0|0
```

Print allelic frequncies instead of counts
```bash
$basecounts -a vars.tsv Normal.bam Tumor.bam | head
chr	pos	refbase	spaceholder	Normal	Tumor
10	100995295	-	-	1|0|0|0|0|0	0.946|0|0|0.054|0|0
10	102506041	-	-	0.02|0|0|0.98|0|0	0|0.343|0|0.657|0|0
10	104182700	-	-	0|0.99|0.01|0|0|0	0.113|0.873|0.014|0|0|0
10	104678771	-	-	0.012|0|0|0.988|0|0	0|0.365|0.01|0.625|0|0
10	134579306	-	-	0|1|0|0|0|0	0|0.87|0|0.13|0|0
```

In case of 4 column data where ref and alt allele are known..
```bash
$basecounts -r -a vars.tsv Normal.bam Tumor.bam | head
chr	pos	genotype	spaceholder	Normal	Tumor
10	100995295	A/C	-	0	0.054
10	102506041	C/T	-	0	0.343
10	104182700	T/A	-	0	0.113
10	104678771	C/T	-	0	0.365
10	134579306	T/C	-	0	0.13
```

This also allows to infer genotypes of known somatic variants without proper variant calling. For example, get allelic frequencies ot readcounts for all known `TP53` loci
```bash
$basecounts -r -a tp53_known_sites.tsv Normal.bam Tumor.bam | head
chr	pos	genotype	spaceholder	Normal	Tumor
17	7579357	C/A	-	0	0.014
17	7577534	C/A	-	0	0.011
17	7577533	C/A	-	0.024	0.011
17	7578175	C/A	-	0	0.01
17	7574020	C/A	-	0	0.007
17	7579363	C/-	-	0	0
```

`basecounts` can be helpful in case of multi-region sequening where one would need to generate VAF table for every somatic variant across all the sequenced samples.
```bash
$basecounts -r -a tp53_known_sites.tsv Germline.bam primary.bam metastasis.bam
chr	pos	genotype	spaceholder	Germline	primary	metastasis
17	7578402	C/A	-	0	0.022	0
17	7577120	G/A	-	0	0.015	0
17	7577123	C/T	-	0	0.008	0
17	7578289	C/T	-	0	0.005	0
17	7574016	C/A	-	0	0.004	0
17	7579357	C/A	-	0	0.003	0
17	7578211	G/A	-	0	0.003	0
17	7576851	C/T	-	0	0.003	0
17	7577533	C/A	-	0	0.002	0.018
17	7578474	G/A	-	0	0.002	0.002
```

### Installation

Install [rust](https://www.rust-lang.org/tools/install) and run the below command.

```
cargo install --git https://github.com/PoisonAlien/basecounts
```