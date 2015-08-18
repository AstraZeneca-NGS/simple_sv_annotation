simple_sv_annotation.py
=======================

A tool for simplifying snpEff annotations

## Table of Contents

1. [Requirements](#requirements)
2. [Usage](#usage)
3. [Alternate Exon Numbers](#alternate-exon-numbers)
4. [Supported SV Types](#supported-sv-types)
5. [Supported SV Callers](#supported-sv-callers)
6. [Example Output](#example-output)

## Requirements

1. python 2.7
2. [PyVcf](http://pyvcf.readthedocs.org/en/latest/) python module
3. [VCF](https://vcftools.github.io/specs.html) file annotated with [snpEff v4.1g+](http://snpeff.sourceforge.net/).

```simple_sv_annotation.py``` is designed around the new [ANN](http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf) annotation field rather than the previous EFF field.

## Usage

```
usage: simple_sv_annotation [options] vcf
```

**Required arguments**

```
vcf FILE - vcf file annotated with snpEff v4.1g+
```

**Optional arguments**

```
--output/-o   FILE - Output file name. Default: <invcf>.simpleann.vcf
--exonNums/-e FILE - List of custom exon numbers (see Alternate Exon Numbers)
-r                 - Replace the ANN field instead of adding SIMPLE ANN (see Example Output)
```

## Alternate Exon Numbers

Occasionally the exon numbering scheme provided by snpEff is incorrect. snpEff
numbers the exons in a transcript sequentially, but sometimes the accepted exon
numbering is not sequential. For example, BRCA1 transcript 1, NM_007294, [does
not have an exon 4](http://www.medscape.com/viewarticle/567639_2).

```simple_sv_annotation.py``` accepts a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
file in which a user can provide custom numbering for a particular transcript. If
a variant is annotated with a transcript listed in this file, the exon numbers
provided by snpEff are replaced with the exon numbers in the file. If a
transcript is not in the file, then the snpEff exon numbers are used. Follow the
format below, separating each field with a tab

```
chr17    41196311    41197819    NM_007294|24
chr17    41199659    41199720    NM_007294|23
chr17    41201137    41201211    NM_007294|22
chr17    41203079    41203134    NM_007294|21
chr17    41209068    41209152    NM_007294|20
chr17    41215349    41215390    NM_007294|19
chr17    41215890    41215968    NM_007294|18
chr17    41219624    41219712    NM_007294|17
chr17    41222944    41223255    NM_007294|16
chr17    41226347    41226538    NM_007294|15
chr17    41228504    41228631    NM_007294|14
chr17    41234420    41234592    NM_007294|13
chr17    41242960    41243049    NM_007294|12
chr17    41243451    41246877    NM_007294|11
chr17    41247862    41247939    NM_007294|10
chr17    41249260    41249306    NM_007294|9
chr17    41251791    41251897    NM_007294|8
chr17    41256138    41256278    NM_007294|7
chr17    41256884    41256973    NM_007294|6
chr17    41258472    41258550    NM_007294|5
chr17    41267742    41267796    NM_007294|3
chr17    41276033    41276132    NM_007294|2
chr17    41277287    41277500    NM_007294|1
```

In the fourth column, provide the transcript name followed by a "```|```"
and then the exon number. Note that the transcript version is not used.
You may have additional fields in the bed file, ```simple_sv_annotation.py```
will only consider the first four.

Note: currently this list of alternate exons is stored in memory because it is
expected to be relatively small. Very large lists of alternate exon numbering
may affect performance.

## Supported SV Types

```simple_sv_annotation.py``` will attempt to simplify interesting and easy
SV types to make the annotation result more interpretable. If you have an 
additional SV type that you want to be able to simplify, please email David
Jenkins, [AZ Email](mailto:david.jenkins1@astrazeneca.com) or [BU Email](mailto:dfj@bu.edu).

1. Intergenic SVs
2. Intronic SVs
3. Whole Exon Loss SVs
4. Gene Fusions Annotated as breakends

Examples of the simplified SV annotations are below.

## Supported SV Callers

```simple_sv_annotation.py``` has been tested on annotated vcf output files from
the following SV callers:

1. [lumpy-sv](https://github.com/arq5x/lumpy-sv)
2. [manta](https://github.com/Illumina/manta)
3. [delly](https://github.com/tobiasrausch/delly)

Additional SV callers will also work with ```simple_sv_annotation.py``` if VCF
specifications are followed and each SV is described with standard SV INFO fields:

1. SVTYPE
2. MATEID (for SVTYPE=BND)
3. END (for whole exon deletions)

## Example Output

There are two different outputs that ```simple_sv_annotation.py``` can produce:

#### 1. Add SIMPLE_ANN field

In the default mode, ```simple_sv_annotation.py``` will not alter the ANN field
provided by snpEff. Instead an additional field called SIMPLE_ANN will be added
to the SV call. A SIMPLE_ANN will only be added to variants that can be
simplified, other variants are not altered.

There are five fields in the SIMPLE_ANN tag separated by "```|```".

1. SV type (deletion, duplication, insertion, breakend)
2. Annotation (fusion, exon loss, intergenic, intronic)
3. Gene name
4. Transcript name
5. For exon loss variants, deleted exon numbers (Exon5del)

example:

```
before:

chr17  41258467  del_5  ATATACCTTTTGGTTATATCATTCTTACATAAAGGACACTGTGAAGGCCCTTTCTTCTGGTTGAGAAGTTTCAGCATGCAAAATCTATA  A  .  .  END=41258555;SVTYPE=DEL;SVLEN=-88;UPSTREAM_PAIR_COUNT=0;DOWNSTREAM_PAIR_COUNT=0;PAIR_COUNT=0;ANN=A|exon_loss_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&splice_region_variant&splice_region_variant&splice_region_variant&intron_variant&intron_variant|HIGH|BRCA1|BRCA1|transcript|NM_007294.3|Coding|4/23|c.135-5_212+5delTATAGATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGTATA||||||

after:

chr17  41258467  del_5  ATATACCTTTTGGTTATATCATTCTTACATAAAGGACACTGTGAAGGCCCTTTCTTCTGGTTGAGAAGTTTCAGCATGCAAAATCTATA  A  .  .  END=41258555;SVTYPE=DEL;SVLEN=-88;UPSTREAM_PAIR_COUNT=0;DOWNSTREAM_PAIR_COUNT=0;PAIR_COUNT=0;ANN=A|exon_loss_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&splice_region_variant&splice_region_variant&splice_region_variant&intron_variant&intron_variant|HIGH|BRCA1|BRCA1|transcript|NM_007294.3|Coding|4/23|c.135-5_212+5delTATAGATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGTATA||||||;SIMPLE_ANN=DEL|EXON_DEL|BRCA1|NM_007294.3|Exon5del
```

#### 2. Replace ANN field 

Optionally, you may choose to replace the 16 column ANN field with the simplified
annotation information. This will roughly preserve the type of information
expected in each field. The five fields of the SIMPLE_ANN annotation will be
placed in the 1st, 2nd, 4th, 7th, and 9th column of the ANN record

```
before:

chr17  41258467  del_5  ATATACCTTTTGGTTATATCATTCTTACATAAAGGACACTGTGAAGGCCCTTTCTTCTGGTTGAGAAGTTTCAGCATGCAAAATCTATA  A  .  .  END=41258555;SVTYPE=DEL;SVLEN=-88;UPSTREAM_PAIR_COUNT=0;DOWNSTREAM_PAIR_COUNT=0;PAIR_COUNT=0;ANN=A|exon_loss_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&splice_region_variant&splice_region_variant&splice_region_variant&intron_variant&intron_variant|HIGH|BRCA1|BRCA1|transcript|NM_007294.3|Coding|4/23|c.135-5_212+5delTATAGATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGTATA||||||

after:

chr17  41258467  del_5  ATATACCTTTTGGTTATATCATTCTTACATAAAGGACACTGTGAAGGCCCTTTCTTCTGGTTGAGAAGTTTCAGCATGCAAAATCTATA  A  .  .  END=41258555;SVTYPE=DEL;SVLEN=-88;UPSTREAM_PAIR_COUNT=0;DOWNSTREAM_PAIR_COUNT=0;PAIR_COUNT=0;ANN=DEL|EXON_DEL||BRCA1|||NM_007294.3||Exon5del|||||||
```
