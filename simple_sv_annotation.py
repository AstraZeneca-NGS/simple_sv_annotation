#!/usr/bin/env python

"""
Name:    simple_sv_annotation.py

Purpose: Simplify SV annotations from snpEff to highlight exon impact.

Input:   vcf file with SVs annotated with snpEff

Output:  vcf file with SV annotation field with simplified annotation

Usage:   simple_sv_annotation.py [OPTIONS]

Authors:  David Jenkins (david.jenkins1@astrazeneca.com/dfj@bu.edu), Miika Ahdesmaki (miika.ahdesmaki @ astrazeneca.com / live.fi)

Notes:
    - Simplified with SnpEff 4.3 updates that add better annotation for fusion events resulting from all BND/INV/DUP/DEL
    - No simplified annotation for within-exon deletions (frameshifts, in-frame dels) as not classified as SVs
    - Annotates fusions that are gene1-gene2 or gene1-intergenic. Soft filters gene-same gene fusions (likely noise)
    - Soft filters intergenic events that are not involved with genes directly

See the provided README.md file for more information
"""

from __future__ import print_function
import argparse, os, vcf, sys
try:
    import pysam
except ImportError:
    pysam = None

def main(vcf_in, outfile, remove_ann, exon_nums, args):
    """Adds additional header information, opens the infile for reading, opens
    the outfile for writing, and iterates through the records to identify
    records that are candidates for adding the simple annotation
    """

    vcf_reader = vcf.Reader(filename=vcf_in)
    # add a new info field into the header of the output file
    vcf_reader.infos['SIMPLE_ANN'] = vcf.parser._Info(id="SIMPLE_ANN", num=".", type="String", desc="Simplified human readable structural variant annotation: 'SVTYPE | ANNOTATION | GENE(s) | TRANSCRIPT | DETAIL'", source=None, version=None)
    # add filters
    vcf_reader.filters["REJECT"] = vcf.parser._Filter(id="REJECT", desc="Rejected due various criteria (e.g. gene-same gene fusion, intergenic deletion)")
    vcf_writer = vcf.Writer(open(outfile, 'w'), vcf_reader) if outfile !="-" else vcf.Writer(sys.stdout, vcf_reader)
    for record in vcf_reader:
        if record.FILTER is None:
            record.FILTER = []

        if 'SVTYPE' in record.INFO and 'ANN' in record.INFO:
            #any(["gene_fusion" in x for x in record.INFO['ANN']])
            vcf_writer.write_record(simplify_ann(record, remove_ann, exon_nums))
        else: 
            record.FILTER.append("REJECT")
            vcf_writer.write_record(record)
    vcf_writer.close()

def switch_to_ann(annotation):
    """Switch SIMPLE_ANN to a more ANN like form"""
    ann_a = annotation.split("|")
    return '%s|%s||%s|||%s||%s|||||||' % (ann_a[0],ann_a[1],ann_a[2],ann_a[3],ann_a[4])

def replace_ann_field(record):
    """Remove the ANN field from a record and replace it with the simplifed 
    record"""
    del record.INFO['ANN']
    record.INFO['ANN'] = record.INFO['SIMPLE_ANN']
    record.INFO['ANN'] = [switch_to_ann(i) for i in record.INFO['ANN']]
    del record.INFO['SIMPLE_ANN']

def simplify_ann(record, remove_ann, exon_nums):
    """Find any annotations that can be simplified and call the method
    to annotate it.
    """
    # marching order is: 'exon_loss_variant', fusions, others (reject)

    # to-do: CNV and INS?
    
    exon_losses = {}
    annotated = False
    is_intergenic = True # is intergenic or otherwise likely rubbish?
    for i in record.INFO['ANN']:
        ann_a = i.split('|')
        if "exon_loss_variant" in [ann_a[1]:
            is_intergenic = False
            try:
                exon_losses[ann_a[6]].append(i)
            except KeyError:
                exon_losses[ann_a[6]] = [i]
        elif "gene_fusion" in [ann_a[1]: 
            # This could be 'gene_fusion', 'bidirectional_gene_fusion' but not 'feature_fusion'
            # 'gene_fusion' could lead to a coding fusion whereas 
            # 'bidirectional_gene_fusion' is likely non-coding (opposing frames, _if_ inference correct)
            annotate_other_var(record, ann_a)
            annotated = True
            is_intergenic = False
    if len(exon_losses) > 0:
        annotate_exon_loss(record, exon_losses, exon_nums)
        annotated = True
    if remove_ann and annotated:
        replace_ann_field(record)
    # REJECT purely intergenic events and other nuisance variants
    if is_intergenic:
        record.FILTER.append("REJECT")
    return record

def uniq_list(inlist):
    """Remove unique elements from a list"""
    inset = set(inlist)
    return list(inset)

def find_deleted_exons(annotations):
    """Take the annotations for a particular transcript
    and parse them to find the numbers of the exons that have
    been deleted
    """
    exons = []
    gene = ''
    for i in annotations:
        ann_a = i.split('|')
        if gene == '':
            gene = ann_a[3]
        try:
            exons.append(int(ann_a[8].split('/')[0]))
        except ValueError:
            pass
    return exons, gene

def find_alt_deleted_exons(record, exon_nums, annotations):
    """In the case where the user has provided a file of alternate
    transcript numbers, use those numbers to find the exons that have been
    deleted
    """
    exons = []
    gene = ''
    for i in annotations:
        ann_a = i.split('|')
        if gene == '':
            gene = ann_a[3]
        else:
            break
    for i in exon_nums:
        if int(i[0]) > record.POS and int(i[1]) <= record.INFO['END']:
            exons.append(int(i[2]))
    return exons, gene

def annotate_exon_loss(record, exon_losses, exon_nums):
    """Create the exon loss simple annotation from the exon dict created
    in simplify_ann
    
    For each transcript with exon losses, find the numbers for each exon
    and create the annotation
    Example: DEL|EXON_DEL|BLM|NM_001287247.1|Exon2-12del
    """
    for transcript in exon_losses:
        #Remove version number if it exists
        if transcript.split('.')[0] in exon_nums:
            #use alternate exon numbers for transcript
            exons, gene = find_alt_deleted_exons(record, exon_nums[transcript.split('.')[0]], exon_losses[transcript])
        else:
            #use snpEff numbers for transcript
            exons, gene = find_deleted_exons(exon_losses[transcript]) 
        exons = uniq_list(exons)
        if len(exons) == 0:
            return None
        if max(exons)-min(exons)+1 == len(exons):
            if len(exons) == 1:
                deleted_exons = "Exon"+str(exons[0])+"del"
            else:
                deleted_exons = "Exon"+str(min(exons))+"-"+str(max(exons))+"del"
        else:
            deleted_exons = "Exon"+str(min(exons))+"-"+str(max(exons))+"del"
        try:
            record.INFO['SIMPLE_ANN'].append("DEL|EXON_DEL|%s|%s|%s" % (gene,transcript,deleted_exons))
        except KeyError:
            record.INFO['SIMPLE_ANN'] = ["DEL|EXON_DEL|%s|%s|%s" % (gene,transcript,deleted_exons)]

def annotate_other_var(record, ann_a):
    """Create a simplified version of the annotation field for non-whole exon loss events

    Regardless of sv type, the simple annotation for an intronic variant
    looks like: SIMPLE_ANN=INV|INTRONIC|ERBB4|NM_005235.2|
    """
    simple_ann = "%s|%s|%s|%s|" % (record.INFO['SVTYPE'], ann_a[1].upper(), ann_a[3], ann_a[6])
    try:
        if simple_ann not in record.INFO['SIMPLE_ANN']: # avoid duplicate entries that bloat the output
            record.INFO['SIMPLE_ANN'].append(simple_ann)
    except KeyError:
        record.INFO['SIMPLE_ANN'] = [simple_ann]

# UNUSED FUNCTION FOR NOW
def annotate_intergenic_var(record, ann_a):
    """Create a simplified version of the annotation field for an intergenic var

    Regardless of sv type, the simple annotation for an intergenic variant
    looks like: SIMPLE_ANN=INV|INTERGENIC|LETM2-FGFR1||
    """
    simple_ann = "%s|INTERGENIC|%s||" % (record.INFO['SVTYPE'],ann_a[4])
    try:
        if simple_ann not in record.INFO['SIMPLE_ANN']: # avoid duplicate entries that bloat the output
            record.INFO['SIMPLE_ANN'].append(simple_ann)
    except KeyError:
        record.INFO['SIMPLE_ANN'] = [simple_ann]

def create_exon_numDict(infile):
    """Create a dictionary of exon numbers based on an input file of alternate
    chromosome numberings.

    example: BRCA1 has no exon 4 (http://www.medscape.com/viewarticle/567639_2),
    making exon numbers from snpEff incorrect.
    """
    exons = {}
    f = open(infile)
    for line in f:
        la = line.rstrip("\n").split("\t");
        name = la[3].split("|")
        transcript = name[0]
        if transcript in exons:
            exons[transcript].append((la[1],la[2],name[1]))
        else:
            exons[transcript] = [(la[1],la[2],name[1])]
    return exons

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Simplify SV annotations from snpEff to highlight exon impact. Requires the pyvcf module.")
    parser.add_argument('vcf', help='VCF file with snpEff annotations')
    parser.add_argument('--output', '-o', help='Output file name (must not exist). Does not support bgzipped output. Use "-" for stdout. [<invcf>.simpleann.vcf]', required=False)
    parser.add_argument('--exonNums', '-e', help='List of custom exon numbers. A transcript listed in this file will be annotated with the numbers found in this file, not the numbers found in the snpEff result')
    parser.add_argument('-r', help="Instead of creating a SIMPLE_ANN field, replace the ANN field with a simplified version that will retain the same inforamtion in the same fields", action='store_true', default=False)

    parser_excl = parser.add_mutually_exclusive_group(required=False)
    args = parser.parse_args() 
    if args.output:
        outfile = args.output
    else:
        outfile = args.vcf.replace("vcf.gz","vcf").replace(".vcf", ".simpleann.vcf")
    if os.path.exists(outfile):
        raise IOError("Output file %s exists" % outfile)
    exonNumDict = {} 
    if args.exonNums:
        exonNumDict = create_exon_numDict(args.exonNums) 
    main(args.vcf, outfile, args.r, exonNumDict, args)
