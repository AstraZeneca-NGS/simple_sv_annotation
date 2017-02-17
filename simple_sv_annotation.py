#!/usr/bin/env python

"""
Name:    simple_sv_annotation.py

Purpose: Simplify SV annotations from snpEff to highlight exon impact.

Input:   vcf file with SVs annotated with snpEff 4.3 or higher, text file with known fusion pairs, text file with genes of interest

Output:  vcf file with SV annotation field with simplified annotation

Usage:   simple_sv_annotation.py [OPTIONS]

Authors:  David Jenkins (david.jenkins1@astrazeneca.com/dfj@bu.edu), Miika Ahdesmaki (miika.ahdesmaki @ astrazeneca.com / live.fi)

Notes:
    - Simplified with SnpEff 4.3 updates that add better annotation for fusion events resulting from all BND/INV/DUP/DEL
    - No simplified annotation for within-exon deletions (frameshifts, in-frame dels) as not classified as SVs
    - Soft filters intergenic events that are not involved with genes directly

Current scheme with priority 1(high)-3(low)
- exon loss 
 . on prioritisation gene list (2)
 . not on gene list (3)
- gene_fusion
. paired (hits two genes)
.. on list of known pairs (1)
.. not on list of known pairs 
... one or two genes on gene list (2)
... neither gene on gene list (3)
. unpaired (hits one gene)
.. on gene list (2)
.. not on gene list (3)
- upstream or downstream of gene list genes (3)
- other variant (REJECT)
- missing ANN or SVTYPE: REJECT

See the provided README.md file for more information
"""

from __future__ import print_function
import argparse, os, vcf, sys
try:
    import pysam
except ImportError:
    pysam = None

def main(vcf_in, outfile, exon_nums, args):
    """Adds additional header information, opens the infile for reading, opens
    the outfile for writing, and iterates through the records to identify
    records that are candidates for adding the simple annotation
    """
    vcf_reader = vcf.Reader(filename=vcf_in)
    #
    # add a new info field into the header of the output file
    #
    vcf_reader.infos['SIMPLE_ANN'] = vcf.parser._Info(id="SIMPLE_ANN", num=".", type="String", desc="Simplified human readable structural variant annotation: 'SVTYPE | ANNOTATION | GENE(s) | TRANSCRIPT | DETAIL (exon losses, KNOWN_FUSION, ON_PRIORITY_LIST, NOT_PRIORITISED) | PRIORITY (1-3) '", source=None, version=None)
    #
    # Add an info field for highest priority (given multiple annotations per entry)
    vcf_reader.infos['SV_HIGHEST_TIER'] = vcf.parser._Info(id="SV_HIGHEST_TIER", num=1, type="Integer", desc="Highest priority tier for the effects of a variant entry", source=None, version=None)
    #

    # Add filters
    #
    vcf_reader.filters["REJECT"] = vcf.parser._Filter(id="REJECT", desc="Rejected due various criteria (missing ANN/BND, purely intergenic, small events)")
    vcf_writer = vcf.Writer(open(outfile, 'w'), vcf_reader) if outfile !="-" else vcf.Writer(sys.stdout, vcf_reader)
    #
    # Read in gene lists
    #
    known_fusions, prioritised_genes = read_gene_lists(args.known_fusion_pairs, args.gene_list)
    #
    #
    #
    for record in vcf_reader:
        if record.FILTER is None:
            record.FILTER = []

        if 'SVTYPE' in record.INFO and 'ANN' in record.INFO:
            #any(["gene_fusion" in x for x in record.INFO['ANN']])
            vcf_writer.write_record(simplify_ann(record, exon_nums, known_fusions, prioritised_genes))
        else: 
            record.FILTER.append("REJECT")
            vcf_writer.write_record(record)
    vcf_writer.close()

def read_gene_lists(known_fusion_pairs, gene_list):
    kfp = []
    gl = []
    if known_fusion_pairs and os.path.isfile(known_fusion_pairs):
        with open(known_fusion_pairs, 'r') as myfhandle:
            for line in myfhandle:
                genes = line.strip().split(",")
                if len(genes) == 2 and len(genes[0]) > 0 and len(genes[1]) > 0:
                    kfp.append([genes[0].strip(), genes[1].strip()])
    if gene_list and os.path.isfile(gene_list):
        with open(gene_list, 'r') as myghandle:
            for line in myghandle:
                gene = line.strip()
                if len(gene) > 0:
                    gl.append(gene)
    return kfp, gl

def simplify_ann(record, exon_nums, known_fusions, prioritised_genes):
    """Find any annotations that can be simplified and call the method
    to annotate it.
    """
    # marching order is: 'exon_loss_variant', fusions, others (reject)

    # to-do: CNV and INS?
    
    exon_losses = {}
    annotated = False
    is_intergenic = True # is intergenic or otherwise likely rubbish?
    record.INFO['SV_HIGHEST_TIER'] = 3
    for i in record.INFO['ANN']:
        ann_a = i.split('|')
        if "exon_loss_variant" in ann_a[1]:
            is_intergenic = False
            try:
                exon_losses[ann_a[6]].append(i)
            except KeyError:
                exon_losses[ann_a[6]] = [i]
        elif "gene_fusion" in ann_a[1] or "downstream" in ann_a[1] or "upstream" in ann_a[1]: 
            # This could be 'gene_fusion', 'bidirectional_gene_fusion' but not 'feature_fusion'
            # 'gene_fusion' could lead to a coding fusion whereas 
            # 'bidirectional_gene_fusion' is likely non-coding (opposing frames, _if_ inference correct)
            annotate_other_var(record, ann_a, known_fusions, prioritised_genes)
            annotated = True
            is_intergenic = False
        elif "downstream" in ann_a[1] or "upstream" in ann_a[1]:
            # get SVs affecting up/downstream regions of prioritised genes
            if len(ann_a[3]) > 0 and ann_a[3] in prioritised_genes:
                var_priority = "3"
                var_detail = "ON_PRIORITY_LIST"
                simple_ann = "%s|%s|%s|%s|%s|%s" % (record.INFO['SVTYPE'], ann_a[1].upper(), ann_a[3], ann_a[6], var_detail, var_priority)
                try:
                    if simple_ann not in record.INFO['SIMPLE_ANN']: # avoid duplicate entries that bloat the output
                        record.INFO['SIMPLE_ANN'].append(simple_ann)
                except KeyError:
                    record.INFO['SIMPLE_ANN'] = [simple_ann]
                annotated = True
                is_intergenic = False
    if len(exon_losses) > 0:
        annotate_exon_loss(record, exon_losses, exon_nums, prioritised_genes)
        annotated = True
    # REJECT purely intergenic events and other nuisance variants
    if is_intergenic:
        record.FILTER.append("REJECT")
        del record.INFO["SV_HIGHEST_TIER"]
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

def annotate_exon_loss(record, exon_losses, exon_nums, prioritised_genes):
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
        var_priority = "2" if gene in prioritised_genes else "3"
        if record.INFO['SV_HIGHEST_TIER'] > int(var_priority):
            record.INFO['SV_HIGHEST_TIER'] = int(var_priority)
        try:
            record.INFO['SIMPLE_ANN'].append("DEL|EXON_DEL|%s|%s|%s|%s" % (gene,transcript,deleted_exons,var_priority))
        except KeyError:
            record.INFO['SIMPLE_ANN'] = ["DEL|EXON_DEL|%s|%s|%s|%s" % (gene,transcript,deleted_exons,var_priority)]

def annotate_other_var(record, ann_a, known_fusions, prioritised_genes):
    """Create a simplified version of the annotation field for non-whole exon loss events
    that are likely to lead to fusions
    Regardless of sv type, the simple annotation for an intronic variant
    looks like: SIMPLE_ANN=INV|GENE_FUSION|ALK&EML4|NM_...&NM_...|KNOWN_FUSION|1
    """
    var_priority = "3"
    var_detail = "NOT_PRIORITISED"
    genes = ann_a[3].strip().split("&")

    if len(genes) > 1 and ([genes[0],genes[1]] in known_fusions or [genes[1],genes[0]] in known_fusions):
        var_priority = "1"
        var_detail = "KNOWN_FUSION"
    elif (len(genes[0]) > 0 and genes[0] in prioritised_genes) or (len(genes) > 1 and len(genes[1])>0 and genes[1] in prioritised_genes):
        var_detail = "ON_PRIORITY_LIST"
        var_priority = "2"
    if record.INFO['SV_HIGHEST_TIER'] > int(var_priority):
        record.INFO['SV_HIGHEST_TIER'] = int(var_priority)
    simple_ann = "%s|%s|%s|%s|%s|%s" % (record.INFO['SVTYPE'], ann_a[1].upper(), ann_a[3], ann_a[6], var_detail, var_priority)
    try:
        if simple_ann not in record.INFO['SIMPLE_ANN']: # avoid duplicate entries that bloat the output
            record.INFO['SIMPLE_ANN'].append(simple_ann)
    except KeyError:
        record.INFO['SIMPLE_ANN'] = [simple_ann]

# UNUSED FUNCTION FOR NOW
#def annotate_intergenic_var(record, ann_a):
#    """Create a simplified version of the annotation field for an intergenic var
#
#    Regardless of sv type, the simple annotation for an intergenic variant
#    looks like: SIMPLE_ANN=INV|INTERGENIC|LETM2-FGFR1||
#    """
#    simple_ann = "%s|INTERGENIC|%s||" % (record.INFO['SVTYPE'],ann_a[4])
#    try:
#        if simple_ann not in record.INFO['SIMPLE_ANN']: # avoid duplicate entries that bloat the output
#            record.INFO['SIMPLE_ANN'].append(simple_ann)
#    except KeyError:
#        record.INFO['SIMPLE_ANN'] = [simple_ann]

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
    parser.add_argument('--gene_list', '-g', help='File with names of genes (one per line) for prioritisation', required=False, default=None)
    parser.add_argument('--known_fusion_pairs', '-k', help='File with known fusion gene pairs, one pair per line delimited by comma', required=False, default=None)
    parser.add_argument('--output', '-o', help='Output file name (must not exist). Does not support bgzipped output. Use "-" for stdout. [<invcf>.simpleann.vcf]', required=False)
    parser.add_argument('--exonNums', '-e', help='List of custom exon numbers. A transcript listed in this file will be annotated with the numbers found in this file, not the numbers found in the snpEff result')
    #parser_excl = parser.add_mutually_exclusive_group(required=False)
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
    main(args.vcf, outfile, exonNumDict, args)
