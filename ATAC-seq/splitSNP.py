#! /usr/bin/env python
def main():
    import argparse
    import pysam
    import os
    import tempfile
    from collections import Counter

    def can_create_file(folder_path):
        try:
            tempfile.TemporaryFile(dir=folder_path)
            return True
        except OSError:
            return False

    parser = argparse.ArgumentParser(description="splitSNP.py - "
        "Extracts reads from a BAM file that cover a specified SNP and "
        "writes the reference and alternate allele containing reads to separate BAM files.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("SNP", help="SNP position, reference and alternate "
        "allele of interest in chr:position:ref:alt format, eg chr21:11106932:A:G. "
        "For deletion analysis the ref should be 'D' and alt be the size of the deletion "
        "in basepairs, eg chr11:67351213:D:64. Coordinates are 1-based.")
    parser.add_argument("--pair_distance", type=int, default=500, help="The distance in "
        "basepairs to search up and downstream from the specified SNP/deletion for "
        "the pair of overlapping reads (default is 500)")
    parser.add_argument("--max_depth", type=int, default=1000000, help="Maximum number "
        "of reads to process at the specified SNP position (default is 1000000)")
    args = parser.parse_args()

    if args.max_depth < 8000:
        print("Specified max_depth is too low - changing to 8000")
        args.max_depth = 8000

    # Check input file exists, and thet output folder is writeable
    if not os.path.isfile(args.input_bam):
        print("Input BAM file %s does not exist!" % args.input_bam)
        return

    if not can_create_file(os.path.dirname(args.output_prefix)):
        print("Output path %s is not writable!" % os.path.dirname(args.output_prefix))
        return

    # Check alleles are valid
    nucleotides = ['A', 'C', 'G', 'T']
    try:
        chrom, pos, ref, alt = args.SNP.split(":")
    except ValueError:
        print("SNP specified '%s' not in chr:position:ref:alt format" % args.SNP)
        return

    ref = ref.upper()
    alt = alt.upper()
    if not ref=="D":
        if ref not in nucleotides:
            print("Reference allele %s is not A, C, G or T" % ref)
            return
        if alt not in nucleotides:
            print("Alternate allele %s is not A, C, G or T" % alt)
            return
    else: # deletion analysis
        try:
            alt = int(alt)
        except ValueError:
            print("Deletion length '%s' is not valid" % alt)
            return

    # Check validity of chrom
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    chrom_no = next((i for i in range(len(samfile.references)) if samfile.references[i]==chrom), None)
    if chrom_no is None:
        print("Chromosome '%s' not in BAM '%s'" % (chrom, args.input_bam))
        return
    
    # Check validity of pos
    try:
        pos = int(pos)
    except ValueError:
        print("Position '%s' is not valid" % pos)
        return
    if pos >= samfile.lengths[chrom_no]:
        print("Position '%s' is out of bounds of chromosome '%s'" % (pos, chrom))
        return
    
    # PASS 1 - find readnames of all reads with ref/alt alleles
    alleles = list()
    ref_readnames = set()
    alt_readnames = set()
    if ref!="D": # SNP analysis
        for pileup in samfile.pileup(chrom, pos-1, pos, max_depth=args.max_depth):
            if pileup.reference_pos == pos-1: # filter for position of interest
                print("Processing %s reads covering SNP position %s:%s in %s" % (
                    len(pileup.pileups), chrom, pos, args.input_bam))
                for read in pileup.pileups:
                    SNP_base = read.alignment.query_sequence[read.query_position]
                    alleles.append(SNP_base)
                    if SNP_base == ref:
                        ref_readnames.add(read.alignment.query_name)
                    elif SNP_base == alt:
                        alt_readnames.add(read.alignment.query_name)
    else: # Deletion analysis
        delbases = range(pos-1, pos+alt-1)
        for pileup in samfile.pileup(chrom, pos-1, pos+alt-1, max_depth=1000000):
            if pileup.reference_pos in delbases: # filter for positions of interest
                for read in pileup.pileups:
                    # Already checked these reads
                    if read.alignment.query_name in ref_readnames or read.alignment.query_name in alt_readnames:
                        continue
                    if read.is_del==0:
                        # Make sure there is evidence for sequencing across the deletion
                        if any(i in delbases for i in read.alignment.get_reference_positions()):
                            ref_readnames.add(read.alignment.query_name)
                    else:
                        # Make sure the whole deletion is present
                        if not any(i in delbases for i in read.alignment.get_reference_positions()):
                            alt_readnames.add(read.alignment.query_name)

    # Remove reads in both
    ref_and_alt = ref_readnames.intersection(alt_readnames)
    if len(ref_and_alt) > 0:
        print("%s reads discarded for being ambiguous" % len(ref_and_alt))
        ref_readnames = ref_readnames.difference(ref_and_alt)
        alt_readnames = alt_readnames.difference(ref_and_alt)

    # PASS 2 - output reads matching above readnames to two new bamfiles
    if ref!="D":
        ref_filename = args.output_prefix+".ref."+ref+".bam"
        alt_filename = args.output_prefix+".alt."+alt+".bam"
        minpos = pos-args.pair_distance
        maxpos = pos+args.pair_distance
    else:
        ref_filename = args.output_prefix+".ref.bam"
        alt_filename = args.output_prefix+".del.%sbp.bam" % alt
        minpos = pos-args.pair_distance
        maxpos = pos+alt+args.pair_distance

    # Handle SNPs near the beginning of a chromosome
    if minpos < 0:
        minpos = 0

    ref_bam = pysam.AlignmentFile(ref_filename, "wb", template=samfile)
    alt_bam = pysam.AlignmentFile(alt_filename, "wb", template=samfile)
    ref_count = 0
    alt_count = 0
    for read in samfile.fetch(chrom, minpos, maxpos): # extra 1bp buffer
        if read.query_name in ref_readnames:
            ref_bam.write(read)
            ref_count += 1
        elif read.query_name in alt_readnames:
            alt_bam.write(read)
            alt_count += 1

    samfile.close()
    ref_bam.close()
    alt_bam.close()
    pysam.index(ref_filename)
    pysam.index(alt_filename)

    # Print a summary
    if ref!="D":
        allele_count = Counter(alleles)
        print("%s read segments with the reference allele '%s' written to %s" % (
            ref_count, ref, ref_filename))
        print("%s read segments with the alternate allele '%s' written to %s" % (
            alt_count, alt, alt_filename))
        for x in allele_count:
            if x != ref and x != alt:
                print("Discarded %s '%s' alleles" % (allele_count[x], x))
    else:
        print("%s read segments with the reference sequence written to %s" % (
            ref_count, ref_filename))
        print("%s read segments with the %sbp deletion written to %s" % (
            alt_count, alt, alt_filename))

if __name__ == '__main__':
    main()
