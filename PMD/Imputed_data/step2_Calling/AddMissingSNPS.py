import gzip
import sys

def add_missing(reference_file, atlas_file, output_file):
    #Reads and store the header from the Atlas file
    atlas_header = []
    with gzip.open(atlas_file, 'rt') as atl:
        for line in atl:
            if line.startswith('#'):
                atlas_header.append(line.strip())
            else:
                break

    #Opens reference file and stores non-header lines
    with gzip.open(reference_file, 'rt') as ref:
        reference_lines = [line.strip() for line in ref if not line.startswith('#')]
    
	#Opens Atlas file and stores non-header lines
    with gzip.open(atlas_file, 'rt') as atl:
        atlas_lines = [line.strip() for line in atl if not line.startswith('#')]

    #Creates a dictionary from Atlas lines (SNPs) with the first two elements (Chromosome, Position) as key
    atlas_dict = {tuple(line.split()[:2]): line for line in atlas_lines}
    merged_lines = []

    #Iterates through all the SNPs of the reference
    for ref_line in reference_lines:
        chrom, pos, id, ref, alt = ref_line.split()[:5]
        key = (chrom, pos)
        #If the SNP is missing from the Atlas file -> creates a "placeholder" SNP using the Chromosome, Position, REF and ALT allele information from the reference
        if key not in atlas_dict:
            placeholder_line = f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t.\t.\t.\tGT\t."
            merged_lines.append(placeholder_line)
        else:
        #If the SNP is present in the Atlas file -> saves this SNP as is from the Atlas file
            merged_lines.append(atlas_dict[key])

    #Writes the header and lines to the output file
    with gzip.open(output_file, 'wt') as out:
        out.write("\n".join(atlas_header) + "\n")
        out.write("\n".join(merged_lines))

#Files: 1:Reference file, 2:Atlas file, 3:Output name
add_missing(merge_vcf(sys.argv[1], sys.argv[2], sys.argv[3]) 
