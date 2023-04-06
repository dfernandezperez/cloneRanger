import gzip
import re
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#------------------------------------------------------------------------------------------------------------------------------------
# Declare input/output files
#------------------------------------------------------------------------------------------------------------------------------------
input_fb  = snakemake.input["fb"]
input_cb  = snakemake.input["cb"]
output_fb = snakemake.output["filt_fb"]
output_cb = snakemake.output["filt_cb"]
barcode_patterns = snakemake.params["barcode_dict"]


#------------------------------------------------------------------------------------------------------------------------------------
# Declare functions
#------------------------------------------------------------------------------------------------------------------------------------
def extract_barcodes(input_file, barcode_patterns):
    """
    This function will take a fastq file from larry barcode enrichment and a dictionary in which the
    keys are the barcode patterns and the corresponding larry library name/color  (i.e: {...TG...AG... : GFP}).

    It will parse the fastq file, look for reads that contain the barcode and extract the barcode sequence
    from the read. All the fastq records containing barcodes will be stored in a dictionary, as a list. The
    keys of the dictionary are the larry colors, and the values are a list of the fastq records containing each
    specific barcode color.
    """

    patterns      = [re.compile(r'{}'.format(barcode)) for barcode in barcode_patterns.keys()]
    patterns_dict = {key: value for key, value in zip(patterns, barcode_patterns.values())}
    barcode_reads = {key: [] for key in barcode_patterns.values()}

    with gzip.open(input_file, 'rt') as input_handle:
        # Iterate over each record in the fastq file
        for title, seq, qual in FastqGeneralIterator(input_handle):
            # Look for barcode patterns
            for pattern in patterns:
                match = pattern.search(seq)
                if match:
                    # Update record to contain just matched sequence 
                    seq  = seq[match.start() : match.end()]
                    qual = qual[match.start() : match.end()]
                    # Save updated records to dictionary. Every key is a different barcode color, the content are all the records
                    # corresponding to that color.
                    barcode_reads[patterns_dict[pattern]] += [title, seq, qual],
                    break

    return(barcode_reads)


def extracted_bc_to_fq(output_fastqs, barcode_reads):
    """
    Take the dictionary of larry colors and fastq records containing barcodes and save it to multiple fastq
    files, one for each larry color.
    """
    # Order of output files is the same as the keys of the dictionary since both of them are taken from the same config variable
    for i, key in enumerate(barcode_reads.keys()):
        with gzip.open(output_fastqs[i], 'wt') as output_handle:
            for title, seq, qual in barcode_reads[key]:
                _ = output_handle.write(f"@{title}\n{seq}\n+\n{qual}\n")


def subset_cb_fastq(input_fastq, output_fastq, barcode_reads):
    """
    """
    read_primer  = re.compile(r"\s.*$")
    filtered_ids = {read_primer.sub('', fastq_record[0]) for larry_color in barcode_reads.values() for fastq_record in larry_color}

    with gzip.open(input_fastq, 'rt') as input_handle, gzip.open(output_fastq, 'wt') as output_handle:
        for title, seq, qual in FastqGeneralIterator(input_handle):
            title_no_read_id = read_primer.sub('', title)
            if title_no_read_id in filtered_ids:
                _ = output_handle.write(f"@{title}\n{seq}\n+\n{qual}\n")
                


#------------------------------------------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------------------------------------------
# Filter reads containing barcodes from feature barcoding read and write to new fastq file
barcode_reads = extract_barcodes(input_fb, barcode_patterns) 
extracted_bc_to_fq(output_fb, barcode_reads)

# Subset the other fastq (the one containing cellular barcodes) based on the read ids filtered before
subset_cb_fastq(input_cb, output_cb, barcode_reads) 