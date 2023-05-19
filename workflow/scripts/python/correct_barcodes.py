import gzip
import re
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

input_file    = snakemake.input[0]
output_fastq  = snakemake.output["corrected_fq"]
feature_ref   = snakemake.output["feature_ref"]
color         = snakemake.wildcards.larry_color
read          = snakemake.wildcards.read_fb

#------------------------------------------------------------------------------------------------------------------------------------
# Declare functions
#------------------------------------------------------------------------------------------------------------------------------------

def create_feature_ref(reference_dict, color, read, feature_ref):
    """
    Generate a dictionary with barcode sequences and their corresponding larry color to 
    generate the feature reference csv file required by cellranger.

    It takes as input a dictionary with read cluster id's as keys and sequences (corresponding to the
    reference barcode of the clusetr) as values.
    """
    extra_nts = r'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
    bc_dict = {str(key): color for key in reference_dict.values()}

    bc_df                 = pd.DataFrame.from_dict(bc_dict, orient='index', columns = ['id'])
    bc_df                 = bc_df.rename_axis("sequence").reset_index()
    bc_df['sequence']     = bc_df['sequence'].str.replace(extra_nts, '', regex=True).astype('str')
    bc_df['name']         = bc_df.groupby('id').cumcount() + 1
    bc_df['name']         = bc_df['id'] + "_" + bc_df['name'].astype(str)
    bc_df['id']           = bc_df['name']
    bc_df['read']         = read
    bc_df['pattern']      = extra_nts + "(BC)"
    bc_df['feature_type'] = "Custom"
    
    bc_df.to_csv(feature_ref, index = False)


def get_reference_barcodes(input_fastq):
    """
    This function will parse a fastq processed by UMIcollapse https://github.com/Daniel-Liu-c0deb0t/UMICollapse.
    By parsing the fastq, it will create a dictionary with read cluster id's as keys and the sequences of the
    reference barcodes of the cluster as values.

    It also returns a list containing all the fastq records to be used later and avoid having to parse again the fastq.
    """

    fastq_entries     = []
    cluster_dict_id   = dict()
    reference_pattern = re.compile(r'cluster_size=')

    with gzip.open(input_fastq, 'rt') as input_handle:
        for title, seq, qual in FastqGeneralIterator(input_handle):
            fastq_entries += [[title, seq, qual]]
            match = reference_pattern.search(title)
            if match:
                # Create a dictionary with the cluster id and the corresponding consensus sequence
                cluster_id = re.sub('^.* cluster_id=([0-9]+).*', r'\1', title)
                cluster_dict_id[cluster_id] = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" + seq

    return(cluster_dict_id, fastq_entries)


def write_corrected_fastq(output_fastq, reference_dict, fastq_entries):
    """
    This function will correct the sequences from the fastq file by substituting them by the reference sequence 
    of every read cluster and write a new fastq file.

    Takes as input a dictionary with read cluster id's as keys and the sequences of the
    reference barcodes of the cluster as values and a list containing all the records of the
    fastq file. Using the dictionary, it will update all the sequences from the list based on the reference sequence
    present in the dictionary.
    """

    with gzip.open(output_fastq, 'wt') as output_handle:
        for title, seq, qual in fastq_entries:
        # Get cluster id for read, use that cluster id to access cluster_dict_id dict
        # which contains the reference sequence. This way I generate a second dict in which
        # Every read id is a associated to a specific sequence.
        # Remove the tag from umicollapse from read_id (title)
            cluster_id = re.sub('^.* cluster_id=([0-9]+).*', r'\1', title)
            seq = reference_dict[cluster_id]
            title = re.sub("\scluster_id.*$", "", title)
            qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" + qual # Fake qscore for the extra 50 fake T's
            _ = output_handle.write(f"@{title}\n{seq}\n+\n{qual}\n")


#------------------------------------------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------------------------------------------
# Create dict with reference barcodes, store fastq in list
reference_dict, fastq_records = get_reference_barcodes(input_file)

# Correct fastq file with reference barcodes and save output fastq
write_corrected_fastq(output_fastq, reference_dict, fastq_records)

# Generate a feature ref csv from the dictionary of reference barcodes
create_feature_ref(reference_dict, color, read, feature_ref)