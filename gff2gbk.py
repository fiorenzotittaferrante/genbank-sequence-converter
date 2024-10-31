from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from BCBio import GFF
import sys, os, warnings

warnings.filterwarnings("ignore") # Some warning were ignored, but they might be important


def find_files(fna_directory, gff_directory):
    """
    Locate and pair FNA and GFF files within specified directories.

    This function scans two directories to find FNA (genomic sequence) and GFF 
    (feature annotation) files. It pairs files by matching names and checks if 
    both formats exist for each genomic dataset. Only paired files are returned.

    Args:
        fna_directory (str): Directory containing .fna files.
        gff_directory (str): Directory containing .gff files.

    Returns:
        zip: An iterator that yields paired paths of FNA and GFF files.
    """

    fasta_files = {}
    gff_files = {}

    for filename in os.listdir(fna_directory):
        name, ext = os.path.splitext(filename)
        if ext == '.fna':
            fasta_files[name] = os.path.join(fna_directory, filename)
        
    for filename in os.listdir(gff_directory):
        name, ext = os.path.splitext(filename)
        if ext == '.gff':
            gff_files[name] = os.path.join(gff_directory, filename)

    tmp = {key.replace('.PATRIC', ''): key for key in gff_files.keys()}
    common_keys = fasta_files.keys() & tmp.keys() # Not all the FNA files have the corresponding GFF files

    paired_fasta_files = [fasta_files[key] for key in common_keys]
    paired_gff_files = [gff_files[tmp[key]] for key in common_keys]

    return zip(paired_fasta_files, paired_gff_files)


def preprocess_gff(gff_file):
    """
    Modify GFF file annotations for compatibility with parsing.

    This function processes a GFF file to modify entries where `Parent=gene-xxx` 
    is replaced with `Parent=gene_xxx` to ensure compatibility with `GFF.parse`.
    Temporary files are created for the modified GFF entries.

    Args:
        gff_file (str): Path to the GFF file to be preprocessed.

    Returns:
        None
    """

    with open(gff_file, 'r') as f:
        lines = f.readlines()

    processed_lines = []
    for line in lines:

        if line.startswith('accn|'):    # To avoid problem with PATRIC
            line = line[len('accn|'):]

        if line.startswith('NZ_'):      # To avoid problem with some ncbi files
            line = line[len('NZ_'):]

        if '\tCDS\t' in line:
            fields = line.strip().split('\t')

            attributes = fields[-1].split(';') # Last field.
            
            for i, attribute in enumerate(attributes):
                if attribute.startswith('Parent='):
                    attributes[i] = attribute.replace('-', '_')
                    break

            fields[-1] = ';'.join(attributes)
            modified_line = '\t'.join(fields)
            processed_lines.append(modified_line + '\n')
        else:
            processed_lines.append(line)

    abs_gff_file = os.path.abspath(gff_file)
    directory = os.path.dirname(abs_gff_file)
    tmp_gff_file = os.path.join(directory, 'tmp_' + os.path.basename(gff_file))

    with open(tmp_gff_file, 'w') as f:
        f.writelines(processed_lines)


def read_fasta(fasta_file):
    """
    Parse a FASTA file to retrieve sequence records.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        list: List of sequence records from the FASTA file.
    """

    with open(fasta_file, "r") as f:
        records = list(SeqIO.parse(f, "fasta"))

    return records


def read_gff(gff_file):
    """
    Parse a temporary GFF file and extract feature annotations.

    This function processes the GFF file, extracting features including 
    type, location, and qualifiers such as locus tags, product, and translation 
    tables. A temporary GFF file is generated and removed after parsing.

    Args:
        gff_file (str): Path to the GFF file.

    Returns:
        dict: Dictionary of parsed features with sequence IDs as keys.
        int: Translation table ID, defaulting to 1 if not specified in qualifiers.
    """

    preprocess_gff(gff_file)

    gff_features = {}

    abs_gff_file = os.path.abspath(gff_file)
    directory = os.path.dirname(abs_gff_file)
    tmp_gff_file = os.path.join(directory, 'tmp_' + os.path.basename(gff_file))

    with open(tmp_gff_file) as f:

        for rec in GFF.parse(f):
            features = []

            for feature in rec.features:
                feature_info = {
                    "type": feature.type,
                    "location": feature.location,
                    "strand": feature.location.strand,
                    "qualifiers": feature.qualifiers
                }

                if 'locus_tag' in feature.qualifiers:
                    feature_info['locus_tag'] = feature.qualifiers['locus_tag'][0]

                if 'product' in feature.qualifiers:
                    feature_info['product'] = feature.qualifiers['product'][0]

                if 'transl_table' in feature.qualifiers:
                    feature_info['transl_table'] = feature.qualifiers['transl_table'][0]
                else:
                    feature_info['transl_table'] = 1 # Standard table

                features.append(feature_info)
            gff_features[rec.id] = features

        os.remove(tmp_gff_file) # Remove temporary file generated.

    return gff_features, feature_info['transl_table']


def translate_cds(seq_record, features, translation_table):
    """
    Translate coding sequences (CDS) into protein sequences.

    This function iterates through CDS features in a sequence record, translating 
    DNA sequences to protein sequences using a specified translation table.

    Args:
        seq_record (SeqRecord): Sequence record containing the DNA sequence.
        features (list): List of features associated with the sequence record, 
                         each containing information such as type, location, 
                         strand, and qualifiers.
        translation_table (int): NCBI translation table ID for translating DNA to proteins.

    Returns:
        list: Updated list of features with translated proteins.
    """

    # Reverse iteration to avoid some problems
    for feature in features:    

        if feature['type'] == 'CDS':

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")

                cds_seq = feature['location'].extract(seq_record.seq).upper()

                # Ensure sequence length is a multiple of 3 to avoid translation warnings/errors.
                decimals = len(cds_seq)/3 - int(len(cds_seq)/3)

                if decimals > 0.5:
                    cds_seq = cds_seq + 'N'
                elif decimals > 0 and decimals <= 0.5:
                    cds_seq = cds_seq + 'NN'
                    
                protein_seq = cds_seq.translate(table=translation_table, to_stop=True)
                feature['qualifiers']['translation'] = str(protein_seq)

    return features


def create_genbank(seq_record, features, output_file):
    """
    Generate a GenBank file with annotated features.

    This function populates a sequence record with features (e.g., genes, CDS), 
    setting relevant qualifiers and annotations. It then writes the annotated 
    sequence record to an output GenBank file.

    Args:
        seq_record (SeqRecord): The sequence record to annotate with features.
        features (list): List of feature information, each containing location, type, and qualifiers.
        output_file (str): Path to the output GenBank (.gbk) file.

    Returns:
        None
    """

    seq_record.features = [] # Contenuto del gff (region, gene, CDS, ...) per quella sequenza

    for feature_info in features:
        location = feature_info['location']
        qualifiers = feature_info['qualifiers']

        # Check if the feature has a translation qualifier, and if so, ensure it is at least 6 characters long
        if 'translation' in qualifiers and len(qualifiers['translation']) < 6:
            print(f"\tSkipping feature with short translation: {qualifiers['translation']}")
            continue  # Skip this feature if the translation is too short

        feature = SeqFeature(
            location=location,
            type=feature_info['type'],
            qualifiers=qualifiers
        )
        seq_record.features.append(feature)

    if "annotations" not in seq_record.__dict__:
        seq_record.annotations = {}

    seq_record.annotations["molecule_type"] = "DNA"

    with open(output_file, "w") as output_handle:
        SeqIO.write(seq_record, output_handle, "genbank")


def main(fna_directory, gff_directory, out_directory):
    """
    Process FNA and GFF files, translate CDS sequences, and save annotations to GenBank files.

    This function scans specified directories to find pairs of FNA and GFF files, 
    extracts DNA sequences and annotations, and translates CDS sequences into proteins 
    based on specified translation tables. It outputs annotated GenBank files in the 
    designated directory and provides information on any CDS translations that resulted 
    in short or incomplete proteins.

    Args:
        fna_directory (str): Directory containing .fna files with genomic sequences.
        gff_directory (str): Directory containing .gff files with feature annotations.
        out_directory (str): Directory to save generated GenBank (.gbk) files.

    Returns:
        None
    """

    if not os.path.exists(out_directory):
        os.makedirs(out_directory)

    files = find_files(fna_directory, gff_directory)

    # For each strain
    for fasta_file, gff_file in files:

        fasta_records = read_fasta(fasta_file)
        gff_features, translation_table = read_gff(gff_file)

        print(f'File:\t{gff_file}\t\tTranslation table: {translation_table}')
        
        # For each DNA sequence
        for seq_record in fasta_records:

            if seq_record.id in gff_features:
                output_file = os.path.basename(gff_file[:-3])+'gbk'
                features = gff_features[seq_record.id]

                translated_features = translate_cds(seq_record, features, translation_table)
                create_genbank(seq_record, translated_features, os.path.join(out_directory, output_file))

            else:
                pass
                # print(f"Warning: No GFF features found for {seq_record.id}")


if __name__ == "__main__":
    fna_directory = sys.argv[1]
    gff_directory = sys.argv[2]
    out_directory = sys.argv[3]
    main(fna_directory, gff_directory, out_directory)
