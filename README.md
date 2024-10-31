# Translate Sequences from FNA and GFF Files

This script processes DNA sequences from `.fna` and `.gff` files, translating coding DNA sequences (CDS) into protein sequences and saving them in GenBank format. It is designed to work with files where GFF annotations might need slight modifications for compatibility. Additionally, the script handles partial codons by adding "N" nucleotides where necessary to ensure proper translation.

## Features
- Identifies coding sequences (CDS) and translates them to proteins based on specified genetic translation tables.
- Cleans and preprocesses GFF files for compatibility, especially when working with PATRIC or NCBI files.
- Filters out any CDS with translations shorter than 6 amino acids.
- Outputs GenBank files with annotated features and translations.

## Requirements
This script requires Python and the following libraries:
- **Biopython**: `pip install biopython`
- **BCBio-gff**: `pip install bcbio-gff`

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/translate_sequences.git
   cd translate_sequences
   ```
2. Install the required Python libraries.

## Usage
To use the script, run the following command:

```bash
python3 gff2gbk.py fna/ gff/ gbk/
```

The script expects three command-line arguments:
1. The path to the directory containing `.fna` files (FASTA format with DNA sequences).
2. The path to the directory containing `.gff` files (GFF annotations for the sequences).
3. The output directory where the GenBank files will be saved.

Each `.fna` file with a corresponding `.gff` file will produce a `.gbk` (GenBank) file in the specified output directory.

### Input/Output Details
- Input files:
  - `.fna` files: FASTA files containing DNA sequences.
  - `.gff` files: GFF files containing CDS annotations, which may need preprocessing for compatibility.
- Output files:
  - `.gbk` files: GenBank files with features and translated protein sequences.

## Code overview
### Functions
- `find_files(fna_directory, gff_directory)`: Identifies and pairs `.fna` and `.gff` files by name.
- `preprocess_gff(gff_file)`: Cleans and preprocesses the GFF file to make it compatible with BioPython.
- `translate_cds(seq_record, features, translation_table, frequencies)`: Translates CDS sequences into proteins and manages translation table usage.
- `create_genbank(seq_record, features, output_file)`: Creates the final GenBank file, including a filter to exclude sequences shorter than 6 amino acids.

### Example Workflow
1. Read the FNA and GFF files.
2. Preprocess the GFF files to standardize the annotation format.
3. Translate CDS sequences, applying the specified translation table and handling partial codons.
4. Filter short translations and save the annotated sequences to GenBank format.
