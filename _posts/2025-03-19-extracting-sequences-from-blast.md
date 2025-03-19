---
layout: post
title: Extracting Sequences from Genomes Using BLAST Results
date: 2025-03-19 15:00:00
description: A practical guide to extracting sequences from genomes based on BLAST search results with proper handling of strand orientation.
tags: bioinformatics python blast genomics sequence-analysis
categories: tutorials
giscus_comments: true
related_publications: false
---

# Extracting Sequences from Genomes Using BLAST Results: A Practical Guide

When working with genomic data, researchers often need to extract specific sequences from a genome based on BLAST search results. This process requires careful handling, especially when dealing with sequence orientation (positive or negative strand). In this post, I'll share a Python script that automates this task while properly accounting for strand orientation.

## The Challenge

BLAST (Basic Local Alignment Search Tool) is commonly used to identify regions of similarity between biological sequences. When you get BLAST results, you'll have information about:

- Which parts of your query matched the genome
- Where these matches are located in the genome
- Whether the matches are on the forward (+) or reverse (-) strand

Manually extracting these sequences would be tedious and error-prone, especially when dealing with multiple hits or when sequences are on the negative strand and need to be reverse-complemented.

## The Solution: Automated Sequence Extraction

Here's a Python script that uses the Biopython library to extract sequences from a genome based on BLAST results:

    python
    #!/usr/bin/env python3
    import argparse
    from Bio import SeqIO
    from Bio.Seq import Seq

    def extract_sequences(blast_file, genome_file, output_file):
    # Parse genome file
    genome_dict = {}
    for record in SeqIO.parse(genome_file, "fasta"):
        genome_dict[record.id] = str(record.seq)
    
    # Store regions to extract
    regions_to_extract = []
    
    # Parse BLAST output file
    with open(blast_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            # Assuming BLAST output format: query_id, subject_id, identity, alignment_length, ...
            # Columns 9 and 10 are subject start and end positions
            # Column 12 is strand information (+/-)
            
            if len(parts) < 12:
                continue
                
            query_id = parts[0]
            subject_id = parts[1]
            start = int(parts[8])
            end = int(parts[9])
            strand = parts[11]
            
            # Ensure start position is less than end position
            if start > end:
                start, end = end, start
            
            regions_to_extract.append({
                'query_id': query_id,
                'subject_id': subject_id,
                'start': start,
                'end': end,
                'strand': strand
            })
    
    # Extract sequences and write to output file
    with open(output_file, 'w') as out:
        for region in regions_to_extract:
            if region['subject_id'] not in genome_dict:
                print(f"Warning: Could not find {region['subject_id']} in the genome")
                continue
            
            # Extract sequence from genome
            sequence = genome_dict[region['subject_id']][region['start']-1:region['end']]
            
            # If on negative strand, reverse complement the sequence
            if region['strand'] == '-':
                sequence = str(Seq(sequence).reverse_complement())
            
            # Write in FASTA format
            header = f">{region['query_id']}|{region['subject_id']}:{region['start']}-{region['end']}:       {region['strand']}"
            out.write(f"{header}\n{sequence}\n")

    def main():
    parser = argparse.ArgumentParser(description='Extract sequences from a genome based on BLAST results')
    parser.add_argument('-b', '--blast', required=True, help='Path to BLAST output file')
    parser.add_argument('-g', '--genome', required=True, help='Path to genome FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Path to output file')
    
    args = parser.parse_args()
    
    extract_sequences(args.blast, args.genome, args.output)
    
    if __name__ == "__main__":
    main()


## How It Works

The script performs the following steps:

1. **Load the genome into memory**: It reads the genome FASTA file and stores it as a dictionary for quick access.
2. **Parse the BLAST results**: It reads through the BLAST output file, extracting key information like sequence IDs, start/end positions, and strand orientation.
3. **Extract and process sequences**: For each hit, it extracts the corresponding sequence from the genome, reverse complementing it if it's on the negative strand.
4. **Output in FASTA format**: It writes the extracted sequences to a FASTA file with informative headers.

## Using the Script

To use this script, you'll need:
- Python 3
- The Biopython library (`pip install biopython`)
- A BLAST output file in tabular format (e.g., obtained with `-outfmt 6` option)
- A genome FASTA file

Run the script with:

bash
python extract_sequences.py -b blast_results.txt -g genome.fasta -o extracted_sequences.fasta


## Example Usage Scenario

Let's walk through a typical use case:

1. **Generate BLAST results**:
   
   bash
   blastn -query my_genes.fasta -db genome.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sstrand" -out blast_hits.txt
   

2. **Examine the BLAST output**:
   
   
   gene1   chromosome1    98.5    200    3    0    1    200    5000    5199    0.0    +
   gene2   chromosome2    99.0    150    1    1    5    154    3000    2851    0.0    -
   

3. **Run the extraction script**:
   
   bash
   python extract_sequences.py -b blast_hits.txt -g genome.fasta -o extracted_sequences.fasta
   

4. **Check the extracted sequences**:
   
   
   >gene1|chromosome1:5000-5199:+
   ATGCCTGAATTAGCTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGCCGTATCGTAGC...
   >gene2|chromosome2:2851-3000:-
   CTAGCTAGCAGCTGCTATTATCGCGCTAGCTAGCGCGCGTATTTTAGCGATCGTTAGCT...
   

Notice how the second sequence has been reverse-complemented because it matched on the negative strand.

## Important Considerations

- **BLAST output format**: The script assumes BLAST output in tabular format with specific column positions. Adjust the parsing logic if your format differs.
- **Coordinate systems**: BLAST uses 1-based coordinates, and the script handles this accordingly.
- **Memory usage**: For very large genomes, you might need to optimize the script to reduce memory consumption.
- **Error handling**: The script includes basic error checking but could be enhanced with more robust validation.

## Final Thoughts

This script addresses a common need in genomic analysis workflows. By automating the extraction of sequences from BLAST results while handling strand orientation correctly, it saves time and reduces the risk of errors.

For more complex scenarios, you might want to extend the script to handle additional features such as:
- Support for different BLAST output formats
- Extraction of flanking regions around hits
- Batch processing of multiple BLAST files
- Integration with other genomic analysis tools

Happy coding and sequence analyzing!
