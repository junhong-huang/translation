import gffutils
import pandas as pd
import subprocess
import os


#### python fetch_gtf_cds.py
gtf_file = "/public1/home/huangjh/annotation/RNAcentral/human/hg38.gencode.v30.annotation.gtf"
org = os.path.basename(gtf_file).split('.')[0]
species="human"
if org == "mm10":
    species = "mouse"
fasta_file = f"/public/home/huangjh/genome/{species}/{org}/WholeGenomeFasta/{org}.fa"
bed_output = f"/public/home/huangjh/perlBioTools/Patho-DBiT/tRNA_viruses/data/{org}_CDS.bed"
final_output = f"/public/home/huangjh/perlBioTools/Patho-DBiT/tRNA_viruses/data/{org}_CDS_sequence.bed"


# Create in-memory database from GTF file
db = gffutils.create_db(gtf_file, ":memory:")
bed12_rows = []
gene_cds = {}

for transcript in db.features_of_type('transcript'):
    # Extract GTF attributes
    transcript_id = transcript.attributes['transcript_id'][0]
    gene_id = transcript.attributes.get('gene_id', ['unknown'])[0]
    gene_name = transcript.attributes.get('gene_name', ['unknown'])[0]  # Include gene_name if available
    # Create name field with full GTF information (e.g., gene_id|transcript_id|gene_name)
    name = f"{gene_id}|{transcript_id}|{gene_name}"
    chrom = transcript.chrom
    strand = transcript.strand
    
    # Get start and stop codons
    start_codon = list(db.children(transcript, featuretype='start_codon'))
    stop_codon = list(db.children(transcript, featuretype='stop_codon'))
    
    if not start_codon or not stop_codon:
        continue
    
    start_codon_start = min([x.start for x in start_codon]) - 1
    stop_codon_end = max([x.end for x in stop_codon])
    
    # Get exons for the transcript
    exons = list(db.children(transcript, featuretype='exon', order_by='start'))
    if not exons:
        continue
    
    # Filter exons to include only CDS regions
    filtered_exons = []
    for exon in exons:
        exon_start = exon.start - 1
        exon_end = exon.end
        cds_start = max(exon_start, start_codon_start)
        cds_end = min(exon_end, stop_codon_end)

        if strand == '+':
            cds_start = max(exon_start, start_codon_start)
            cds_end = min(exon_end, stop_codon_end)
            if cds_end > cds_start:
                filtered_exons.append((cds_start, cds_end))
        else:
            cds_start = max(exon_start, stop_codon_end)
            cds_end = min(exon_end, start_codon_start)
            if cds_end > cds_start:
                filtered_exons.append((cds_start-3, cds_end+3))
    
    if not filtered_exons:
        continue
    
    # Calculate CDS length
    cds_length = sum(e[1] - e[0] for e in filtered_exons)
    
    # Keep the transcript with the longest CDS per gene
    if gene_id not in gene_cds or cds_length > gene_cds[gene_id][0]:
        chrom_start = min([e[0] for e in filtered_exons])
        chrom_end = max([e[1] for e in filtered_exons])
        block_count = len(filtered_exons)
        block_sizes = [e[1] - e[0] for e in filtered_exons]
        block_starts = [e[0] - chrom_start for e in filtered_exons]
        
        # Create BED12 row with comprehensive name field
        bed12_row = [
            chrom,
            chrom_start,
            chrom_end,
            name,  # Use combined name: gene_id|transcript_id|gene_name
            0,  # Score
            strand,
            start_codon_start,
            stop_codon_end,
            0,  # RGB
            block_count,
            ",".join(map(str, block_sizes)) + ",",
            ",".join(map(str, block_starts)) + ","
        ]
        gene_cds[gene_id] = (cds_length, bed12_row)

# Extract BED12 rows from gene_cds
bed12_rows = [row for _, row in gene_cds.values()]
bed12_df = pd.DataFrame(bed12_rows)
bed12_df.to_csv(bed_output, sep="\t", header=False, index=False)

# Run bedtools getfasta to extract sequences
bedtools_cmd = [
    "bedtools2", "getfasta",
    "-split", "-s", "-bedOut", "-nameOnly",
    "-fi", fasta_file,
    "-bed", bed_output
]
print(f"Running command: {' '.join(bedtools_cmd)}")

with open(final_output, "w") as outfile:
    subprocess.run(bedtools_cmd, stdout=outfile, check=True)