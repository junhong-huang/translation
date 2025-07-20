import pandas as pd
import numpy as np
from collections import Counter

# Define codon order once at the module level
CODON_ORDER = [
    'TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG',
    'TAT', 'TAC', 'TGT', 'TGC', 'TGG',
    'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG',
    'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG',
    'ATT', 'ATC', 'ATA', 'ACT', 'ACC', 'ACA', 'ACG',
    'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG',
    'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG',
    'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG'
]

def load_codons_table(codons_file):
    codons_df = pd.read_csv(codons_file, sep='\t')
    codon_to_aa = dict(zip(codons_df['Codon'], codons_df['AA']))
    return codons_df['Codon'].tolist(), codon_to_aa

def load_expression_matrix(expression_file):
    expr_df = pd.read_csv(expression_file, sep='\t', index_col=0)
    
    # Filter out rows with '+' in index before processing
    expr_df = expr_df[~expr_df.index.str.contains('\+', regex=True)]
    
    gene_df = expr_df[expr_df.index.str.contains('__exon')].copy()
    gene_df.index = gene_df.index.str.replace('__exon', '')
    
    trna_df = expr_df[expr_df.index.str.contains('__tRNA')].copy()
    extracted = trna_df.index.str.extract(r'tRNA-([A-Za-z]+)-([A-Z]{3})-(\d+-\d+)')
    trna_df.index = [f'tRNA-{row[0]}-{row[2]}-{row[1]}' for row in extracted.values]
    
    # print("Gene Expression Matrix:\n", gene_df)
    # print("tRNA Expression Matrix:\n", trna_df)
    return gene_df, trna_df

def load_sequences(sequence_file):
    seq_df = pd.read_csv(sequence_file, sep='\t', header=None, 
                         names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 
                                'strand', 'thickStart', 'thickEnd', 'itemRgb', 
                                'blockCount', 'blockSizes', 'blockStarts', 'sequence'])
    seq_df['name'] = seq_df['name'].str.split('|').str[-1]
    return dict(zip(seq_df['name'], seq_df['sequence']))

def count_codons(sequence):
    sequence = sequence.upper()
    if len(sequence) % 3 != 0:
        sequence = sequence[:-(len(sequence) % 3)]
    return Counter(sequence[i:i+3] for i in range(0, len(sequence), 3))

def calculate_codon_demand(codons, codon_to_aa, expr_df, sequences):
    samples = expr_df.columns
    codon_demand = pd.DataFrame(0.1, index=CODON_ORDER, columns=samples)
    
    for gene in expr_df.index:
        if gene not in sequences:
            continue
        codon_counts = count_codons(sequences[gene])
        expr_values = expr_df.loc[gene]
        for codon in codon_counts:
            if codon in CODON_ORDER:
                codon_demand.loc[codon] += codon_counts[codon] * expr_values
    
    aa_groups = {}
    for codon in CODON_ORDER:
        aa = codon_to_aa.get(codon)
        aa_groups.setdefault(aa, []).append(codon)
    
    normalized_demand = codon_demand.copy()
    for aa, aa_codons in aa_groups.items():
        if aa == 'Stop':
            continue
        if aa_codons:
            max_cu = codon_demand.loc[aa_codons].max()
            for codon in aa_codons:
                normalized_demand.loc[codon] = codon_demand.loc[codon] / max_cu.replace(0, 1)
    
    return normalized_demand

def calculate_codon_frequency(codons, expr_df, sequences):
    genes = expr_df.index
    codon_freq = pd.DataFrame(0, index=genes, columns=CODON_ORDER, dtype=int)
    
    for gene in genes:
        if gene not in sequences:
            continue

        codon_counts = count_codons(sequences[gene])
        total_codons = sum(codon_counts.values())
        for codon in CODON_ORDER:
            codon_freq.loc[gene, codon] = codon_counts[codon]
        
        # print(codon_freq)
    
    return codon_freq

def main(codons_file, expression_file, sequence_file, output_file, freq_output_file, trna_output_file):
    codons, codon_to_aa = load_codons_table(codons_file)
    gene_df, trna_df = load_expression_matrix(expression_file)
    sequences = load_sequences(sequence_file)
    
    codon_demand = calculate_codon_demand(codons, codon_to_aa, gene_df, sequences)
    codon_demand.to_csv(output_file, sep='\t', float_format='%.1f')
    
    codon_freq = calculate_codon_frequency(codons, gene_df, sequences)
    codon_freq.to_csv(freq_output_file, sep='\t', float_format='%.1f')
    
    trna_df.to_csv(trna_output_file, sep='\t', float_format='%.1f')

if __name__ == "__main__":
    outdir="/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um"

    codons_file = "/public/home/huangjh/perlBioTools/Patho-DBiT/tRNA_viruses/data/codons_table.txt"
    sequence_file = "/public/home/huangjh/perlBioTools/Patho-DBiT/tRNA_viruses/data/hg38_CDS_sequence.bed"
    expression_file = f"{outdir}/expmat.tsv"
    
    output_file = f"{outdir}/translation/demand_codon_by_spot.tsv"
    freq_output_file = f"{outdir}/translation/frequency_mRNA_by_codon.tsv"
    trna_output_file = f"{outdir}/translation/expression_tRNA_by_spot.tsv"
    main(codons_file, expression_file, sequence_file, output_file, freq_output_file, trna_output_file)