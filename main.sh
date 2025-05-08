#!/bin/bash

set -euo pipefail

# --- Configuration ---
gene_name="sxl" # replace with your gene of interest
protein_fasta="${gene_name}_protein.fasta" # protein sequence fasta file for gene of interest
fasta_dir="../zipped_fasta"
gtf_dir="../gtf_files"
output_dir="results/${gene_name}"

mkdir -p ./{bed_genes,fna_sizes,500_upstream/{plus/{fasta,canonical_tagteam/csv},minus/{fasta,canonical_tagteam/csv}}}
# The directory is as follows in which this is the file main.sh in the workding directory \gene-name
#   gene_name/
#      500_upstream/plus/
#          plus/
#              fasta/
#              canonical_tagteam/
#                csv/
#           minus/
#              fasta/
#              canonical_tagteam/
#                csv/
#      bed_genes/
#      fna_sizes/
#

export PATH=/programs/bedtools2-2.29.2/bin:$PATH
export PATH=/programs/bedops-2.4.35/bin:$PATH
export PATH=/programs/seqkit-0.15.0:$PATH


# --- Step 1: Run tblastn against all fna files ---
# Description: Finds the best hit for the gene of interest for each species. Determines where in the genome the best hit corresponds to.

echo "Running tblastn searches..."
for fasta_file in "$fasta_dir"/*.fna; do
    # Run tblastn and process the results in a single pipeline
    tblastn -query "$protein_fasta" -subject "$fasta_file" -outfmt 6 -evalue 1e-5 -max_target_seqs 10 |
        sort -k 11,11 -g |       # Sort by bitscore
        head -n 1 |              # Extract top hit
        awk '{
            OFS="\t";
            strand = ($9 <= $10) ? "+" : "-";
            start = ($9 <= $10) ? $9 : $10;
            end = ($9 <= $10) ? $10 : $9;
            print $2, "tblastn", "similarity", start, end, $12, strand, ".", "ID="$1";evalue="$11";bitscore="$12
        }' > $(basename "$fasta_file" .fna)_bh.gff

    # Perform intersection with the corresponding GTF file
    bedtools intersect -a "$gtf_dir"/$(basename "$fasta_file" .fna).gtf -b $(basename "$fasta_file" .fna)_bh.gff -wa \
        > $(basename "$fasta_file" .fna)_bh_overlap

    # Clean up intermediate files
    rm $(basename "$fasta_file" .fna)_bh.gff
done

echo "Isolating mRNA strands..."
# Isolate mRNA strands from the overlap files
for overlap_file in ./*_overlap; do
    grep -w "mRNA" "$overlap_file" > "${overlap_file}_mRNA"
    rm "$overlap_file"  # Remove the original overlap file after processing
done

echo "Modifying annotations..."
# Modify annotations in the mRNA files
for mRNA_file in ./*_mRNA; do
    file_a=""$gtf_dir"/$(basename "$mRNA_file" _bh_overlap_mRNA).gtf"
    file_b="$mRNA_file"
    
# Create a temporary file for modified output
    temp_file=$(mktemp)

    # Use awk to process both files, passing temp_file as a variable
    awk -v gene_name="$gene_name" -v temp_file="$temp_file" 'BEGIN { OFS="\t" }  # Set output field separator to tab
        # Step 1: Read all lines of file_b and store the 9th column values in an array
        NR == FNR { 
            match_column_b[$9]; 
            next 
        }

        # Step 2: For lines in file_a, check if the 9th column exists in file_b
        {
#	If	the	9th	column	matches,	add	the	suffix suffix ";name=real-gene-name;"
            if ($9 in match_column_b) {
            $9	=	$9";name=" gene_name ";"						
            }

            # Print the modified or original line to the temporary file
            print > temp_file
        }
    ' "$file_b" "$file_a"

    # Replace file_a with the modified content
    mv "$temp_file" "$file_a"
    rm "$mRNA_file"

done


for filename in "$gtf_dir"/*.gtf; do
    output_file="bed_genes/$(basename "$filename" .gtf)_genes.gtf"
    > "$output_file"
    echo "$filename"
    
    awk -v gene_name="$gene_name" '
    BEGIN { in_gene = 0; strand = ""; cds_line = ""; last_cds_line = "" }
    {
        if ($3 == "mRNA") {
            # Output the last CDS line if previous mRNA was on the - strand
            if (in_gene && strand == "-" && last_cds_line != "") {
                print last_cds_line > output_file
            }

            # Reset variables for a new mRNA block
            in_gene = 0
            cds_line = ""
            last_cds_line = ""
            strand = ""

#	Check	if	this	mRNA	is	"real-gene-name"		
if	(index($0,	"gene_name")	!=	0)	{			
                print "gene_name" > "/dev/stderr"
                in_gene = 1
                strand = $7
            }
        }
        
#	If	inside	a	"real-gene-name"	mRNA	annotation	and	the the line is a CDS
        else if (in_gene && $3 == "CDS") {
            if (strand == "+") {
                # For + strand, output the first CDS line found
                if (cds_line == "") {
                    cds_line = $0
                    print cds_line > output_file
                }
            } else if (strand == "-") {
                # For - strand, store the last CDS line encountered
                last_cds_line = $0
            }
        }
    }
    END {
        # After the final mRNA block, append the last CDS line if strand is -
        if (in_gene && strand == "-" && last_cds_line != "") {
            print last_cds_line > output_file
        }
    }
    ' output_file="$output_file" "$filename"
done

for gtf_file in bed_genes/*.gtf; do
    base=$(basename "$gtf_file" .gtf)
    
    # Convert GTF to BED, remove duplicates, and extract 500bp upstream regions
    gff2bed < "$gtf_file" | awk -v base="$base" '
        BEGIN { OFS = "\t" }
        !seen[$2, $3]++ {
            if ($6 == "-") {
                start = $3
                end = start + 500
                print $1, start, end, $4, $5, $6 > "500_upstream/minus/" base "_minus_500_upstream"
            } else if ($6 == "+") {
                end = $2
                start = end - 500
                if (start < 0) start = 0
                print $1, start, end, $4, $5, $6 > "500_upstream/plus/" base "_plus_500_upstream"
            }
        }
    '

    # Remove the original GTF file
    rm "$gtf_file"
done

for filename in "$fasta_dir"/*.fna; do
    faidx "$filename" -i chromsizes -o fna_sizes/$(basename "$filename" .fna).fna_sizes
done

# Define function to process strand
process_strand() {
    strand=$1  # "plus" or "minus"
    input_dir="500_upstream/$strand"

    for filepath in "$input_dir"/*_upstream; do
        echo "Processing $filepath"

        base=$(basename "$filepath")
        fasta_file="$input_dir/fasta/${base}.fasta"
        genome_fasta="$fasta_dir/$(basename "$filepath" _genes_${strand}_500_upstream).fna"
        tagteam_output="$input_dir/canonical_tagteam/${base}_canonical_tagteam"
        final_output="${tagteam_output}_fin"
        csv_output="$input_dir/canonical_tagteam/csv/${base}.csv"

        # Get FASTA
        bedtools getfasta -fo "$fasta_file" -fi "$genome_fasta" -bed "$filepath" -nameOnly

        # Locate motifs
        seqkit locate --ignore-case -p "CAGGTAG" -p "tAGGTAG" -p "CAGGcAG" "$fasta_file" > "$tagteam_output"

        # Remove duplicate coordinates
        awk '!seen[$2, $3]++' "$tagteam_output" > "$final_output"

        # Convert to CSV format
        sed 's/ \+/,/g' "$final_output" > "$csv_output"
    done
}

# Run for both strands
process_strand "minus"
process_strand "plus"

for filename in "$gtf_dir"/*; do
    input_file="$filename"
    output_file="$filename".2

    awk 'BEGIN { OFS="\t" }
    {
        # Check if "real-gene-name" is in any column
        for (i=1; i<=NF; i++) {
            if ($i ~ /real-gene-name/) {
                # Combine all columns from 9 onwards into column 9
                for (j=9; j<=NF; j++) {
                    $9 = $9 " " $j;
                }
                # Set NF to 9, effectively removing all columns after the 9th
                NF = 9;
                break;  # Exit the loop after modifying the columns
            }
        }
        print $0
    }' "$input_file" > "$output_file"
    mv "$output_file" "$input_file"
done

