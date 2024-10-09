#!/bin/bash
# The directory structure is as follows and this is the for_loop file in /workdir/$USER. Set the working directory to /workdir/$USER/sxl.
#$USER/
#   sxl/
#      500_upstream/plus/
#          plus/
#              fasta/
#              canonical_tagteam/
#                csv/
#           minus/
#              fasta/
#              canonical_tagteam/
#                csv/
#      5000_upstream/
#          plus/
#              fasta/
#              canonical_tagteam/
#                csv/
#           minus/
#              fasta/
#              canonical_tagteam/
#                csv/
#      bed/
#      bed_genes/
#      fna_sizes/
#      genes/

export PATH=/programs/bedops-2.4.35/bin:$PATH
for filename in bed/*.gtf; do
    gff2bed < "$filename" > $(basename "$filename" .gtf).bed
    mv $(basename "$filename" .gtf).bed bed_genes
done



# #download Dmel ref genome fasta and annotations (GTF) from genbank; convert from GTF -> BED
# export PATH=/programs/bedops-2.4.35/bin:$PATH
# for filename in *.gtf; do
# # account for some .gtf files missing "transcript_id", add empty string, then convert GTF -> BED
#     #awk '{ if ($0 ~ "gene_id") print $0; else print $0 "; gene_id \"" $1 "\";"}' "$filename" > "$filename".txt
#     awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' "$filename" > "$filename".txt
# done
# | gtf2bed - > bed/$(basename "$filename" .gtf).bed
#get bed file for only start codons and sxl gene

for filename in bed/*.bed ; do
    #isolate x chromosome
    grep -E "NW_022587371.1|NC_057931.1|NC_004354.4|NC_046683.1" "$filename" > genes/$(basename "$filename" .bed)_genes_og.bed
    grep -E genes/$(basename "$filename" .bed)_genes_og.bed > bed_genes/$(basename "$filename" .bed)_genes.bed
done

#get lengths of contigs
#.fna.fai column 1: name of contig column 2: number of bases in contig column 3: byte index of file where contig sequence begins column 4: bases per line in fasta file column 5: bytes per line in fasta file
for filename in *.fna; do
    faidx "$filename" -i chromsizes -o fna_sizes/$(basename "$filename" .fna).fna_sizes
done

#get 500 and 5000 bp upstream of start codon (takes into account negative or positve strand)
export PATH=/programs/bedtools2-2.29.2/bin:$PATH
#loop over .bed files in bed_genes folder
for filename in bed_genes/*.bed; do
    bedtools slop -i "$filename" -g fna_sizes/$(basename "$filename" _genes.bed).fna_sizes -s -l 500 -r 0 > 500_upstream/minus/$(basename "$filename" .bed)_minus_500_upstream
    bedtools slop -i "$filename" -g fna_sizes/$(basename "$filename" _genes.bed).fna_sizes -s -l 5000 -r 0 > 5000_upstream/minus/$(basename "$filename" .bed)_minus_5000_upstream
    #seperate plus and minus
    grep '+' 500_upstream/minus/$(basename "$filename" .bed)_minus_500_upstream > 500_upstream/plus/$(basename "$filename" .bed)_plus_500_upstream
    sed -i '/+/d' 500_upstream/minus/$(basename "$filename" .bed)_minus_500_upstream
    #5000
    grep '+' 5000_upstream/minus/$(basename "$filename" .bed)_minus_5000_upstream > 5000_upstream/plus/$(basename "$filename" .bed)_plus_5000_upstream
    sed -i '/+/d' 5000_upstream/minus/$(basename "$filename" .bed)_minus_5000_upstream
done

#loop over files in 500_upstream directory
#ending in _upstream as to not include folders
for filename in 500_upstream/minus/*_upstream; do
#get fasta for 500 bp windows
    bedtools getfasta -fo "$filename".fasta -fi $(basename "$filename" _genes_minus_500_upstream).fna -bed "$filename" -nameOnly
    mv "$filename".fasta 500_upstream/minus/fasta
#find motif
    export PATH=/programs/seqkit-0.15.0:$PATH
    seqkit locate --ignore-case -p "CAGGTAG" -p "tAGGTAG" -p "CAGGcAG" 500_upstream/minus/fasta/$(basename "$filename").fasta >"$filename"_canonical_tagteam
    mv "$filename"_canonical_tagteam 500_upstream/minus/canonical_tagteam
#remove duplicates in file
    awk '!visited[$0]++' 500_upstream/minus/canonical_tagteam/$(basename "$filename")_canonical_tagteam > 500_upstream/minus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin
#turn file into .csv format
    sed 's/ \+/,/g' 500_upstream/minus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin > 500_upstream/minus/canonical_tagteam/csv/$(basename "$filename").csv
done

for filename in 500_upstream/plus/*_upstream; do
#get fasta for 500 bp windows
    bedtools getfasta -fo "$filename".fasta -fi $(basename "$filename" _genes_plus_500_upstream).fna -bed "$filename" -nameOnly
    mv "$filename".fasta 500_upstream/plus/fasta
#find motif
    export PATH=/programs/seqkit-0.15.0:$PATH
    seqkit locate --ignore-case -p "CAGGTAG" -p "tAGGTAG" -p "CAGGcAG" 500_upstream/plus/fasta/$(basename "$filename").fasta >"$filename"_canonical_tagteam
    mv "$filename"_canonical_tagteam 500_upstream/plus/canonical_tagteam
#remove duplicates in file
    awk '!visited[$0]++' 500_upstream/plus/canonical_tagteam/$(basename "$filename")_canonical_tagteam > 500_upstream/plus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin
#turn file into .csv format
    sed 's/ \+/,/g' 500_upstream/plus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin > 500_upstream/plus/canonical_tagteam/csv/$(basename "$filename").csv
done


#loop over files in 5000_upstream directory
#ending in _upstream as to not include folders
for filename in 5000_upstream/minus/*_upstream; do
#get fasta for 5000 bp windows
    bedtools getfasta -fo "$filename".fasta -fi $(basename "$filename" _genes_minus_5000_upstream).fna -bed "$filename" -nameOnly
    mv "$filename".fasta 5000_upstream/minus/fasta
#find motif
    export PATH=/programs/seqkit-0.15.0:$PATH
    seqkit locate --ignore-case -p "CAGGTAG" -p "tAGGTAG" -p "CAGGcAG" 5000_upstream/minus/fasta/$(basename "$filename").fasta >"$filename"_canonical_tagteam
    mv "$filename"_canonical_tagteam 5000_upstream/minus/canonical_tagteam
#remove duplicates in file
    awk '!visited[$0]++' 5000_upstream/minus/canonical_tagteam/$(basename "$filename")_canonical_tagteam > 5000_upstream/minus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin
#turn file into .csv format
    sed 's/ \+/,/g' 5000_upstream/minus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin > 5000_upstream/minus/canonical_tagteam/csv/$(basename "$filename").csv
done

for filename in 5000_upstream/plus/*_upstream; do
#get fasta for 5000 bp windows
    bedtools getfasta -fo "$filename".fasta -fi $(basename "$filename" _genes_plus_5000_upstream).fna -bed "$filename" -nameOnly
    mv "$filename".fasta 5000_upstream/plus/fasta
#find motif
    export PATH=/programs/seqkit-0.15.0:$PATH
    seqkit locate --ignore-case -p "CAGGTAG" -p "tAGGTAG" -p "CAGGcAG" 5000_upstream/plus/fasta/$(basename "$filename").fasta >"$filename"_canonical_tagteam
    mv "$filename"_canonical_tagteam 5000_upstream/plus/canonical_tagteam
#remove duplicates in file
    awk '!visited[$0]++' 5000_upstream/plus/canonical_tagteam/$(basename "$filename")_canonical_tagteam > 5000_upstream/plus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin
#turn file into .csv format
    sed 's/ \+/,/g' 5000_upstream/plus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin > 5000_upstream/plus/canonical_tagteam/csv/$(basename "$filename").csv
done