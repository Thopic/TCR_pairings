#!/usr/bin/env bash
# Behaviour described in the README attached, in short create pairs from list of wells

cutoff=20
fdr=0.01
nbwells=5
logpvalue=-2.5

echo "Copy and extraction"
dir="$1"
cp src/{list_sequences.py,recreate_cells.py,pair} "$dir"
cp "$2"*.tsv.gz "$dir"
cd "$dir"
gunzip *.tsv.gz


echo "Removing unwanted columns"
awk -F "\t" -i inplace 'BEGIN {OFS = FS} {print $0, FILENAME}' *.tsv
awk -F"\t" 'BEGIN {OFS = FS} {print $1,$2,$3,$9,$11,$17,$23,$28,$29,$30,$31,$32,$33,$43,$NF}' TCRA.* > TCRA_all.tsv
awk -F"\t" 'BEGIN {OFS = FS} {print $1,$2,$3,$9,$11,$17,$23,$28,$29,$30,$31,$32,$33,$43,$NF}' TCRB.* > TCRB_all.tsv

echo "Sort"
mkdir tmp
sort -T tmp/ TCRA_all.tsv -o TCRA_all.tsv
sort -T tmp/ TCRB_all.tsv -o TCRB_all.tsv
rm -r tmp/

echo "Parsing sequences alpha"
./list_sequences.py TCRA_all.tsv sequences_alpha.tsv -w="$nbwells"
echo "Parsing sequences beta"
./list_sequences.py TCRB_all.tsv sequences_beta.tsv -w="$nbwells" --p_values="$logpvalue"

echo "Find the pairs ..."
awk -F"\t" 'BEGIN {OFS = FS} {print $1,$NF}' sequences_alpha.tsv > TCRA_short.tsv
awk -F"\t" 'BEGIN {OFS = FS} {print $1,$NF}' sequences_beta.tsv > TCRB_short.tsv
echo "    alpha-alpha"
./pair dvalues.tsv pvalues.tsv TCRA_short.tsv TCRA_short.tsv pairs_aa.tsv "$fdr"
echo "    beta-beta"
./pair dvalues.tsv pvalues.tsv TCRB_short.tsv TCRB_short.tsv pairs_bb.tsv "$fdr"
echo "    alpha-beta"
./pair dvalues.tsv pvalues.tsv TCRA_short.tsv TCRB_short.tsv pairs_ab.tsv "$fdr"

echo "Recreate the cells"
./recreate_cells.py sequences_alpha.tsv sequences_beta.tsv pairs_ab.tsv pairs_aa.tsv pairs_bb.tsv cells.tsv --cutoff="$cutoff"

# echo "Remove intermediate files"
# rm TCR*
# rm {list_sequences.py,recreate_cells.py,pair}

echo "Done."
