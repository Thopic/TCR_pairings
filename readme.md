# Supplementary to the paper XXXXX


## Find the pairs using a variation on PairSeq

CPU and memory expensive, can be skipped if you are only interested in the data analysis.
Necessit a recent version of Python 3.x and numpy/scipy

For a given experiment, in the current directory, 
Directly use the bash file (default options)
```bash
mkdir dir_exp/
./pair_sequences.sh dir_exp/ /path/to/sequences/
```
Or:
- Copy and extract the files associated with the different wells of the experiment.
  There is a different file for each well and each sequence type, of the form "TCRA...5" or "TCRB...87"
```bash
mkdir dir_exp/ 
cp {list_sequences.py,recreate_cells.py,pair}  dir_exp/
cp /path/to/sequences/*.tsv.gz dir_exp/
cd dir_exp
gunzip *.tsv.gz
```
- Remove unwanted columns:
[Keep: sequence, amino-acid, number of copies, length of cdr3, V/D/Jgenes names, number of insertions/deletions, status]
```bash
awk -F "\t" -i inplace 'BEGIN {OFS = FS} {print $0, FILENAME}' *.tsv
awk -F"\t" 'BEGIN {OFS = FS} {print $1,$2,$3,$9,$11,$17,$23,$28,$29,$30,$31,$32,$33,$43,$NF}' TCRA.* > TCRA_all.tsv
awk -F"\t" 'BEGIN {OFS = FS} {print $1,$2,$3,$9,$11,$17,$23,$28,$29,$30,$31,$32,$33,$43,$NF}' TCRB.* > TCRB_all.tsv
```
- Sort the file and reunite the sequences (-w : minimal number of wells to keep a sequence). 
The argument -w indicate the minimal number of wells a sequence needs to be in to be conserved.
The argument -pvalues create a file necessary for pre-pairing, here it is possible to indicate
the wanted precision.
```bash
sort TCRA_all.tsv -o TCRA_all.tsv
sort TCRB_all.tsv -o TCRB_all.tsv
./list_sequences.py TCRA_all.tsv sequences_alpha.tsv -w4
./list_sequences.py TCRB_all.tsv sequences_beta.tsv -w4 --p_values -4.0
```
- Find the potential pairs
```bash
awk -F"\t" 'BEGIN {OFS = FS} {print $1,$NF}' sequences_alpha.tsv > TCRA_short.tsv
awk -F"\t" 'BEGIN {OFS = FS} {print $1,$NF}' sequences_beta.tsv > TCRB_short.tsv
./pair dvalues.tsv pvalues.tsv TCRA_short.tsv TCRB_short.tsv pairs_ab.tsv 0.01
./pair dvalues.tsv pvalues.tsv TCRA_short.tsv TCRA_short.tsv pairs_aa.tsv 0.01
./pair dvalues.tsv pvalues.tsv TCRB_short.tsv TCRB_short.tsv pairs_bb.tsv 0.01
```
- Final vetting (will ask for the cutoff value each time)
```bash
./vetting.py sequences_alpha.tsv sequences_beta.tsv TCRAB_pairs.tsv pairs_ab.tsv --fdr=0.01
./vetting.py sequences_alpha.tsv sequences_alpha.tsv TCRAA_pairs.tsv pairs_aa.tsv --fdr=0.01 --cutoff=5
./vetting.py sequences_beta.tsv sequences_beta.tsv TCRBB_pairs.tsv pairs_bb.tsv --fdr=0.01 --cutoff=5
```
- Recreate the cells
```bash
./recreate_cells.py sequences_alpha.tsv sequences_beta.tsv pairs_ab.tsv pairs_aa.tsv pairs_bb.tsv cells.tsv --cutoff=5
```
- Remove intermediate files
``` bash
rm TCR*; rm TCR*;
```


## Analysis of the results

The different figures appearing in the paper and their generating code can be found in the *.ipynb (python notebook) files in the folder Notes.
