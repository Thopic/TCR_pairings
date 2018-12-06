# Supplementary to the paper "Genesis of αβ T-cell receptor"

This github is the companion to the paper entitled "Genesis of αβ T-cell receptor".

[Link to the paper (biorxiv)](https://www.biorxiv.org/content/early/2018/06/28/353128)


## Data 

We use the data obtained in the paper [Howie et al](https://www.ncbi.nlm.nih.gov/pubmed/26290413). Their dataset is accessible at [this address](http://s3-us-west-2.amazonaws.com/publishedproject-supplements/howie-2015-pairseq/index.html). 



The final results of the pairing algorithm, as well as some data needed for the analysis can be found in the `Data/` folder. 



## Analysis and generation of the figures of the paper

The different figures appearing in the paper and their generating code can be found in the *.ipynb (python notebook) files in the `Notes/` folder. All notebooks are jupyter notebooks and work with Python 3.6+ (and probably 3+). They can be displayed (but not run) by github. 

More precisely:
- `correlations_VJ.ipynb`: Contains the correlations between the V and J genes fragments for the different pairings, as well as the model used to model them. 
- `distances_entropies.ipynb`: Describes the distribution of CDR3's length for different type of sequences
- `mutual_information.ipynb`: Computes the mutual information between the different features, and analyzes the difference with the null. Also contains the selection model described in Methods.
- `pgen_distributions.ipynb`: Contains the probability of generation distribution for different types of sequences. The probability of generation are already computed in `../Datas/Pgen`, but they can be re-obtained using [IGoR](https://github.com/qmarcou/IGoR) if needed. 
- `proba_recombination.ipynb`: Generates the figures used in "Probability of recombination of ..." and "Fraction of cells with ...". 
- `proportion_shared_beta`: Show the proportion of shared betas.
- `sharing.ipynb`: Compute the selection factors, the number of sequences shared between individuals, the number of beta present in multiple clonotypes... Warning, generated sequences not included for size reasons. To obtain those one needs to download/install [IGoR](https://github.com/qmarcou/IGoR) and follow the instructions in the notebook.


## Pairing algorithm

The code use a variation on the PairSeq algorithm of [Howie et al.](https://www.ncbi.nlm.nih.gov/pubmed/26290413). It's contained in the `src/` folder.

This part is CPU and memory expensive, can be skipped if you are only interested in the data analysis.
Require a recent version of Python 3.x and numpy/scipy

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
