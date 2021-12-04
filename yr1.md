# Overview

The HPRC YR1 sample selection was designed to find representative candidates of 1000 Genome Project (abbreviated as 1KG from now on) subpopulations to be sequenced.

"Representative sample" is defined as being genetically similar to other samples of the same subpopulation, and differ from samples from other subpopulations. 
"Candidate" refers to the child in a 1KG trio, where both parents are present in the phase 3 genotyping release (i.e. the n=2504 vcf release). The ranking of a candidate is inferred using the parental data, specifically, the 1KG phase 3 genotyping results. 

# Implementation

## Pipeline input

Genotype - [genotyping of 1KG](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). Only autosomes are retained. 

Phenotype - [1KG sample info (spreadsheet)](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx) (`20130606_sample_info.xlsx`).

Auxiliary - subpopulation and superpopulation relationships from IGSR (`pheno_superpopulation.tsv`).

Auxiliary - cell line information from Coriell, July 2019 (`Coriell_Sample_List_072919.xls`). 

Auxiliary - [previous sequencing projects](https://docs.google.com/spreadsheets/d/1SBaJqLkQRRB92DEoxPHi7pVH9AESxHCbfYGWPXG7tLg/edit#gid=100551345) (`Previous Sequencing Projects.xlsx`)

Auxiliary - [trio existing data](https://docs.google.com/spreadsheets/d/1MJHYIqkH9vJ_1xEnFn-occJvgbcZENN_5SKDpoQKpOM/edit#gid=0) (`Trios Existing Data.xlsx`).

Auxiliary - previous sequencing project (`prevSeqProj.txt`)

The above mentioned files and some reformatted ones are available in the `input` folder, except for the vcf files. 

## Pipeline output

A list of samples sorted by their scores, along with parental heterozygosity and the sample's expected heterozygosity. See `output/het_with_ranking.1KGchildren.xlsx` (downloaded from the shared google sheet on Dec 12 2021).

## The pipeline

### Step 1: preprocessing

Input: genotype vcf files described above.

Output: 22 matrices of size (N, Mi), where N is the 1KG phase 3's sample size (i.e. 2504) and Mi is the number of loci that pass the minor allele frequency (MAF) requirement described below in the `i`th autosome.

Do: 

1. data lines of the input vcf file are retained if the MAF is at least 0.05. MAF is calculated as `min(AF, 1-AF)` where `AF` is extracted from the "AF" field of the vcf `INFO` field. Since 1KG phase 3 release's vcf files each containing all 2504 samples, the denominator of `AF` is the whole cohort. String to float conversion is python's `float()`.

2. if a data line passes the MAF requirement, it is converted to a integer array of length N. Each entry is calculated as the following: if the minor allele is REF, let the value be `genotype.count('0')`; otherwise, `2-genotype.count('0')`. To determine whether the minor allele is REF: true if `all_genotypes.count('0') < 2*N/2`, false otherwise.

### Step 2: dimension reduction

1. Each of the 22 feature matrices is normalized, goes through PCA and is reduced to size (N, 100). The normalization is implemented as:

```
import numpy as np
def normalize(C):
    # C is matrix of shape (nb_samples, nb_features)
    u = np.mean(C, axis=0,keepdims=True)
    p = u/2
    M = (C - u)/(p*(1-p))**0.5
    return M
```

PCA used the following implementation and parameters:

```
from sklearn.decomposition import PCA
PCA(n_components=100, svd_solver='arpack')
```

2. The feature matrices are concatenated to one matrix of shape (N, 2200) and is further reduced to (N, 100):

```
PCA(n_components=100, svd_solver='full', random_state=0)
```

The choice of solver and random states should have little impact on the final results.

3. As a sanity check, we can do a scatter plot using the first two dimensions of this features, or do a tSNE to visualize:

```
from sklearn.manifold import TSNE
tsne = TSNE(n_components=2, init='random', random_state=0,
           perplexity=20,
           early_exaggeration=30, 
           learning_rate=300,
           n_iter=2000)
```

The top-ranked samples in step3 should show up close to the center of each subpopulation clusters in these plots; 
note that step3 is done in the feature space, not the 2D representations.

We also do another sanity check by estimating heterozygosity of parents and children. 
Briefly, for samples where genotyping is available (i.e. the parents), 
we count heterozygous loci and divide the count by the total number of loci.
To infer heterozygosity for a child sample, we use the expectation to heterozygous loci counts. 
For example, we count 0.5 when parental genotypes are (a,a) and (a,b).

### Step 3: candidate selection

1. Ranking samples. We would try to find the most representative sample of a subpopulation by `argmax(i) (10*Dist(i, p, X[p]) - 1*Dist(i, ^p, X[^p]))`, where:

- `i` is the sample ID

- `X` is the feature matrix obtained from step2 (also, assume we know the subpopulation labels of each entry)

- `p` is the subpopulation label of `i`; `^p` is all subpopulation labels that are not `p`

- `X[p]` gives a subset of the X where all remaining samples have subpopulation labels belonging to `p`. Similar goes for `X[^p]`.

- `argmax(i)` finds the i that maximizes the target function

- `Dist(i, p, mat)` gives the average of the pairwise distance between `i` and all other data points in `mat`, normalized by `len(p)`.

The implementation: 

```
import pickle
import numpy as np
with open('100features_by_2ndPCA.pkl', 'rb') as file:
    X_pca = pickle.load(file)
nb_samples = X_pca.shape[0]

## load `phenos`, which is a list of non-redundant collection of subpopulation labels
## load `labels`, which is a list of length nb_samples

for pheno in phenos:
    mask = np.array([i for i in range(nb_samples) if labels[i]==pheno])
    mask_other = np.array([i for i in range(nb_samples) if labels[i]!=pheno])
    ins, outs, ratio = [], [], []
    mean_of_subpopulation = np.mean(X_pca[mask, :], axis=0).reshape(1, X_pca.shape[1])
    for idx_indv in mask:
        dist_in = np.linalg.norm(X_pca[idx_indv, :] - mean_of_subpopulation)
        dist_out = np.linalg.norm(X_pca[idx_indv, :].reshape(1, X_pca.shape[1]) - X_pca[mask_other, :])
        ratio.append(10*dist_in/1 + 1*dist_out/(len(phenos)-1) )
    result_scores = np.sort(ratio)
    ## then argsort the sample IDs using `result_scores` to get the ranked samples.
```

2. Ranking children of the trios:
a child is ranked as `best(parental_rank, maternal_rank)`.
Note that this might be ambiguous if samples of a trio have different subpopulation labels,
which was not an issue for the cohort we used here.

3. Apply the cell line optimal passage requirement:
remove candidates that are not `expansion==0` and `freeze_passage==2` in the `Coriell_Sample_List_072919.xls` data sheet.
If `freeze_passage` is empty, it is treated as N/A and dropped. 
We have 274 trios (822 samples) from 13 subpopulations after filtering. 

4. Remove candidates with existing or ongoing sequencing efforts.

5. Selection: ideally, we would like to select the same number of candidates from each subpopulation, and have equal 
number of candidates from both genders. We applied the following manual interventions: 
1) when gender is unbalanced (i.e. off by more than 1), try to swap in the next-best candidate of the less represented gender; do nothing if this is not possible.
2) if a subpopulation has less individuals than the desired sample selection size (i.e. all candidates are selected), their unused slots will be distributed to other unsaturated subpopulations. The latter choice is arbitrary but should have little impact on the overall results.

## Side notes as of 2021

1. The sample selection was done and maintained using the n=2504 1KG release; [the newer set that expanded the collection to 3202 samples](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/) was not considered here.

1. Due to the cell line optimal passage requirement and also the trio requirement, European samples have little representation in the final candidate list of this selection. Most of them had only suboptimal passages.

1. This approach focuses on selecting the "representative" candidates. It does not maximize diversity in terms of heterozygosity. Our justification of this choice (i.e. focusing on subpopulations), back in 2019, was the following (quote from the initial writeup): for the initial pilot, we should also make sure our methods work for populations with low diversity such as native Americans and Australians (these populations have much lower diversity than all 1000g samples). A simple solution is to choose one or two representative samples from each subpopulation. This will give us a spectrum of diversity. At the same time, from the popgen point of view, a sample we select [in this way could] better represents the larger subpopulation. This is critical when we make claims like "population X has more and longer SVs". 

1. I think it might have problems if intra-subpopulation diversity is very high (in other words, if a subpopulation could be further divided into clearly separated clusters, we should not fully rely on single subpopulation label and agglomerate everyone into one lump; 1KG cohorts have no such challenge, though). 

1. This approach replies on variant calling / genotyping. It would be difficult to recruite samples with different data types or processed by different pipelines.

# Change history of the candidate list (chronological)

Sep 2019: Initial candidate lists (size of 40 or 60) proposed for approval.

Jun 2020: Due to cell line issues, we expanded the sample selection from 60 to 120. Subpopulations ESN (AFR), MSL (AFR) and STU (SAS) had insufficient sample size.  Most subpopulations with sufficient sample size  were represented by 12~13 samples (a few 13s to fill up the selection to 120).

July 2020: Sample replacement: HG02068 and HG03008 were reported to have existing sequencing efforts. We proposed to replace them with HG02077 (different gender, but no more female samples were available) and HG03909 (HG03834 scored the next best choice, but HG03909 and HG03008 are of the same gender), respectively.

May 2021: We expanded the sample selection size from 120 to 160. Two more subpopulations now have all their samples selected: KHV (super: EAS), BEB (super: SAS). We still tried to maintain the gender balance when possible. As a side note, many populations currently have close to 20 samples, gender balance tends to need no intervention.

