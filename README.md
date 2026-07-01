# generalized-independence-test

This repository contains a lightweight R implementation of the generalized independence test, together with a real data example using GTEx normalized expression data.

## Testing Problem

Given two random variables, $X$ and $Y$, with marginal distributions $P_X$ and $P_Y$ over the spaces $\mathcal{X}$ and $\mathcal{Y}$, respectively, and their joint distribution $P_{XY}$ on $\mathcal{X} \times \mathcal{Y}$, the goal is to test the hypothesis:

$$
H_0: P_{XY}=P_X P_Y \quad \text{ versus } \quad H_1: P_{XY} \neq P_X P_Y
$$

using paired samples $\{(X_i, Y_i)\}_{i=1}^n$ drawn independently and identically from $P_{XY}$.

In the code, `X` and `Y` should be paired data matrices with the same number of rows:

- `X` is an `n x p1` matrix.
- `Y` is an `n x p2` matrix.
- Rows correspond to the `n` paired samples.
- Columns correspond to dimensions or features.

## Large File Storage

The GTEx data files in `data/gtex/` are tracked with Git LFS because one file exceeds GitHub's normal file size limit.

If you clone this repository and need the full data locally, make sure Git LFS is installed and then run:

```bash
git lfs pull
```

## Files

- `R/GIT_test.R`: implementation of the generalized independence test.
- `examples/real_data_example.R`: real data example using breast mammary tissue and EBV-transformed lymphocyte expression data.
- `data/gtex/Breast_Mammary_Tissue_normalized_expression_final.txt`
- `data/gtex/Cells_EBV-transformed_lymphocytes_normalized_expression_final.txt`

## Code

The main user-facing function is:

```r
GIT_test(X, Y, k = floor(sqrt(nrow(X))), robustgraph = TRUE, lambda = 0.3)
```

The function uses only base R. No additional R packages are required.

### Arguments

- `X`, `Y`: paired numeric data matrices, or objects coercible to matrices, with the same number of rows. `X` is `n x p1` and `Y` is `n x p2`; rows are paired samples and columns are dimensions/features.
- `k`: number of nearest/farthest neighbors used to build the graph. The default is `floor(sqrt(nrow(X)))`.
- `robustgraph`: whether to use the robust graph construction. The default is `TRUE`.
- `lambda`: degree-regularization parameter for the robust graph. The default is `0.3`.
- `symmetric`: whether to symmetrize directed graph ranks. The default is `TRUE`.

### Output

`GIT_test()` returns a named numeric vector:

```r
overall        G1        G2        G3        G4
```

`overall` is the generalized p-value, and `G1`, `G2`, `G3`, and `G4` are the four component p-values.

## Simulated Example

```r
source("R/GIT_test.R")

set.seed(1)
X <- matrix(rnorm(100), ncol = 2)
Y <- matrix(rnorm(100), ncol = 2)

out <- GIT_test(X, Y)
print(out)
```

## Real Data Example

Run the included real data example from the repository root:

```r
source("examples/real_data_example.R")
```

The example reads the two GTEx normalized expression files, keeps the shared sample columns, removes the first four metadata columns, transposes the expression matrices, and then runs:

```r
set.seed(1)
result <- GIT_test(X, Y)
print(result)
```

## Paper

This implementation is for the method introduced in:

Mingshuo Liu, Doudou Zhou, and Hao Chen. *Generalized Independence Test for Modern Data*. arXiv:2409.07745. https://arxiv.org/abs/2409.07745

## Robust Graph Construction

The robust graph construction used in this implementation follows:

Yejiong Zhu and Hao Chen. *Mitigating dimensionality effects with robust graph constructions for testing*. arXiv:2307.15205. https://arxiv.org/abs/2307.15205

The default value `lambda = 0.3` is used as suggested by the robust graph authors.

## Data Source

The data processing in this repository is related to the following study:

Khunsriraksakul C, McGuire D, Sauteraud R, Chen F, Yang L, Wang L, Hughey J, Eckert S, Weissenkampen JD, Shenoy G, et al. *Integrating 3D genomic and epigenomic data to enhance target gene discovery and drug repurposing in transcriptome-wide association studies*. Nature Communications. 2022;13(1):3258.

The published study used GTEx v7 data. The files in this repository are GTEx v8 data processed by the same research group.

## Acknowledgment

We thank Dr. Lida Wang for kindly providing the processed GTEx v8 data used in this repository.

## Data Structure

Each data file keeps the original GTEx-style tabular structure.

- Column 1: `gene_id`
- Column 2: `chromosome`
- Column 3: `start`
- Column 4: `end`
- Columns 5 and onward: normalized expression values for samples

In the analysis code, the first four columns are metadata columns and are removed before constructing the expression matrices. The analysis uses the expression values from column 5 onward.
