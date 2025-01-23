# Data Files
This directory contains the synthetic data used.

Due to the large file sizes of real Visium data, we do not upload the Visium datasets to github. However, they are publically available at:

- Mouse Cortex stab wound: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226208`
- Spinal Cord Injury dataset: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184369`

For the synthetic data, there are 4 different setups.
- `checkerboard_simple_mvn_expr.csv` and `checkerboard_simple_mvn_meta.csv` contain the expression and metadata for the simulation of a checkerboard pattern where the underlying distribution is a multivariate normal.
- `checkerboard_spiral_expr.csv` and `checkerboard_spiral_meta.csv` contain the expression and metadata for the simulation of a checkerboard pattern where the underlying distribution is a low dimensional spiral pattern.
- `scurve_expr.csv` and `scurve_meta.csv` contain the expression and metadata for the simulation of a checkerboard pattern where the underlying distribution is a S-curve.
- `spiral_shape_mvn_synthetic_expression.csv` and `spiral_shape_mvn_synthetic_meta_df.csv` contain the expression and metadata for the simulation of a spiral pattern where the underlying distribution is a multivariate normal.
