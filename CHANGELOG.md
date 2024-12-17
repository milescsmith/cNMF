# 1.9.999: (2.0 beta)

## Added:
- Partially type hinted

## Changed:
- Replacing as much as I can of numpy with jax.numpy
- Simplified a number of functions
    - Removing unnecessary conversions of arrays to `pandas.Series`
- Combined `cnmf.get_highvar_genes` and `cnmf.get_highvar_genes_sparse`
