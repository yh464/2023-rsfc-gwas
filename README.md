This is the section of project-specific scripts for the resting-state functional MRI GWAS paper. Please combine this with the [GWAS pipeline](https://github.com/yh464/gwas_pipeline).
Usage:
1. Use `subj_list.py` to extract individuals with valid fMRI data.
2. `pheno.py` should be put into the same folder as `pheno_batch.py` in `gwas_pipeline` and batch run over all subjects found in `subj_list.py`.
3. `pheno_decomp.py` should be used to calculate the number of effectively independent phenotypes for multiple testing correction.