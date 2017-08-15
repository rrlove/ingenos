def filter_and_convert_genotypes(boolean_filter,genotypes,max_alleles=2,min_count=3):
    '''Take a scikit-allel genotype array and a filter saying which positions to include and return a set of allele counts ready for PCA.
    Max_alleles and min_count determine which genotypes are ultimately used in the PCA. The default selects biallelic positions where the minor allele is present at least 3 times.'''
    
    genotypes_subset = genotypes.subset(boolean_filter,)
    allele_count_subset = allel.AlleleCountsArray(genotypes_subset.count_alleles())
    pca_selection_bool = (allele_count_subset.max_allele() == max_alleles-1) & (allele_count_subset[:, :2].min(axis=1) > min_count)
    number_of_alternate_alleles = genotypes_subset.subset(pca_selection_bool).to_n_alt()[:]
    
    return number_of_alternate_alleles
