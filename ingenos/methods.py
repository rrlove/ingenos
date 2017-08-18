def construct_filter_expression(name,inversion_dict,
                                buffer=1500000,whole_inversion=False):
    '''Construct an expression describing the desired SNPs.
    
    buffer describes the # of bp around each breakpoint to include.
    
    whole_inversion specifies whether only the breakpoints and the buffer,
    or the actual inversion, should be used.'''
    
    proximal_start = inversion_dict[name].proximal_start - buffer
    proximal_end = inversion_dict[name].proximal_end + buffer
    distal_start = inversion_dict[name].distal_start - buffer
    distal_end = inversion_dict[name].distal_end + buffer
    
    if not all(item > 0 for item in [proximal_start,proximal_end,distal_start,distal_end]):
        raise ValueError("An output coordinate is negative, please adjust buffer")

    if whole_inversion:
        
        expression = '( POS > {prox_start} & POS < {dist_end} )'.format(
            prox_start=proximal_start,
            dist_end=distal_end)
        
    else:
        
        expression = '( ( POS > {prox_start} & POS < {prox_end} ) | ( POS > {dist_start} & POS < {dist_end} ) )'.format(
            prox_start=proximal_start,
            prox_end=proximal_end,
            dist_start=distal_start,
            dist_end=distal_end)
        
    return(expression)

def filter_and_convert_genotypes(boolean_filter,genotypes,max_alleles=2,min_count=3):
    '''Take a scikit-allel genotype array and a filter saying which positions to include and return a set of allele counts ready for PCA.
    Max_alleles and min_count determine which genotypes are ultimately used in the PCA. The default selects biallelic positions where the minor allele is present at least 3 times.'''
    
    genotypes_subset = genotypes.subset(boolean_filter,)
    allele_count_subset = allel.AlleleCountsArray(genotypes_subset.count_alleles())
    pca_selection_bool = (allele_count_subset.max_allele() == max_alleles-1) & (allele_count_subset[:, :2].min(axis=1) > min_count)
    number_of_alternate_alleles = genotypes_subset.subset(pca_selection_bool).to_n_alt()[:]
    
    return number_of_alternate_alleles

def prune_by_LD(number_of_alternate_alleles,window_size=1000,step_size=100,r2=0.2):
    '''Take an array of the number of alternate alleles and return a smaller array pruned to remove SNPs in LD with each other.'''

    pruned_Bool = allel.locate_unlinked(number_of_alternate_alleles, window_size, step_size, r2)
    pruned = number_of_alternate_alleles[pruned_Bool]
    
    return pruned
