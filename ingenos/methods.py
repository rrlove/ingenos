def construct_filter_bool(variants,expression):
    '''Take a set of variants and a positional expression on which to filter them, return a Boolean as to whether each variant/row of genotypes should be included post-filtering.'''

    selection_bool = variants.eval(expression)[:]
    return selection_bool

def construct_filter_expression(inversion,inversion_dict,region='expanded_breakpoints'):
    
    '''Construct an expression containing the desired SNPs.
    
    region takes three options: 
    
    inversion means SNPs are extracted from the whole inversion.
    
    breakpoints_only means SNPs are extracted from only between the breakpoints (this results in 0 SNPs for some inversions!)
    
    expanded_breakpoints means SNPs are extracted from the breakpoints + 1.5 Mb on either side.'''
    
    valid_region = {'inversion','breakpoints','expanded_breakpoints'}
    
    if region not in valid_region:
        raise ValueError("Region parameter must be one of the following: " % valid_region)
    
    proximal_start = inversion_dict[inversion].proximal_start
    proximal_end = inversion_dict[inversion].proximal_end
    distal_start = inversion_dict[inversion].distal_start
    distal_end = inversion_dict[inversion].distal_end
    
    if region == "breakpoints_only":
        
        expression = '( ( (POS > {prox_start}) & (POS < {prox_end}) ) | ( (POS > {dist_start} ) & (POS < {dist_end} ) ) )'.format(prox_start=proximal_start, prox_end=proximal_end, dist_start=distal_start, dist_end=distal_end)
        
    elif region == "inversion":
        
        expression = '( (POS > {prox_start}) & (POS < {dist_end}) )'.format(prox_start=proximal_start, dist_end=distal_end)
    
    elif region == "expanded_breakpoints":
        
        offset = 1500000
        
        expression = '( ( (POS > {prox_start}) & (POS < {prox_end}) ) | ( (POS > {dist_start} ) & (POS < {dist_end} ) ) )'.format(prox_start=proximal_start-offset, prox_end=proximal_end+offset, dist_start=distal_start-offset, dist_end=distal_end+offset)
    
    return expression

def downsample_genotypes(boolean_filter,downsample_rate=10):
    '''Downsample a set of genotypes for computational tractability.
    To be precise, this function stochastically shrinks the Boolean filter later used to pull out the genotypes.
    Shrinkage is approximate.'''
    
    downsampled_boolean = [random.randrange(100) < downsample_rate if value else value for value in boolean_filter]
    return downsampled_boolean
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
