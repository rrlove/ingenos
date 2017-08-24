import allel
import warnings

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
        
        expression = '( (POS > {prox_start}) & (POS < {dist_end}) )'.format(
            prox_start=proximal_start,
            dist_end=distal_end)
        
    else:
        
        expression = '( ( (POS > {prox_start}) & (POS < {prox_end}) ) | ( (POS > {dist_start}) & (POS < {dist_end}) ) )'.format(
            prox_start=proximal_start,
            prox_end=proximal_end,
            dist_start=distal_start,
            dist_end=distal_end)
        
    return(expression)

def filter_and_convert_genotypes(genotypes,sites_boolean=None,samples_boolean=None,max_alleles=2,min_count=3):
    '''Filter a genotype array based on booleans of sites and samples to include.
    
    Further filter genotypes based on allele count data.
    
    Return a set of alternate allele counts ready for PCA, and the allele counts filter.
    '''
    
    if not all(item > 0 for item in [max_alleles,min_count]):
        raise ValueError("Max alleles and minimum allele count must be greater than 0 ")
    
    if sites_boolean is not None:
        
        if not len(sites_boolean) == genotypes.shape[0]:
            raise ValueError("Length of sites filter does not match length of genotypes")
    
    if samples_boolean is not None:
        
        if not len(samples_boolean) == genotypes.shape[1]:
            raise ValueError("Length of samples filter does not match length of genotypes")
        
    if sites_boolean is not None and samples_boolean is None:

        genotypes_subset = genotypes.subset(sel0=sites_boolean)
        
    if samples_boolean is not None and sites_boolean is None:
        
        genotypes_subset = genotypes.subset(sel1=samples_boolean)
        
    if sites_boolean is not None and samples_boolean is not None:
        
        genotypes_subset = genotypes.subset(sel0=sites_boolean,sel1=samples_boolean)
                
    else:
        
        raise ValueError("Either a samples or a sites filter must be passed")
        
    allele_counts = allel.AlleleCountsArray(genotypes_subset.count_alleles())
    
    allele_counts_boolean = (allele_counts.max_allele() == max_alleles - 1) & (
        allele_counts[:, :2].min(axis=1) > min_count)
    
    num_alt_alleles = genotypes_subset.subset(allele_counts_boolean).to_n_alt()[:]
    
    return num_alt_alleles, allele_counts_boolean

def prune_by_LD(number_of_alternate_alleles,window_size=1000,step_size=100,r2=0.2):
    '''Take an array of the number of alternate alleles and return a smaller array pruned
    to remove SNPs in LD with each other above a specified threshold as well as the Boolean
    used for filtering.'''
    
    if not all(item > 0 for item in [window_size,step_size,r2]):
        raise ValueError("All numeric parameters must be positive")
        
    pruned_Bool = allel.locate_unlinked(number_of_alternate_alleles, window_size, step_size,
                                       r2)
    
    pruned = number_of_alternate_alleles[pruned_Bool]
    
    if not len(pruned) <  len(number_of_alternate_alleles):
        warnings.warn("Warning, no pruning occurred!")
        
    return pruned, pruned_Bool
