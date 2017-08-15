def prune_by_LD(number_of_alternate_alleles,window_size=1000,step_size=100,r2=0.2):
    '''Take an array of the number of alternate alleles and return a smaller array pruned to remove SNPs in LD with each other.'''

    pruned_Bool = allel.locate_unlinked(number_of_alternate_alleles, window_size, step_size, r2)
    pruned = number_of_alternate_alleles[pruned_Bool]
    
    return pruned
