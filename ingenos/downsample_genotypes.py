def downsample_genotypes(boolean_filter,downsample_rate=10):
    '''Downsample a set of genotypes for computational tractability.
    To be precise, this function stochastically shrinks the Boolean filter later used to pull out the genotypes.
    Shrinkage is approximate.'''
    
    downsampled_boolean = [random.randrange(100) < downsample_rate if value else value for value in boolean_filter]
    return downsampled_boolean
