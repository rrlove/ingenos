import allel
import warnings
import h5py
import pandas as pd

def import_data(filepath='/afs/crc.nd.edu/group/BesanskyNGS/data05/comp_karyo/data/ag1000g.phase1.ar3.pass.2R.h5', chrom_name='2R'):
    '''Take the path to a well-formed h5py file and return a VariantTable and a GenotypeArray.'''
    
    ##to-do: check that h5py file is well-formed
    
    callset_handle = filepath
    callset = h5py.File(callset_handle, mode='r')
    
    variants = allel.VariantChunkedTable(callset[chrom_name]['variants'],names=['POS','REF','ALT','DP','MQ','QD','num_alleles',], index='POS')
    
    genotypes = allel.GenotypeChunkedArray(callset[chrom_name]['calldata']['genotype'])
    
    if not genotypes.shape[0] == variants.shape[0]:
        raise ValueError("Genotypes and variant table must contain the same number of positions")

    return variants, genotypes

def import_metadata(md_filepath='/afs/crc.nd.edu/group/BesanskyNGS/data05/comp_karyo/metadata/Ag1K_phase1_metadata.txt',
                    kt_filepath='/afs/crc.nd.edu/group/BesanskyNGS/data05/comp_karyo/metadata/Cameroon_phase1_karyotypes.xlsx',
                    training_set_filepath='/afs/crc.nd.edu/group/BesanskyNGS/data05/comp_karyo/metadata/Cameroon_phase1_karyotyped_training_set.txt'):
    
    '''Import the AG1K metadata, the file with karyotypes given, and the list of samples in
    the training set given the correct paths.'''
    
    metadata = pd.DataFrame.from_csv(md_filepath, sep='\t', dtype={'year':int})
    
    karyotypes = pd.read_excel(kt_filepath, skiprows=6, na_values=['#N/R'],
                           converters={'2Rj':int,'2Rb':int,'2Rc':int,'2Rd':int,'2Ru':int,'2La':int})
    
    training = pd.read_table(training_set_filepath, header=None)
    
    return metadata, karyotypes, training

def process_metadata(metadata, karyotypes, training):
    
    '''Return:
    1.) Karyotypes for the samples in the training set that are also in the metadata
    (samples not in the metadata were not sequenced).
    2.) A boolean the length of the metadata (and therefore the genotype array)
    indicating which samples are in the training set.'''
    
    sequenced_training_set_IDs = training[training[0].isin(metadata["src_code"])]
    
    karyotyped_training_samples = karyotypes[karyotypes["Sample ID"].isin(sequenced_training_set_IDs[0])]
    
    pull_out_training_set = metadata["src_code"].isin(training[0])
    
    return karyotyped_training_samples, pull_out_training_set

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
