'''This module contains the ingenos methods.'''

import datetime
import warnings

import allel
import h5py
import numpy as np
import pandas as pd

def compute_concordance(alleles, is_called, karyos):
    '''For each SNP, calculate the % of samples in which
    karyotype and genotype are concordant.'''

    assert len(alleles) == len(is_called) == len(karyos), \
    "Differing lengths of genotypes, is_called, karyotypes"

    assert np.min(karyos) >= 0, \
    "Improper karyotype {karyo}".format(karyo=np.min(karyos))

    assert np.max(karyos) <= 2, \
    "Improper karyotype {karyo}".format(karyo=np.max(karyos))

    inv_karyos = abs(2 - karyos)

    assert np.min(inv_karyos) >= 0, \
    "Improper karyotype {karyo}".format(karyo=np.min(karyos))

    assert np.max(inv_karyos) <= 2, \
    "Improper karyotype {karyo}".format(karyo=np.max(karyos))

    score = np.sum(alleles[is_called] ==\
                   karyos[is_called])/len(karyos[is_called])

    inv_score = np.sum(alleles[is_called] ==\
                       inv_karyos[is_called])/len(karyos[is_called])

    return np.max([score, inv_score])

#for now, I'm dropping the calculation of inverse karyotypes
def compute_concordance_strat(alleles, is_called, karyos):
    '''For each SNP, calculate the % of samples in each karyotype
    in which karyotype and genotype are concordant.'''

    assert len(alleles) == len(is_called) == len(karyos),\
    "Differing lengths of genotypes, is_called, karyotypes"

    assert np.min(karyos) >= 0, \
    "Improper karyotype {karyo}".format(karyo=np.min(karyos))

    assert np.max(karyos) <= 2, \
    "Improper karyotype {karyo}".format(karyo=np.max(karyos))
    
    assert np.min(alleles) >= 0, \
    "Improper alt allele count {count}".format(count = np.min(alleles))
    
    assert np.max(alleles) <= 2, \
    "Improper alt allele count {count}".format(count = np.max(allelesmaaaaa))

    #homo standard
    bool_0 = (karyos == 0)

    score_0 = np.sum(alleles[bool_0][is_called[bool_0]] == 0) /\
    len(karyos[bool_0][is_called[bool_0]])

    called_0 = np.sum(is_called[bool_0]) / np.sum(bool_0)

    #hets
    bool_1 = (karyos == 1)

    score_1 = np.sum(alleles[bool_1][is_called[bool_1]] == 1) /\
    len(karyos[bool_1][is_called[bool_1]])

    called_1 = np.sum(is_called[bool_1]) / np.sum(bool_1)

    #homo inv
    bool_2 = (karyos == 2)

    score_2 = np.sum(alleles[bool_2][is_called[bool_2]] == 2) /\
    len(karyos[bool_2][is_called[bool_2]])

    called_2 = np.sum(is_called[bool_2]) / np.sum(bool_2)

    #all
    score = np.sum(alleles[is_called] ==\
                   karyos[is_called])/len(karyos[is_called])

    called = np.sum(is_called)/len(karyos)

    return(score_0, score_1, score_2, score,
           called_0, called_1, called_2, called)

def construct_filter_expression(name, inversion_dict,
                                buffer=1500000, whole_inversion=False):
    '''Construct an expression describing the desired SNPs.
    Buffer describes the # of bp around each breakpoint to include.
    Whole_inversion specifies whether only the breakpoints and the buffer,
    or the actual inversion, should be used.'''

    proximal_start = inversion_dict[name].proximal_start - buffer
    proximal_end = inversion_dict[name].proximal_end + buffer
    distal_start = inversion_dict[name].distal_start - buffer
    distal_end = inversion_dict[name].distal_end + buffer

    if not all(item > 0 for item in [proximal_start, proximal_end,
                                     distal_start, distal_end]):
        raise ValueError("An output coordinate is negative,\
                         please adjust buffer")

    if whole_inversion:

        expression = '( (POS > {prox_start}) & (POS < {dist_end}) )'.format(
            prox_start=proximal_start,
            dist_end=distal_end)

    else:

        expression = '( ( (POS > {prox_start}) & (POS < {prox_end}) ) |\
        ( (POS > {dist_start}) & (POS < {dist_end}) ) )'.format(
            prox_start=proximal_start,
            prox_end=proximal_end,
            dist_start=distal_start,
            dist_end=distal_end)

    return expression

def filter_and_convert_genotypes(genotypes, sites_boolean=None,
                                 samples_boolean=None, max_alleles=2,
                                 min_count=3, variance_threshold=0.15):
    '''Filter a genotype array based on booleans of sites
    and samples to include.

    Further filter genotypes based on allele count data.

    Return a set of alternate allele counts ready for PCA,
    and the allele counts filter.
    '''

    if not all(item > 0 for item in [max_alleles]):
        raise ValueError("Max alleles must be greater than 0 ")

    if sites_boolean is not None:

        if not len(sites_boolean) == genotypes.shape[0]:
            raise ValueError("Length of sites filter \
                             does not match length of genotypes")

    if samples_boolean is not None:

        if not len(samples_boolean) == genotypes.shape[1]:
            raise ValueError("Length of samples filter \
                             does not match length of genotypes")

    if sites_boolean is not None and samples_boolean is None:

        genotypes_subset = genotypes.subset(sel0=sites_boolean)

    elif samples_boolean is not None and sites_boolean is None:

        genotypes_subset = genotypes.subset(sel1=samples_boolean)

    elif sites_boolean is not None and samples_boolean is not None:

        genotypes_subset = genotypes.subset(sel0=sites_boolean,
                                            sel1=samples_boolean)

    else:

        raise ValueError("Either a samples or a sites filter must be passed")

    allele_counts = allel.AlleleCountsArray(genotypes_subset.count_alleles())

    allele_counts_boolean = (allele_counts.max_allele() <= max_alleles - 1) &\
    (allele_counts[:, :2].min(axis=1) > min_count) & (
        allele_counts.to_frequencies()[:,1] > variance_threshold)

    num_alt_alleles = genotypes_subset.subset(
        allele_counts_boolean).to_n_alt()[:]

    return num_alt_alleles, allele_counts_boolean

def import_data(filepath, chrom_name,
                names=['POS', 'REF', 'ALT', 'DP', 'MQ', 'QD', 'num_alleles']):
    '''Take the path to a well-formed h5py file and return a VariantTable
    and a GenotypeArray.'''

    ##to-do: check that h5py file is well-formed

    callset_handle = filepath
    callset = h5py.File(callset_handle, mode='r')

    variants = allel.VariantChunkedTable(callset[chrom_name]['variants'],
                                         names=names, index='POS')

    genotypes = allel.GenotypeChunkedArray(
        callset[chrom_name]['calldata']['genotype'])

    if not genotypes.shape[0] == variants.shape[0]:
        raise ValueError("Genotypes and variant table must contain the\
                         same number of positions")

    return variants, genotypes

def import_metadata(md_filepath=
                    '/afs/crc.nd.edu/group/BesanskyNGS/data05/comp_karyo/\
                    metadata/Ag1K_phase1_metadata.txt',
                    kt_filepath=
                    '/afs/crc.nd.edu/group/BesanskyNGS/data05/comp_karyo/\
                    metadata/Cameroon_phase1_karyotypes.xlsx',
                    training_set_filepath=
                    '/afs/crc.nd.edu/group/BesanskyNGS/data05/comp_karyo/\
                    metadata/Cameroon_phase1_karyotyped_training_set.txt'):

    '''Import the AG1K metadata, the file with karyotypes given,
    and the list of samples in the training set, given the correct paths.'''

    metadata = pd.read_csv(md_filepath, sep='\t', dtype={'year':int})

    karyotypes = pd.read_excel(kt_filepath, skiprows=6, na_values=['#N/R'],
                               converters={'2Rj':int,
                                           '2Rb':int,
                                           '2Rc':int,
                                           '2Rd':int,
                                           '2Ru':int,
                                           '2La':int})

    training = pd.read_table(training_set_filepath, header=None)

    return metadata, karyotypes, training

def make_column_tuples(metadata, should_print=True):
    '''From a pandas DataFrame metadata object, make a multi-index for the
    columns of a DataFrame used to look at by-SNP concordance.'''

    new_tuples = []

    for sample in metadata["ox_code"]:

        country = metadata.loc[metadata["ox_code"] == sample, \
                               "country"].values[0]

        species = metadata.loc[metadata["ox_code"] == sample, \
                               "species"].values[0]

        id_string = '_'.join([country, species.lstrip('An. ')])

        idx = metadata.loc[metadata["ox_code"] == sample].index[0]

        if should_print:

            print(species, sample, id_string, idx)

        new_tuples.append((sample, id_string, idx))

    return new_tuples

def make_date_stamp(filename, extension=".pdf"):
    '''Generate filename paths with the current date.'''

    date_stamp = filename + '.' + datetime.date.today().strftime("%m%d%y") + \
    extension

    return date_stamp

def make_snp_index(positions, vtbl):
    '''Make a multi-index for the rows of a DataFrame used to look at by-SNP
    concordance. The resulting multi-index gives the position of each SNP
    in genome coordinates as well as the index of that SNP in the genotype
    object.'''

    indices = []

    for position in positions:

        where = np.where(vtbl["POS"] == position)

        if len(where[0]) > 0:
            indices.append(where[0][0])

    return indices

def plot_pca_coords(coords, model, pc1, pc2, ax, metadata, inversion_color,
                    region=None, alpha=0.5, title=None, buffer=False):
    ##adapted from Alistair's walkthrough
    '''Show two principal components of a PCA.'''

    marker_list = ['o', 'v', 'X', 's', 'p', 'd', '^', '_', '*', 'H']

    x = coords[:, pc1]
    y = coords[:, pc2]
    for index, status in enumerate(metadata[inversion_color].unique()):
        flt = (metadata[inversion_color] == status).values
        ax.plot(x[flt], y[flt], marker=marker_list[index], linestyle=' ',
                label=status, markersize=6, alpha=alpha)

    ax.set_xlabel('PC%s (%.1f%%)' % \
                  (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % \
                  (pc2+1, model.explained_variance_ratio_[pc2]*100))

    if title is None and buffer is True:
        ax.set_title('PC%s vs. PC%s for %s, colored by %s' % \
                     (pc1+1, pc2+1, region, inversion_color))
    elif title is None and buffer is False:
        ax.set_title('PC%s vs. PC%s for %s, colored by %s, no buffer' % \
                     (pc1+1, pc2+1, region, inversion_color))
    else:
        ax.set_title(title)

def process_metadata(metadata, karyotypes, training_set, inversion,
                     md_label="src_code", kt_label="Sample ID",
                     tr_label="ox_code"):

    '''Return:
    1.) Karyotyped samples in both the metadata and the training set
    (samples not in the metadata were not sequenced)
    2.) A boolean the length of the metadata
    (and therefore also of the genotype array)
    indicating which samples are in 1.) above
    '''

    kt_md = pd.merge(metadata, karyotypes, left_on=md_label, right_on=kt_label)

    kt_md_tr = kt_md[~kt_md[inversion].isnull()]\
    [kt_md[~kt_md[inversion].isnull()][md_label].isin(training_set[0])]

    filter_genotypes = metadata[tr_label].isin(kt_md_tr[tr_label])

    return kt_md_tr, filter_genotypes

def prune_by_ld(number_of_alternate_alleles, window_size=1000,
                step_size=100, r2=0.2):
    '''Take an array of the number of alternate alleles and return a smaller
    array pruned to remove SNPs in LD with each other above a specified
    threshold as well as the Boolean used for filtering.'''

    if not all(item > 0 for item in [window_size, step_size, r2]):
        raise ValueError("All numeric parameters must be positive")

    pruned_bool = allel.locate_unlinked(number_of_alternate_alleles,
                                        window_size, step_size, r2)

    pruned = number_of_alternate_alleles[pruned_bool]

    if not len(pruned) < len(number_of_alternate_alleles):
        warnings.warn("Warning, no pruning occurred!")

    return pruned, pruned_bool
