def construct_filter_bool(variants,expression):
    '''Take a set of variants and a positional expression on which to filter them, return a Boolean as to whether each variant/row of genotypes should be included post-filtering.'''

    selection_bool = variants.eval(expression)[:]
    return selection_bool
