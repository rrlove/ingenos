
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

