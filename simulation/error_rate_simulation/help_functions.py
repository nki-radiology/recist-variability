import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector

def install_R_packages():
    
    ''' Install the necessary R packages'''
    package_names = ('stats', 'tmvtnorm', 'lme4', 'magic','statmod')
       
    if all(rpackages.isinstalled(x) for x in package_names):
        have_package = True
    
    else: 
        have_package = False
    
    if not have_package: 
        utils = rpackages.importr('utils')
        utils.chooseCRANmirror(ind=1)
        
        packnames_to_install = [x for x in package_names if not rpackages.isinstalled(x)]
        
        if len(packnames_to_install) > 0: 
            utils.install_packages(StrVector(packnames_to_install))


def pool_of_measurable_lesions(lesions, per):
    
    ''' Create a pool of "lesions" (all lesions that have a size <=per% the size of the biggest lesions/organ)'''
    ''' order by organ and by size '''
    ordered_lesions = lesions.sort_values(by=['organ', 'size_BL'], ascending = (True, False))
    
    ''' identify the largest lesion of all organs '''
    largest_lesion_per_organ = ordered_lesions.groupby('organ').head(1).reset_index(drop=True)
    
    ''' select largest lesions per organ (<=20% of the largest lesion of that organ)'''
    largest_lesions = pd.DataFrame(columns = largest_lesion_per_organ.columns)

    for this_organ in np.unique(largest_lesion_per_organ.organ):
        diam_largest_lesion = largest_lesion_per_organ[largest_lesion_per_organ.organ == this_organ].reset_index(drop=True).at[0,'size_BL']
        all_lesions_from_this_organ = ordered_lesions[ordered_lesions.organ == this_organ].reset_index(drop=True)
        largest_lesions_from_this_organ = all_lesions_from_this_organ[((all_lesions_from_this_organ.size_BL-diam_largest_lesion)*-100)/diam_largest_lesion <= per].reset_index(drop=True)

        ''' Even if it is not within 20% of the size of the biggest lesion, we still want two max/organ, so let's choose the second biggest as well (and all same sized lesions (<=20%))'''
        if len(all_lesions_from_this_organ) > 1 and len(largest_lesions_from_this_organ) < 2:
            all_but_the_biggest = all_lesions_from_this_organ.tail(len(all_lesions_from_this_organ)-1).reset_index(drop=True)
            diam_2nd_largest_lesion = all_but_the_biggest.at[0,'size_BL']
            second_largest_lesions_from_this_organ = all_but_the_biggest[((all_but_the_biggest.size_BL-diam_2nd_largest_lesion)*-100)/diam_2nd_largest_lesion <= per]
            frames = [largest_lesions_from_this_organ, second_largest_lesions_from_this_organ]
            largest_lesions_from_this_organ = pd.concat(frames)
        
        frames = [largest_lesions, largest_lesions_from_this_organ]
        largest_lesions = pd.concat(frames)
    
    return largest_lesions



def simulate_readers(n_readers, lesions, percents, categ, pat_ind):

    for reader in range(n_readers):
        ''' make sure we don't choose more than 2 lesions/organ. Sample shuffles the organs before sampling'''
        remove_repeated_organs = lesions.sample(frac=1).groupby('organ').tail(2) 
        
        ''' make sure we don't choose more than 5 lesions in total '''
        N = len(remove_repeated_organs) if len(remove_repeated_organs)<5 else 5 
        
        ''' Choose 5 random lesions. Sample shuffles the lesions before sampling '''
        final_5_lesions = remove_repeated_organs.sample(n = N).reset_index(drop=True)
        
        ''' Calculate percent growth'''
        size_FU = (final_5_lesions.size_BL*final_5_lesions.growth)/100 + final_5_lesions.size_BL
        final_5_lesions['size_FU'] = size_FU
        final_5_lesions_sum = final_5_lesions.sum()
        abs_diff = final_5_lesions_sum.size_FU - final_5_lesions_sum.size_BL
        percent_growth = (final_5_lesions_sum.size_FU - final_5_lesions_sum.size_BL)*100/final_5_lesions_sum.size_BL
        percents[pat_ind,reader] = percent_growth
        
        ''' Assign categories '''
        if (percent_growth > 20 and abs_diff > 5):
            categ[pat_ind,reader] = 1 # progressive disease
        elif percent_growth == -100:
            categ[pat_ind,reader] = 2 # complete response
        elif percent_growth < -30:
            categ[pat_ind,reader] = 3 # partial response
        else:
            categ[pat_ind,reader] = 4 # stable disease
            
    return percents, categ


def get_non_target(all_lesions_growth, target_growth):
    
    non_target_index = []

    for k, lesion_growth in enumerate(all_lesions_growth):
        
        if lesion_growth not in target_growth:
            non_target_index.append(k)

    return non_target_index



def simulate_readers_nontarget(n_readers, lesions, percents, categ, pat_ind, thr):

    #################################
    # RECIST 1.0
    # lesions_total = 10
    # lesions_organ = 5
    #################################
    # RECIST 1.1
    lesions_total = 5
    lesions_organ = 2


    for reader in range(n_readers):
        ''' make sure we don't choose more than 2 lesions/organ. Sample shuffles the organs before sampling'''
        remove_repeated_organs = lesions.sample(frac=1).groupby('organ').tail(lesions_organ) 
        
        ''' make sure we don't choose more than 5 lesions in total '''
        N = len(remove_repeated_organs) if len(remove_repeated_organs)<lesions_total else lesions_total 


        ''' Choose 5 random lesions. Sample shuffles the lesions before sampling '''
        final_5_lesions = remove_repeated_organs.sample(n = N).reset_index(drop=True)
        
        ''' Calculate percent growth'''
        size_FU = (final_5_lesions.size_BL*final_5_lesions.growth)/100 + final_5_lesions.size_BL
        final_5_lesions['size_FU'] = size_FU

        final_5_lesions_sum = final_5_lesions.sum()
        abs_diff = final_5_lesions_sum.size_FU - final_5_lesions_sum.size_BL
        percent_growth = (final_5_lesions_sum.size_FU - final_5_lesions_sum.size_BL)*100/final_5_lesions_sum.size_BL
        percents[pat_ind,reader] = percent_growth
        
        ''' Assign categories '''
        if (percent_growth > 20 and abs_diff > 5):
            categ[pat_ind,reader] = 1 # progressive disease
        elif percent_growth == -100:
            categ[pat_ind,reader] = 2 # complete response
        elif percent_growth < -30:
            categ[pat_ind,reader] = 3 # partial response
        else:
            categ[pat_ind,reader] = 4 # stable disease

        ''' analyze non-target to see if we should change the category '''

        # reset index to select non-target
        lesions = lesions.reset_index(drop=True)
        # get non target index by matching the growth
        nontarg_ind = get_non_target(list(lesions.growth), list(final_5_lesions.growth))
        non_target_lesions = lesions.iloc[nontarg_ind]

        if len(non_target_lesions) != 0:
            ''' Calculate percent growth'''
            size_FU = (non_target_lesions.size_BL*non_target_lesions.growth)/100 + non_target_lesions.size_BL
            non_target_lesions['size_FU'] = size_FU
            non_target_lesions = non_target_lesions.sum()
            # abs_diff_nt = non_target_lesions.size_FU - non_target_lesions.size_BL
            percent_growth_nt = (non_target_lesions.size_FU - non_target_lesions.size_BL)*100/non_target_lesions.size_BL
            # percents[pat_ind,reader] = percent_growth
            if percent_growth_nt > thr:
                categ[pat_ind,reader] = 1 # progressive disease

    return percents, categ


def reader_loop(lesions, thr = None):
    
    #################################
    # RECIST 1.0
    # lesions_total = 10
    # lesions_organ = 5
    #################################
    # RECIST 1.1
    lesions_total = 5
    lesions_organ = 2
    
    ''' make sure we don't choose more than 2 lesions/organ. Sample shuffles the organs before sampling'''
    remove_repeated_organs = lesions.sample(frac=1).groupby('organ').tail(lesions_organ) 
    
    ''' make sure we don't choose more than 5 lesions in total '''
    N = len(remove_repeated_organs) if len(remove_repeated_organs)<lesions_total else lesions_total 

    ''' Choose 5 random lesions. Sample shuffles the lesions before sampling '''
    final_5_lesions = remove_repeated_organs.sample(n = N).reset_index(drop=True)
    
    ''' Calculate percent growth'''
    size_FU = (final_5_lesions.size_BL*final_5_lesions.growth)/100 + final_5_lesions.size_BL
    final_5_lesions['size_FU'] = size_FU

    final_5_lesions_sum = final_5_lesions.sum()
    abs_diff = final_5_lesions_sum.size_FU - final_5_lesions_sum.size_BL
    percent_growth = (final_5_lesions_sum.size_FU - final_5_lesions_sum.size_BL)*100/final_5_lesions_sum.size_BL
    p = percent_growth
    
    ''' Assign categories '''
    if (percent_growth > 20 and abs_diff > 5):
        c = 1 # progressive disease
    elif percent_growth == -100:
        c = 2 # complete response
    elif percent_growth < -30:
        c = 3 # partial response
    else:
        c = 4 # stable disease

    ''' analyze non-target to see if we should change the category '''
    if thr is not None:
        # reset index to select non-target
        lesions = lesions.reset_index(drop=True)
        # get non target index by matching the growth
        nontarg_ind = get_non_target(list(lesions.growth), list(final_5_lesions.growth))
        non_target_lesions = lesions.iloc[nontarg_ind]

        if len(non_target_lesions) != 0:
            ''' Calculate percent growth'''
            size_FU = (non_target_lesions.size_BL*non_target_lesions.growth)/100 + non_target_lesions.size_BL
            non_target_lesions['size_FU'] = size_FU
            non_target_lesions = non_target_lesions.sum()
            # abs_diff_nt = non_target_lesions.size_FU - non_target_lesions.size_BL
            percent_growth_nt = (non_target_lesions.size_FU - non_target_lesions.size_BL)*100/non_target_lesions.size_BL
            # percents[pat_ind,reader] = percent_growth
            if percent_growth_nt > thr:
                c = 1 # progressive disease

    return p, c


def simulate_readers_nontarget_adjudicator(lesions, percents, categ, pat_ind, thr):

    classes = []
    # This is done for 2 readers and potentially and adjudicator
    for reader in range(2):
        p,c = reader_loop(lesions, thr)
        classes.append(c)

    if classes[0]!=classes[1]: # get a third reader as an adjudicator
        p,c = reader_loop(lesions, thr)
        print('readers',classes)
        print('adj',c)


    percents[pat_ind,reader] = p
    categ[pat_ind,reader] = c
    if classes[0]!=classes[1]: 
        print('fianl', c)

    return percents, categ

def plot_discordances(discordance, xrange, xlabel):
    
    res = pd.DataFrame(discordance[-1])
    ranges = list(xrange)
    
    lesion_rows = []
    lesion_index = []
    for k in range(0, len(res)):
        lesion_rows.append(res.iloc[k].values)
        lesion_index.append(ranges)
    
    lesion_rows = np.array([item for sublist in lesion_rows for item in sublist])
    lesion_index = np.array([item for sublist in lesion_index for item in sublist])
    lesions_long_df = pd.DataFrame([lesion_rows, lesion_index]).T
    lesions_long_df.columns = ['Discordance (%)', xlabel]
    plt.figure()
    sns.lineplot(data=lesions_long_df, x=xlabel, y="Discordance (%)", errorbar = 'sd')
    plt.savefig(xlabel + '.png')