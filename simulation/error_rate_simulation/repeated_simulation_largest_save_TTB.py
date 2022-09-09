import numpy as np
import pandas as pd
import random
import math
from growth_simulation_fixed_lesions import dataset_simulation
from help_functions import pool_of_measurable_lesions, simulate_readers, plot_discordances
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina')

def simulation(n_patients = 100, n_readers = 100, reps = 100, 
               n_organs_range = range(1,11), n_lesions_range = range(1,21), mean_growth_range = range(-100,310,10), 
               var_pat_range = range(10,110,10), var_org_range = range(10,110,10), var_res_range = range(10,110,10),
               df_organs = 3, df_lesions = 10, df_mean_growth = 12.8, 
               df_var_pat = 44.55**2, df_var_org = 11.60**2, df_var_res = 36.62**2, 
               df_per = 20, which_var = [0,1,2,3,4,5], verb = True, plot_disc = True):

    ranges = [n_organs_range, n_lesions_range, mean_growth_range, var_pat_range, var_org_range, var_res_range]
    xaxixes = ['Max. No. organs', 'Max. No. lesions', 'Mean Growth', 'Patients Variance ', 'Organs Variance', ' Residuals Variance']
    
    discord_choice = []
    discord_categ = []
    categs = []
    percents = []
    ttbs_BL, ttbs_FU = [], []

    ''' Looping through each of the variables (organs, lesions, growth, and variances)'''
    for i, which_range in enumerate(ranges):

        ''' Fix every other variable that is not being changed '''        
        max_organs, max_lesions, growth = df_organs, df_lesions, df_mean_growth
        var_pat, var_org, var_res =  df_var_pat, df_var_org, df_var_res
        
        ''' Store the percentages of discordance '''
        percent_discordance_choice_of_lesions = np.zeros((reps,len(which_range)))
        percent_discordance_categories = np.zeros((reps,len(which_range)))
        
        if i not in which_var:
            continue
        else:
        
            print('\n', 'Varying', str(xaxixes[i]), 'within range', str(list(which_range)),'...')
    
            for k, var in enumerate(which_range):
        
                if i == 0:
                    max_organs = var
                elif i == 1:
                    max_lesions = var
                elif i == 2:
                    growth = var
                elif i == 3:
                    var_pat = var**2
                elif i == 4:
                    var_org = var**2
                elif i == 5:
                    var_res = var**2
    
                for rep in range(reps):
                    
                    ''' store individual percent growths | npatients X n_readers '''
                    percent = np.zeros((n_patients, n_readers))
                    ''' store individual categories | npatients X n_readers '''
                    categ = np.zeros((n_patients, n_readers))        
                    
                    ''' Simulate n_patients and n_readers'''        
                    for pat_ind in range(n_patients):
                        
                        ''' Create a lesions dataframe for this patient'''
                        
                        lesions = dataset_simulation(n_organs_max = max_organs, 
                                                            n_lesions_max = max_lesions, 
                                                            n_patients = 1, # we are simulating for just one patient in this step
                                                            mean_growth = growth, 
                                                            var_patients = var_pat, 
                                                            var_organs = var_org, 
                                                            var_residuals = var_res, 
                                                            seed = 42)
                        
                        ''' randomly define baseline size of each lesion '''
                        lesions['size_BL'] = [random.randint(10,100) for x in range(len(lesions))]
                        lesions["size_FU_for_TTB"] = lesions["size_BL"] * lesions["growth"]/100 + lesions["size_BL"]

                        ''' Define the lesions readers can choose from (eg: all lesions; only the largest by 20%)'''   
                        largest_lesions = pool_of_measurable_lesions(lesions, df_per)
                        ''' Simulate 'n_readers' and determine which categories each reader attributes to each patient''' 
                        percent, categ = simulate_readers(n_readers, largest_lesions, percent, categ, pat_ind)

                        ttbs_BL.append(lesions.sum()['size_BL'])
                        ttbs_FU.append(lesions.sum()['size_FU_for_TTB'])
                        
                    ''' Find those patients which have ≠ choice of lesions and  ≠ categories'''
                    diff_percent = pd.DataFrame(percent)[pd.DataFrame(percent).diff(axis=1).fillna(0).ne(0).any(axis=1)]
                    diff_categories = pd.DataFrame(categ)[pd.DataFrame(categ).diff(axis=1).fillna(0).ne(0).any(axis=1)]
                
                    if verb:
                        print('\n','Organs: ', max_organs,' | Lesions: ', max_lesions,' | growth: ',
                              growth,' | Var pat: ', math.sqrt(var_pat),' | Var org: ', 
                              math.sqrt(var_org),' | Var res: ', math.sqrt(var_res))
                        print('% of Discordance due to different choice of lesions', 
                              len(diff_percent)*100/n_patients)
                        print('% of Discordance due to different categories', 
                              len(diff_categories)*100/n_patients)
                    
                    ''' Store percentages to plot them later'''
                    percent_discordance_choice_of_lesions[rep,k] = len(diff_percent)*100/n_patients
                    percent_discordance_categories[rep,k] = len(diff_categories)*100/n_patients
                    categs.append(categ)
                    percents.append(percent)
            
        discord_choice.append(percent_discordance_choice_of_lesions)
        discord_categ.append(percent_discordance_categories)
            
        ''' Plot percentages against variable range '''
        if plot_disc:
            plot_discordances(discord_categ, which_range, xaxixes[i])

    return discord_choice, discord_categ, categs, percents, ttbs_BL, ttbs_FU