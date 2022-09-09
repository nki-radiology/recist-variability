import numpy as np
import itertools
import pandas as pd
import random
import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri

import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina')

def rpy2_to_pandas(r_df):
    ''' 
    To convert R variables to python
    '''
    with localconverter(robjects.default_converter + pandas2ri.converter):
        return robjects.conversion.rpy2py(r_df)

def dataset_simulation(n_organs_max, 
                       n_lesions_max, 
                       n_patients, 
                       mean_growth, 
                       var_patients, 
                       var_organs, 
                       var_residuals, 
                       seed):
    
    data_simulated = pd.DataFrame()
    
    for patient in range(n_patients):
        n_lesions = n_lesions_max # Different from the previous simulation
        n_organs = random.randint(1,n_organs_max)

        # R. Numerical variables are converted to strings
        robjects.r('''
            library(tmvtnorm)
            library(lme4)
            library(nlme)
            library(magic)
            n_lesions <- '''  + str(n_lesions) + '''
            n_organs <- '''  + str(n_organs) + '''
            mean_growth <- '''  + str(mean_growth) + '''
            var_patients <- '''  + str(var_patients) + '''
            var_organs <- '''  + str(var_organs) + '''
            var_residuals <- '''  + str(var_residuals) + '''
            n_lesions_organ <- rmultinom(1,n_lesions, prob = rep(1/n_organs, n_organs))
            matrices <- lapply(n_lesions_organ, FUN = function(x){matrix(1,nrow = x,  ncol = x)})
            matrix_organs <- Reduce(adiag,matrices)
            matrix_patient <- matrix(1,ncol = n_lesions, nrow = n_lesions)
            matrix_residuals <- diag(1,n_lesions)
            matrix <- matrix_organs * var_organs + matrix_patient * var_patients + matrix_residuals  * var_residuals           
            growth_patient <- rtmvnorm(1, mean = rep(mean_growth, n_lesions), sigma = matrix, lower = rep(-100,n_lesions))''')

        # Access R code variables
        rkernel = rpy2_to_pandas(robjects.r.get('growth_patient'))
        n_lesion_organ = np.squeeze(rpy2_to_pandas(robjects.r.get('n_lesions_organ')))
        if n_organs == 1:
            n_lesion_organ = [n_lesion_organ]
        
        # Create patient dataframe
        patient_id = np.array([patient]*n_lesions)
        organ_ids = list(range(len(n_lesion_organ)))
        organ = list(itertools.chain(*(itertools.repeat(elem, n) for elem, n in zip(organ_ids, n_lesion_organ))))
        A = np.array([patient_id, organ, rkernel[0]])
        results_patient = pd.DataFrame(A).T
        results_patient.columns = ['patient_id', 'organ', 'growth']
        
        # Concatenate patient dataframe to final dataframe
        data_simulated = pd.concat([data_simulated,results_patient])
        
    return data_simulated