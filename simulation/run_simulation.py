from help_functions import install_R_packages
from repeated_simulation_largest import simulation
import pandas as pd
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina')

''' install the necessary R packages ''' 
install_R_packages()

''' Define number of patients, n_readers and variables' ranges'''
n_patients = 1
n_readers = 1
reps = 1

''' Define default values'''
df_organs = 8 
df_lesions = 15
which_var = [0,1]

''' Run Simulation '''
discord_choice_of_lesions, discord_category, categs, percents = simulation(
    n_patients = n_patients, 
    n_readers = n_readers, 
    reps = reps, 
    df_organs = df_organs, 
    df_lesions = df_lesions,
    verb = True, 
    plot_disc = True,
    which_var = which_var,
    )
