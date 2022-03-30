from help_functions import install_R_packages
from repeated_simulation_largest import simulation

''' install the necessary R packages ''' 
# install_R_packages()

''' Define number of patients, n_readers and variables' ranges'''
n_patients = 5
n_readers = 5
reps = 5

''' Define default values'''
df_organs = 8 
df_lesions = 15

''' Run Simulation '''
discord_choice_of_lesions, discord_category, categs, percents = simulation(
    n_patients = n_patients, 
    n_readers = n_readers, 
    reps = reps, 
    df_organs = df_organs, 
    df_lesions = df_lesions,
    verb = True, 
    plot_disc = True
    )
