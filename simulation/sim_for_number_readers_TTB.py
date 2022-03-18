from help_functions import install_R_packages
from repeated_simulation_largest import simulation
import numpy as np
import pandas as pd
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina')

# =============================================================================
# CHANGE MAX NUMBER OF LESIONS AND ONLY VARIANT READERS CODE BEFORE RUNNING THIS
# =============================================================================


''' install the necessary R packages ''' 
install_R_packages()

''' Define number of patients, n_readers and variables' ranges'''
n_patients = 100
n_readers = 100
reps = 1

''' Define default values'''
df_organs = 3 # 8

n_lesions_range = range(3,21)

''' Run Simulation '''
discord_choice_of_lesions, discord_category, categs, percents, ttbs_BL, ttbs_FU = simulation(
    n_patients = n_patients, n_readers = n_readers, reps = reps, 
    df_organs = df_organs, n_lesions_range = n_lesions_range,  
    which_var = [1], verb = True, plot_disc = False, df_per = 20)

ttbs_percent = [(f-b)*100/b for f,b in zip(ttbs_FU, ttbs_BL)]

p = np.vstack(percents)

numb_lesions = list(n_lesions_range)*n_patients*reps
numb_lesions.sort()

final = np.zeros((len(numb_lesions),n_readers + 2))
final[:, 0] = numb_lesions
final[:, 1] = ttbs_percent
final[:, 2:] = p


df = pd.DataFrame(final)
# df.to_excel(excel_writer = "/home/teresa/stefanosidea/simulation/largest_var_readers.xlsx")
    