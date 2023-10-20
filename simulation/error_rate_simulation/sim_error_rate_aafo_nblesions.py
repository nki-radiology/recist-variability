from repeated_simulation_largest_save_TTB import simulation
import numpy as np
import pandas as pd
from IPython.display import set_matplotlib_formats
import matplotlib.pyplot as plt
set_matplotlib_formats('retina')
import seaborn as sns
from matplotlib.pyplot import figure
from sklearn.metrics import accuracy_score
import warnings
warnings.filterwarnings("ignore")

''' Define number of patients, n_readers and variables' ranges'''
n_patients = 10
n_readers = 2
same_lesion_number_reps = 2
thr = None # if it is not None, then non-target disease is also assessed for PD, based on the threshold defined here
adj = 1 # whether to do adjudication or not. If yes, then two or three readers evaluate each patient

''' Define default values'''
df_per = 20
# immuno + chemo
df_var_pat = 38.0**2
df_var_org = 11.29**2
df_var_res = 34.45**2
df_mean_growth = 9.064
df_organs = 4

for df_organs in [2,4,8]:

file_name = str(n_patients) + '_n_readers_' + str(n_readers) + '_thr_' + str(thr) + '_df_organs_' + str(df_organs) + '_adj_' + str(adj)

n_lesions_range = range(3,21)

all_accs = np.zeros((len(list(n_lesions_range)), same_lesion_number_reps))

for rep in range(same_lesion_number_reps):
    
    print('#################### ' + str(rep))
    
    ''' Run Simulation '''
    discord_choice_of_lesions, discord_category, categs, percents, ttbs_BL, ttbs_FU = simulation(
        n_patients = n_patients, n_readers = n_readers, reps = 1,
        df_organs = df_organs, 
        n_lesions_range = n_lesions_range,  
        df_var_pat = df_var_pat,
        df_var_org = df_var_org,
        df_var_res = df_var_res,
        df_mean_growth = df_mean_growth,
        which_var = [1],
        thr = thr,
        adj = adj,
        verb = False, plot_disc = False, df_per = df_per)

    ttb_categs = []
    ttbs_percent = []
    
    ''' Determine the TTB category ('Ground Truth')'''
    for ttb_BL, ttb_FU in zip(ttbs_BL, ttbs_FU):
        
        abs_diff = abs(ttb_FU-ttb_BL)
        percent_growth = (ttb_FU-ttb_BL)*100/ttb_BL
        ttbs_percent.append(percent_growth)
    
        ''' Assign categories '''
        if (percent_growth > 20 and abs_diff > 5):
            ttb_categs.append(1) # progressive disease
        elif percent_growth == -100:
            ttb_categs.append(2) # complete response
        elif percent_growth < -30:
            ttb_categs.append(3) # partial response
        else:
            ttb_categs.append(4) # stable disease
    
    numb_lesions = list(n_lesions_range)*n_patients
    numb_lesions.sort()
    
    percen = np.zeros((len(numb_lesions), n_readers + 2))
    percen[:, 0] = numb_lesions
    percen[:, 1] = ttbs_percent
    percen[:, 2:] = np.vstack(percents)
    
    categ = np.zeros((len(numb_lesions), n_readers + 2))
    categ[:, 0] = numb_lesions
    categ[:, 1] = ttb_categs
    categ[:, 2:] = np.vstack(categs)
    
    percen = pd.DataFrame(percen)
    categ = pd.DataFrame(categ)

    categ = categ.loc[:, (categ != 0).any(axis=0)]        
    accs = []
    
    ''' iterate the lesions'''
    for l in n_lesions_range:
            
        ''' iterate the readers. Compute accuracy per reader '''
        data_categs = categ[categ.iloc[:,0] == l] # 0 is the position of the number of lesions
        ttb = data_categs.iloc[:,1] # 0 is the position of the number of ttb
        ttb = np.array(ttb)
        readings = data_categs.iloc[:,2:].T # 2 is where the readings start
        readings = readings.reset_index(drop = True)
        
        acc_per_reader = []
        
        for r, read in enumerate(readings.iterrows()): # each row is a reader
            y_pred = np.array(read[1]) # as long as the nb of patients
            acc_per_reader.append(accuracy_score(ttb, y_pred))
        
        accs.append(np.mean(acc_per_reader))
        
    nblesions = np.unique(numb_lesions)
    all_accs[:,rep] = accs

correct_estimations = pd.DataFrame(all_accs)
correct_estimations.insert(0,'nb_lesions', nblesions)
long_df = correct_estimations.iloc[0:len(correct_estimations),:2]
long_df.columns = ['nb_lesions', 'acc']

for k in range(2,correct_estimations.shape[1]-1):
    pp = pd.DataFrame(correct_estimations.iloc[0:len(correct_estimations),[0,k]])
    pp.columns = long_df.columns
    long_df = long_df.append(pp, ignore_index=True)

long_df['error'] = 1 - long_df['acc']

# Plot error rate as a function of the number of lesions
fsize = 20
figure(figsize=(8, 6))
g = sns.lineplot(data=long_df, 
                x="nb_lesions", 
                y="error", 
                errorbar = 'sd',
                markers = True,
                marker='o',
                #  palette=['r', 'tab:blue'],
                )

g.set_xlabel(r'$L$', fontsize=fsize+1, fontweight = 'bold')
g.set_ylabel('Error Rate', fontsize=fsize+1, fontweight = 'bold')
g.set_xticks([5, 10, 15, 20])
g.set_xticklabels([5, 10, 15, 20], fontsize = fsize)

s = g.get_yticks()
s = [round(a, 2) for a in s]
s = range(0, int(s[-1]*100), 4) #/100
s = [round(a/100, 2) for a in s]
g.set_yticks(list(s))
g.set_yticklabels(list(s), fontsize = fsize)
g.set_ylim([0,0.5])

long_df.to_csv(file_name + '.csv')
