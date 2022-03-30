# Studying RECIST Variability with Computer Simulated Models 💻📈

Managed by Teresa Tareco-Bucho (t.tareco.bucho@nki.nl) and Stefano Trebeschi (s.trebeschi@nki.nl).

RECIST is grounded on the assumption that target lesion selection is objective and representative of total tumor burden’s response to therapy. A computer simulation model was designed to investigate the impact of target lesion selection on response assessment. Readers’ disagreement as a function of the total number of lesions, affected organs and lesion growth was analyzed. Disagreement aggravates when the number of lesions increases, when lesions are concentrated on few organs and when lesion growth borders the thresholds of progressive disease and partial response. In a metastatic setting, RECIST displays a highly non-linear, unpredicatable behaviour, not reliable for response assessment due to target lesion selection variability. Including more (if not all) lesions in the quantitative analysis of tumor burden is desirable.

In this repo, we share the code of the simulation model, and the statistical analysis performed in the publication. 
You can find the full publication pre-print at the following link: `TO ADD``

## 1. Requirements

1. Create and activate a conda environment with Python
``conda create -n recist_sim python``

``conda activate recist_sim``

2. Install R, numpy, pandas, seaborn, rpy2 (linux only) and r-packages tmvtnorm and lme4
``conda install -c r r``

``conda install numpy pandas seaborn``

``pip install rpy2``

``conda install -c conda-forge r-tmvtnorm``

``conda install -c conda-forge r-lme4``

## 2. How to run the Simulation Model
In the ``run_simulation.py`` file, specify the inputs to the simulation model. Particularly, you can specify the number of readers, patients and repetitions, the ranges for Lmax, Omax, miu, and variances, among others. Check the arguments of the ``simulation`` function. 

``python run_simulation.py``

If ``plot_disc = True``, the discordance plots will be save in your working directory.

If you encounter a ``libRlapack.so``/``libRblas.so`` error, a workaround is: 

``cd /path_to_your_env/lib/``

``mv liblapack.so libRlapack.so``

``mv libblas.so libRblas.so``

## 3. Contribution

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are greatly appreciated. If you want to add your analysis, or have a suggestion that would make this better, please fork the 