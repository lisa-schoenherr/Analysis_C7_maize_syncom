# Code Repository - Master Thesis - Lisa Schönherr

This repository contains all of Lisa's work on her master thesis project about Analysing Plant-Microbe Interactions via community modelling of metabolic interactions in the C7 Maize root SynCom.

There are 5 main folders: 
1. Datatsets 
   * protein fasta files (.faa) for the seven bacterial species; used to create the metabolic model with carveme
   * csv files for creating a community model with MICOM
   * csv files that contain metabolite and reaction info that I downloaded from BIGG
   * csv files that have different media (mainly maize root exudates)
2. Figures 
3. Models, i.e. metabolic models as .xml files after different stages of curation 
4. Reports, i.e. MEMOTE and MACAW reports to check model quality 
5. Scripts, i.e. Jupyter Notebooks that contain my work in progress as well as a python file that contains all functions that are used in multiple notebooks


More explaination about my scripts: 
1. create_metabolic_models.ipynb 
   * CarveMe is used to create metabolic models
   * MEMOTE reports
2. balance_metabolic_models.ipynb
   * mass and charge balancing
3. check_duplicates.ipynb
   * MACAW; finding and cleaning up duplicates
   * investigation about dead-ends and other topological issues
4. simulations.ipynb 
   * simulations with the individual models: growth assays, metabolic niches
   * gap-filling of models
5. community.ipynb 
   * create community model with MICOM (C7)
   * simulations with this model: growth & niches, cross-feeding
6. dropouts.ipynb 
   * creating drop-out communities (C6)
   * simulations with these models: growth assays, niches, PCA
7. vitaminB.ipynb
   * adding B vitamins to the medium: V8 and V1/V7 experiments
   * gap-filling vitamin synthesis pathways
   * adding possible vitamin auxotrophies



I used Scripts from Frowin Reichhardt (AG Töpfer) to mass balance my models and later on for fixes based on macaw testing. I fully integrated parts that I took from his scripts and integrated them into mine for an easier workflow on my end and credited the parts accordingly. I used Scripts from Tiago Machado (AG Töpfer) to create Budget Plots for visualisation.
 


