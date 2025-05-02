# Lisa's Code Repository

This repository contains all of Lisa's work on her master thesis project about Analysing Plant-Microbe Interactions via community modelling of metabolic interactions in the C7 Maize root SynCom.

There are 4 main folders:
(1) Datatsets, i.e. the protein fasta files (.faa) for every of the seven bacterial species that I work with and that are used to create the metabolic model with carveme \
(2) Models, i.e. metabolic models as .xml files after different stages of curation \
(3) Reports, i.e. memote and macaw reports to check model quality 
(4) Scripts, i.e. Jupyter Notebooks that contain my work in progress 


More explaination about my scripts: \
(01) create_metabolic_models.ipynb \
(02) balance_metabolic_models.ipynb (this will focuses on charge balancing; mass balance was done after creating the models with Frowins Scripts: Identify_imbalanced_reactions and apply_mass_balance_functions) \
(03) check_duplicates_and_loops.ipynb (more in-depth analysis about models and their structure with the goal to clean up topology and getting rid of energy-generating cycles) \

Additionally I used Scripts from Frowin (that can also be found in his Github repo) to mass balance my models and later on for fixes based on macaw testing. 
 


