# Lisa's Code Repository

This repository contains all of Lisa's work on her master thesis project about Analysing Plant-Microbe Interactions via community modelling of metabolic interactions in the C7 Maize root SynCom.

There are 5 main folders: \
(1) Datatsets, i.e. the protein fasta files (.faa) for every of the seven bacterial species that I work with and that are used to create the metabolic model with carveme as well as csv files that contain metabolite and reaction info that I downloaded from BIGG. \
(2) Figures \
(3) Models, i.e. metabolic models as .xml files after different stages of curation \
(4) Reports, i.e. memote and macaw reports to check model quality 
(5) Scripts, i.e. Jupyter Notebooks that contain my work in progress 


More explaination about my scripts: \
(01) create_metabolic_models.ipynb \
(02) balance_metabolic_models.ipynb (mass and charge balancing) \
(03) check_duplicates.ipynb (more in-depth analysis about models and their structure with the goal to clean up topology) \
(04) simulations.ipynb \
(05) community.ipynb \
(06) dropouts.ipynb \
(07) vitaminB.ipynb


I used Scripts from Frowin (that can also be found in his Github repo) to mass balance my models and later on for fixes based on macaw testing. I fully integrated parts that I took from his scripts and integrated them into mine for an easier workflow on my end and credited the parts accordingly.
 


