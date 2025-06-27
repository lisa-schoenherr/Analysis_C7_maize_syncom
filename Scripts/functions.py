"""
This python file contains functions that I use within my Jupyter Notebooks.
Especially, this python file focuses on functions that are used by multiple notebooks and
where it is easier to import them than to multiple copy and paste them into the notebooks.
Every Jupyter Notebooks still contains functions, but these are then more unique to the notebook.
"""

###
# PART 0 - Imports
###

import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.exceptions import Infeasible


###
# PART 1 - Get Metabolite & Reaction Info
###

def get_rxn(model, rxn_id, bounds = False, mass = False, GPR=False, charge=True):
    rxn = model.reactions.get_by_id(rxn_id)
    charges = {met.id: met.charge for met in rxn.metabolites}
    masses = {met.id: met.formula for met in rxn.metabolites}
    gpr = rxn.gene_reaction_rule
    up = rxn.upper_bound
    lb = rxn.lower_bound
    print(rxn)

    if charge:
        print(charges)
    if mass:
        print(masses)
    if GPR:
        print(gpr)
    if bounds:
        print(f"from {lb} to {up}")


def get_rxn_unknown(models, rxn_id, bounds=False, mass = False, GPR=False, charge=True):
    model_list = []
    gpr_list= []
    ub = set()
    lb = set()
    for model in models.values():
        if rxn_id in model.reactions:
            model_list.append(model)
            gpr_list.append(model.reactions.get_by_id(rxn_id).gene_reaction_rule)
            lb.add(model.reactions.get_by_id(rxn_id).lower_bound)
            ub.add(model.reactions.get_by_id(rxn_id).upper_bound)

    print(f"{rxn_id} is found in: {[model.id for model in model_list]}")

    if model_list:
        get_rxn(model_list[0], rxn_id, mass=mass, charge=charge)

    if bounds:
        print(f"lower bound: {lb}; upper bound: {ub}")

    if GPR:
        print(f"all GPRs: {set(gpr_list)}")


def get_met(model, met_id):
    if met_id in model.metabolites:
        met = model.metabolites.get_by_id(met_id)
        rxns = {rxn.id:model.reactions.get_by_id(rxn.id).reaction for rxn in met.reactions}
        print(f"{met.name} ({met.formula})")
        print(rxns)
        return

    else:
        query_mets = set()
        if "_" in met_id:
            met_id = met_id.split("_")[0]

        mets = [met.id for met in model.metabolites.query(met_id)]
        query_mets.update(set(mets))

        if query_mets:
            print(f"No exact match for your ID was found, but \"{met_id}\" got the following metabolites: {query_mets}")
            return

    print("No match.")


def get_met_unknown(models, met_id):
    model_list = []
    rxn_list = []
    for model in models.values():
        if met_id in model.metabolites:
            model_list.append(model)
            rxns = model.metabolites.get_by_id(met_id).reactions
            for rxn in rxns:
                rxn_list.append(rxn.id)

    met_name, met_formula, met_charge = ["" for _ in range(3)]
    if model_list:
        met_name = model_list[0].metabolites.get_by_id(met_id).name
        met_formula = model_list[0].metabolites.get_by_id(met_id).formula
        met_charge = model_list[0].metabolites.get_by_id(met_id).charge
        print(f"{met_id} ({met_name} ({met_formula}, {met_charge})) is found  models: {[model.id for model in model_list]} and in reactions: {set(rxn_list)}")
        return

    if len(model_list) == 0:
        query_mets = set()
        if "_" in met_id:
            met_id = met_id.split("_")[0]
        for model in models.values():
            mets = [met.id for met in model.metabolites.query(met_id)]
            query_mets.update(set(mets))

        if query_mets:
            print(f"No exact match for your ID was found, but \"{met_id}\" got the following metabolites: {query_mets}")
            return

    print("No match.")


###
# PART 2 - Simulations: change and test different media
###

def safe_parse(x):
    if x.startswith('['):
        # Add quotes around items inside brackets if they're not already quoted
        items = re.findall(r'\w+__?\w*', x)
        return items
    else:
        return x


def change_medium(model, medium_dict):

    # when I read the csv files with medium info it is saved as a dataframe, but i want a dict
    if isinstance(medium_dict, pd.DataFrame):
        medium_dict = dict(zip(medium_dict.reaction, medium_dict.bound))

    # Only include reactions that are in the model
    valid_medium = {
        rxn: bound for rxn, bound in medium_dict.items()
        if rxn in model.reactions
    }
    model.medium = valid_medium


# changes medium, does pfba and returns growth rate
def test_medium(model, medium_dict):
    with model:
        change_medium(model, medium_dict)
        try:
            solution = pfba(model)
            growth = solution.fluxes["Growth"]
            return growth
        except Infeasible:
            return None


# takes result dict from "create medium" to create a heatmap
def visualise_heatmap_medium(results):
    # Convert the results dictionary to a DataFrame
    df = pd.DataFrame(results).T  # Transpose the DataFrame

    # Create a heatmap
    plt.figure(figsize=(10, 6))
    sns.heatmap(df, annot=True, cmap="YlOrRd", cbar_kws={'label': 'Growth Value'})
    plt.title("Growth Values Heatmap")
    plt.xlabel("Model ID")
    plt.ylabel("Carbon Source")
    plt.show()


# takes a list with different carbon sources (EX reactions) and a list with reactions that together form a minimal medium
# each carbon source is coupled with the minimal medium one at the time and these media are then testes for growth
# model_dict contain the models where the media are tested on;
# it is visualised through a heatmap
def create_medium(carbon_list, minimal_list, model_dict, medium_uptake_bound):
    # Flatten the carbon_list for consistent keys
    flattened_sources = [item[0] if isinstance(item, list) else item for item in carbon_list]
    results = {carbon[3:-2]: {model.id: None for model in model_dict.values()} for carbon in flattened_sources}

    for carbon in carbon_list:
        # Standardize the carbon key (used to index the results dict)
        if isinstance(carbon, list):
            carbon_key = carbon[0][3:-2]
            new_medium = minimal_list + carbon
        else:
            carbon_key = carbon[3:-2]
            new_medium = minimal_list + [carbon]

        # define uptake bound here (this has big impact on how much biomass can be produced)
        med_dict = {new_medium[i]: medium_uptake_bound for i in range(len(new_medium))}

        for model in model_dict.values():
            growth_val = test_medium(model, med_dict)
            # print(growth_val)
            results[carbon_key][model.id] = growth_val
        #print("-----")

    visualise_heatmap_medium(results)