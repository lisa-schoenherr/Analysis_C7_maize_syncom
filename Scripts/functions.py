"""
This python file contains functions that I use within my Jupyter Notebooks.
This python file focuses on functions that are used by multiple notebooks and
where it is easier to import them than to multiple copy and paste them into the notebooks.
Every Jupyter Notebooks still contains functions, but these are then more unique to the notebook.
"""

###
# PART 0 - Imports
###

# basics
import re
import random
import numpy as np
import pandas as pd

# plots
import matplotlib.pyplot as plt
from cobra import Reaction
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import seaborn as sns

# cobra
from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.exceptions import Infeasible, OptimizationError
from cobra.medium import minimal_medium

# scikit-learn (PCA)
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

###
# PART 1 - Get Metabolite & Reaction Info
###

def get_basic_amounts(models):
    """
    Get the number of reactions and metabolites for all models in a dictionary.
    :param models: dict with model_id as a key and model as value
    :return: print number of metabolites and reactions per model and average amount.
    """
    mean_met = 0
    mean_rxn = 0
    for model in models.values():
        rxns = len(model.reactions)
        mets = len(model.metabolites)
        print(f"{model.id} has {rxns} reactions and {mets} metabolites")
        mean_rxn += rxns
        mean_met += mets
    print(f"Mean number of reactions: {mean_rxn / len(models)}")
    print(f"Mean number of metabolites: {mean_met / len(models)}")


def get_rxn(model, rxn_id, bounds = False, mass = False, GPR=False, charge=False):
    """
    This functions returns information about a reaction in one model. The amount of information is controlled by the arguments.
    :param model: Cobrapy model object.
    :param rxn_id: reaction ID as string.
    :param bounds: bool, if True, the lower and upper bounds are returned.
    :param mass: bool, if True, the formula of each metabolite in the reaction is returned.
    :param GPR: bool, if True, the GPR of the reaction is returned.
    :param charge: bool, if True, the charge of each metabolite in the reaction is returned.
    """

    if rxn_id in model.reactions:
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
    else:
        print("No match.")


def get_rxn_unknown(models, rxn_id, bounds=False, mass = False, GPR=False, charge=False):
    """
    This functions returns information about a reaction for a dict of models. The amount of information is controlled by the arguments.
    :param models: dict with keys being names of models and the values being a COBRApy model object.
    :param rxn_id: reaction ID as string.
    :param bounds: bool, if True, the lower and upper bounds are returned.
    :param mass: bool, if True, the formula of each metabolite in the reaction is returned.
    :param GPR: bool, if True, the GPR of the reaction is returned.
    :param charge: bool, if True, the charge of each metabolite in the reaction is returned.
    """
    model_list = []
    gpr_list = []
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
    """
    This functions returns information about a metabolite in one model. The amount of information is controlled by the arguments.
    :param model: Cobrapy model object.
    :param met_id: metabolite ID as string.
    """
    if met_id in model.metabolites:
        met = model.metabolites.get_by_id(met_id)
        # rxns = {rxn.id:model.reactions.get_by_id(rxn.id).reaction for rxn in met.reactions}
        rxns = {
            rxn.id: "" if "Growth" == rxn.id else model.reactions.get_by_id(rxn.id).reaction
            for rxn in met.reactions
        }
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
    """
   This functions returns information about a metabolite for a dict of models. The amount of information is controlled by the arguments.
   :param models: dict with keys being names of models and the values being a COBRApy model object.
   :param met_id: metabolite ID as string.
   """
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


def get_identity_model(model):
    """
    Gets you the info if we have community or individual model.
    :param model: cobrapy model
    :return: "ind" or "com"
    """
    ind_model_ids = [f"AA{i}" for i in range(1, 8)] + [f"AA{i}f" for i in range(1, 8)]
    if model.id in ind_model_ids:
        return "ind"
    else:
        return "com"


def check_early_biomass_component(model, medium_dict, check_rxn_id):
    rxn_id = "objective_check"
    sad_mets = []

    check_mets = model.reactions.get_by_id(check_rxn_id).metabolites
    #check_mets = model.reactions.check_rxn_id.metabolites
    for met, flux in check_mets.items():
        if flux < 0: # only look at mets that are consumed
            #print(met)
            stoich = {met: flux}

            if rxn_id in model.reactions: # biomass test reaction is updated
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.subtract_metabolites(rxn.metabolites)
                rxn.add_metabolites(stoich)
            else: # first time reaction is created
                new_rxn = Reaction(id=rxn_id, name="objective reaction", lower_bound=0, upper_bound=1000)
                new_rxn.add_metabolites(stoich)
                model.add_reactions([new_rxn])

            model.objective = rxn_id

            with model:
                change_medium(model, medium_dict)
                try:
                    pfba_flux = pfba(model)
                    if pfba_flux[rxn_id] == 0:
                        sad_mets.append(met.id)
                except Infeasible:
                    print("Cannot get result because pfba is infeasible")

    if rxn_id in model.reactions:
        rxn = model.reactions.get_by_id(rxn_id)
        model.remove_reactions([rxn])

    print(sad_mets)

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


def e_to_m(medium):
    """
    Changes the suffix of the reaction IDs in a medium dict/df from '_e' to '_m'.
    :param medium: dict or dataframe with Exchange reaction IDs & uptake bounds or a list with only reaction IDs
    :return: medium dict with modified reaction IDs or a list if only list of IDs were given without bounds
    """
    if isinstance(medium, dict):
        # Convert reaction IDs ending in '_e' to '_m'
        medium_dict = {k.removesuffix('_e') + '_m' if k.endswith('_e') else k: v for k, v in medium.items()}
        return medium_dict

    elif isinstance(medium, pd.DataFrame):
        # Convert DataFrame to a dict with converted reaction IDs
        medium_dict = {rxn.removesuffix('_e') + '_m' if rxn.endswith('_e') else rxn: bound for rxn, bound in zip(medium["reaction"], medium["bound"])}
        return medium_dict

    elif isinstance(medium, list):
        medium_list = [rxn[:-2] + "_m" for rxn in medium]
        return medium_list

    else:
        return "Medium is wrong data type"


def change_medium(model, medium_dict):

    # when I read the csv files with medium info it is saved as a dataframe, but i want a dict
    if isinstance(medium_dict, pd.DataFrame):
        medium_dict = dict(zip(medium_dict.reaction, medium_dict.bound))

    #ind_model_ids = [f"AA{i}" for i in range(1, 8)] + [f"AA{i}f" for i in range(1, 8)] # old version because if this ever changes i need to change it in every function
    #if model.id not in ind_model_ids:  # community models
    model_ident = get_identity_model(model)
    if model_ident == "com":
        # adjust the medium suffices for community models with medium compartment
        medium_dict = {k.removesuffix('_e') + '_m' if k.endswith('_e') else k: v for k, v in medium_dict.items()}

    #print(medium_dict)

    # Only include reactions that are in the model
    valid_medium = {
        rxn: bound for rxn, bound in medium_dict.items()
        if rxn in model.reactions
    }
    model.medium = valid_medium


# changes medium, does pfba and returns growth rate
def test_medium(model, medium_dict, frac=1, min_growth=0):
    """
    Test a certain medium for one model and only returns the growth value.
    :param model:
    :param medium_dict:
    :param frac: Optional for community models (fraction for MICOM's cooperative tradeoff function)
    :param min_growth:
    :return: Growth value [int] or None if Infeasible
    """
    with model:
        change_medium(model, medium_dict)
        try:
            model_ident = get_identity_model(model)
            if model_ident == "com": # community models
                solution = model.cooperative_tradeoff(fluxes=True, pfba=True, fraction=frac, min_growth=min_growth).fluxes.transpose()
                growth = solution[solution.index.str.contains("Growth")].transpose()
                growth = growth[~ growth.index.str.contains("medium")]
                growth.index.name = "model"
                growth = growth["Growth"]
                return growth
            else: # individual models
                solution = pfba(model)
                growth = solution.fluxes["Growth"]
                return growth
        except Infeasible:
            return None
        except OptimizationError as e:
            return None


# takes result dict from "create medium" to create a heatmap
def visualise_heatmap_medium(results, save_path=None):
    # Convert the results dictionary to a DataFrame
    df = pd.DataFrame(results).T  # Transpose the DataFrame
    df = df.apply(pd.to_numeric, errors='coerce') # take care of non-numerical values
    # Create a heatmap
    fig = plt.figure(figsize=(10, 6))
    sns.heatmap(df, annot=True, cmap="YlOrRd", cbar_kws={'label': 'Growth Value'})
    plt.title("Growth Values Heatmap")
    plt.xlabel("Model ID")
    plt.ylabel("Carbon Source")

    if save_path:
        fig.savefig(save_path, format="svg", bbox_inches="tight")
        print(f"Heatmap saved to {save_path}")

    plt.show()
    return df


# takes a list with different carbon sources (EX reactions) and a list with reactions that together form a minimal medium
# each carbon source is coupled with the minimal medium one at the time and these media are then testes for growth
# model_dict contain the models where the media are tested on;
# it is visualised through a heatmap
def create_medium(carbon_list, minimal_list, model_dict, medium_uptake_bound, carbon_only = False, save_path=None, average="yes"):
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

        # define uptake boundaries
        if carbon_only:
            carbon_sources = carbon if isinstance(carbon, list) else [carbon]
            med_dict = {
                met: medium_uptake_bound if met in carbon_sources else 1000
                for met in new_medium}
        else:
            med_dict = {met: medium_uptake_bound for met in new_medium}

        for model in model_dict.values():
            growth_val = test_medium(model, med_dict)
            if isinstance(growth_val, pd.Series): # check if we have a community aka multiple growth values
                if average=="yes":
                    # average growth
                    if model.id == "C7":
                        growth_val= pd.DataFrame(growth_val).iloc[:, 0].sum()/7
                    else: # drop_out
                        growth_val= pd.DataFrame(growth_val).iloc[:, 0].sum()/6
                else:
                    # total growth
                    growth_val= pd.DataFrame(growth_val).iloc[:, 0].sum()
            results[carbon_key][model.id] = growth_val
            #print(carbon_key, growth_val)
        #print("-----")
    return visualise_heatmap_medium(results, save_path=save_path)


###
# PART 3 - Simulations: Uptakes and Secretions (Heatmaps)
###

def convert_cooptradeoff_into_fluxes(model, medium=None, frac=1):
    """
    Convert cooperative tradeoff into fluxes for a given community model.

    :param model: The metabolic model to compute cooperative
        tradeoff fluxes on (must be community model).
    :param medium: Optional medium for the model to consider
        during computation. Default is ``None``.
    :param frac: Fraction to modulate the cooperative tradeoff
        computation. Default is ``1``.
    :return: Flux values for reactions in the same style as the normal pfba function [pandas.Series].
    """

    if hasattr(model, "modification") and model.modification is not None:
        model.modification = None

    with model:
        if medium is not None:
            change_medium(model, medium)
        pfba_fluxes = model.cooperative_tradeoff(fluxes=True, pfba=True, fraction=frac).fluxes.transpose()
        df_long = pfba_fluxes.melt(ignore_index=False, var_name='model', value_name='flux')
        df_long = df_long.reset_index()

        df_long['reaction_id'] = np.where(
            df_long['model'] != 'medium',
            df_long['reaction'] + "__" + df_long['model'],
            df_long['reaction'])

        rxns_in_syncom = set(r.id for r in model.reactions)
        df_long = df_long[df_long['reaction_id'].isin(rxns_in_syncom)]
        pfba_fluxes = df_long.set_index('reaction_id')['flux']
        return pfba_fluxes


# this gets you long dataframe for a community
def get_pfba_fluxes(model, medium, frac=1):
    """
    Compute pfba fluxes for a given COBRApy model for a given medium.
    :param model: COBRApy model (individual or community)
    :param medium: dict or df
    :param frac: optional for community models, between 0 and 1; 1 is standard
    :return: series with reaction ids and fluxes
    """
    with model:
        change_medium(model, medium)
        model_ident = get_identity_model(model)

        try:
            if model_ident == "ind":
                pfba_fluxes = pfba(model).fluxes
            else:
                pfba_fluxes = convert_cooptradeoff_into_fluxes(model, frac=frac)

            return pfba_fluxes

        except Infeasible:
            return None


###
# PART 3.1 - Uptakes and Secretions for multiple individual models + 1 syncom
###

def get_fluxes_for_heatmap(model, medium, type, c7_all_ex, dict_with_rxn="no"):
    fluxes = get_pfba_fluxes(model, medium)

    if fluxes is None: # Infeasible solution
        print("Infeasible:", model.id)
        return {}

    if type == "uptake":
        filtered_fluxes = fluxes[(fluxes.index.str.startswith('EX_')) &
                                 (fluxes < 0)]
    elif type == "secret":
        filtered_fluxes = fluxes[(fluxes.index.str.startswith('EX_')) &
                                 (fluxes > 0)]
    else:
        return "Specify either \"uptake\" or \"secret\" as type"

    # this decides if we just include medium EX or all EX between the bacteria
    model_ident = get_identity_model(model)
    if model_ident == "com" and c7_all_ex == "no":
        filtered_fluxes = filtered_fluxes[(filtered_fluxes.index.str.endswith('_m'))]
        # so here we update filtered fluxes and kick EX_.*_e reactions out and only keep EX_.*_m

    # i need that option for later; i need the EX reaction names and not the metabolite names
    if dict_with_rxn == "yes":
        return filtered_fluxes

    flux_dict = {}
    for rxn_id, flux in filtered_fluxes.items():
        rxn = model.reactions.get_by_id(rxn_id)
        # Exchange reactions should have exactly one metabolite on one side
        met_id = list(rxn.metabolites.keys())[0].name  # us .id or use .name
        flux_dict[met_id] = flux

    #print(flux_dict.keys())
    return flux_dict


def uptake_secret_heatmap(all_models, type, initial_medium=None, c7_all_ex="no", vmin=1e-5, save_path=None, width=12):
    data = {}

    for model_name, model_obj in all_models.items():
        if initial_medium is None:
            medium = minimal_medium(model_obj, 2)
            medium = medium.rename("bound").reset_index().rename(columns={"index": "reaction"})
            flux_dict = get_fluxes_for_heatmap(model_obj, medium, type, c7_all_ex)
        else:
            flux_dict = get_fluxes_for_heatmap(model_obj, initial_medium, type, c7_all_ex)
        data[model_name] = flux_dict  # keys: model names, values: dict of metabolite: flux

    df = pd.DataFrame(data).fillna(0)
    df = df.abs() # absolute values, so even negative fluxes have right colour scale that |-1000| is bigger than 0
    df = df.replace(0, np.nan)

    #vmin = threshold: if all flux values for a metabolite are below, they dont show up in the plot to not screw the color scale too much
    df_filtered = df[df.max(axis=1) >= vmin]

    # Depending on Uptake or Secretion, choose different color scale
    color_map = {"uptake": "Blues", "secret": "Reds"}
    base_cmap_name = color_map[type]

    # Build custom colormap that starts with white
    base = sns.color_palette(base_cmap_name, 256).as_hex()
    custom_colors = ["#ffffff"] + base[1:]  # replace first entry with white
    custom_cmap = LinearSegmentedColormap.from_list(f"{base_cmap_name}_custom", custom_colors)

    # Plot heatmap
    plt.figure(figsize=(width, min(len(df_filtered) * 0.4, 32)))
    plt.grid(False)

    ax = sns.heatmap(
        df_filtered,
        cmap=custom_cmap,
        norm=LogNorm(vmin=vmin, vmax=df_filtered.max().max()),
        cbar_kws={"label": "Flux (log scale)"},
        mask=df_filtered == 0,  # hide zero-flux cells
        linewidths=0.01,
        linecolor="whitesmoke"
    )

    # Clean up axes and background grid artifacts
    ax.set_facecolor("white")  # fill background of masked cells with white
    sns.despine(left=True, bottom=True)  # remove spines
    ax.tick_params(length=0)  # remove tick marks

    if c7_all_ex == "no":
        plt.title(f"{type.capitalize()} fluxes across models (only sink/source)")
    if c7_all_ex == "yes":
        plt.title(f"{type.capitalize()} fluxes across models (sinks/sources and exchanges)")
    plt.xlabel("Models")
    plt.ylabel(f"{type.capitalize()} reactions")
    plt.tight_layout()
    #plt.show()

    if save_path is not None:
        plt.savefig(save_path, format="svg", bbox_inches="tight")

    """ without logarithmic scale
    # Plot heatmap
    plt.figure(figsize=(12, len(df) * 0.4))
    sns.heatmap(
        df,
        cmap=custom_cmap,
        vmin=0,  # ensures white is only used at exactly 0
        vmax=df.max().max(),  # max flux defines the top of the gradient
        linewidths=0.01,
        linecolor="whitesmoke",
        cbar_kws={"label": "Flux"}
    )
    plt.title(f"{type.capitalize()} fluxes across models")
    plt.xlabel("Models")
    plt.ylabel(f"{type.capitalize()} reactions")
    plt.tight_layout()
    plt.show()
    """

    return df_filtered



###
# PART 3.2 - Uptakes and Secretions for 1 community model
###

def sort_reactions_by_model_presence(df):
    seen = set()
    sorted_reactions = []

    for model_id in df.columns:
        present_reactions = df[df[model_id].notna()].index.tolist()
        new_reactions = [rxn for rxn in present_reactions if rxn not in seen]
        sorted_reactions.extend(new_reactions)
        seen.update(new_reactions)

    return df.loc[sorted_reactions]

def make_zero_white_cmap(base_cmap='Blues', n=256, mcolors=None):
    cmap = plt.get_cmap(base_cmap, n)
    colors = cmap(np.linspace(0, 1, n))
    mid = n // 2
    colors[mid] = np.array([1, 1, 1, 1])  # white for zero flux
    return mcolors.ListedColormap(colors)

def clean_ex_reaction_ids(index):
    """Remove EX_ prefix and _e suffix from exchange reaction IDs."""
    return index.str.replace(r"^EX_", "", regex=True).str.replace(r"_e$", "", regex=True)

def cleaner_ex_reaction_ids(index, model):
    return pd.Index([
        next(iter(model.reactions.query(rxn)[0].metabolites)).name
        for rxn in index])


def get_uptakes_within_community(model, all_ex_rxns, com_fluxes, epsilon, model_ids, sort_by_model=False):
    uptake_data = {model_id: [] for model_id in model_ids}

    flux_dict = dict(zip(com_fluxes['reaction_id'], com_fluxes['flux']))

    for base_rxn_id in all_ex_rxns:
        for model_id in model_ids:
            com_rxn_id = f"{base_rxn_id}__{model_id}"
            flux = flux_dict.get(com_rxn_id, None)

            # Uptake: flux significantly < 0
            if flux is not None and flux < -epsilon:
                uptake_data[model_id].append(flux)
            else:
                uptake_data[model_id].append(None)

    uptake_df = pd.DataFrame(uptake_data, index=all_ex_rxns)
    uptake_df.index.name = "reaction"

    # Prepare data for uptake
    filtered_uptake = uptake_df.dropna(how='all')
    if sort_by_model:
        sorted_uptake = sort_reactions_by_model_presence(filtered_uptake)
        plot_uptake = sorted_uptake.fillna(0)
    else:
        plot_uptake = filtered_uptake.fillna(0)

    # Clean reaction names for y-axis
    #plot_uptake.index = clean_ex_reaction_ids(plot_uptake.index.to_series())
    plot_uptake.index = cleaner_ex_reaction_ids(plot_uptake.index, model)

    return plot_uptake


def get_secretions_within_community(model, all_ex_rxns, com_fluxes, epsilon, model_ids, sort_by_model=False):
    secretion_data = {model_id: [] for model_id in model_ids}

    flux_dict = dict(zip(com_fluxes['reaction_id'], com_fluxes['flux']))

    for base_rxn_id in all_ex_rxns:
        for model_id in model_ids:
            com_rxn_id = f"{base_rxn_id}__{model_id}"
            flux = flux_dict.get(com_rxn_id, None)

            # Secretion: flux significantly > 0
            if flux is not None and flux > epsilon:
                secretion_data[model_id].append(flux)
            else:
                secretion_data[model_id].append(None)

    secretion_df = pd.DataFrame(secretion_data, index=all_ex_rxns)
    secretion_df.index.name = "reaction"


    # Prepare data for secretion
    filtered_secretion = secretion_df.dropna(how='all')
    if sort_by_model:
        sorted_secretion = sort_reactions_by_model_presence(filtered_secretion)
        plot_secretion = sorted_secretion.fillna(0)
    else:
        plot_secretion = filtered_secretion.fillna(0)

    # Clean reaction names for y-axis
    plot_secretion.index = cleaner_ex_reaction_ids(plot_secretion.index, model)

    return plot_secretion


def get_fluxes_within_community(com_model, medium, type, dict_with_ind_models, sort_by_model):
    model_ids = sorted(list(set([rxn.id[-3:] for rxn in com_model.reactions if rxn.id[-3:].startswith("AA")]))) # which models are present in community
    epsilon = 0.0001 # threshold for fluxes, below they are 0

    # Step 0: pFBA with Community Model
    with com_model:
        change_medium(com_model, medium)
        com_fluxes = convert_cooptradeoff_into_fluxes(com_model)
        com_fluxes = com_fluxes.rename("flux").reset_index().rename(columns={"index": "reaction_id"})

        com_fluxes["flux"] = com_fluxes["flux"].apply(lambda x: 0 if abs(x) < epsilon else x)

    # Convert index to strings explicitly
    com_fluxes.index = com_fluxes.index.astype(str)

    # Step 1: Collect all EX reaction IDs from individual models (without AA suffix)
    all_ex_rxns = set()
    for model_id in dict_with_ind_models.keys():
        model = dict_with_ind_models[model_id]
        ex_rxns = [rxn.id for rxn in model.reactions if rxn.id.startswith("EX_") and rxn.id.endswith("_e")]
        all_ex_rxns.update(ex_rxns)

    all_ex_rxns = sorted(all_ex_rxns)

    if type == "uptake":
        plot_uptake = get_uptakes_within_community(com_model, all_ex_rxns, com_fluxes, epsilon, model_ids, sort_by_model)
        return plot_uptake
    elif type == "secret":
        plot_secretion = get_secretions_within_community(com_model, all_ex_rxns, com_fluxes, epsilon, model_ids, sort_by_model)
        return plot_secretion
    else:
        return None # type was falsely specified


# OG Aufruf um Heatmap zu generieren
def heatmap_fluxes_withinCommunity(com_model, dict_with_ind_models, medium, type, vmin=0.001, save_path=None, width=12, sort_by_model=False):
    df_fluxes = get_fluxes_within_community(com_model, medium, type, dict_with_ind_models, sort_by_model)

    df = pd.DataFrame(df_fluxes).fillna(0)
    df = df.abs() # absolute values, so even negative fluxes have right colour scale that |-1000| is bigger than 0
    df = df.replace(0, np.nan)

    #vmin = 1e-5  # threshold: if all flux values for a metabolite are below, they dont show up in the plot to not screw the color scale too much
    df_filtered = df[df.max(axis=1) >= vmin]

    color_map = {"uptake": "Blues", "secret": "Reds"}
    base_cmap_name = color_map[type]

    # Build custom colormap that starts with white
    base = sns.color_palette(base_cmap_name, 256).as_hex()
    custom_colors = ["#ffffff"] + base[1:]  # replace first entry with white
    custom_cmap = LinearSegmentedColormap.from_list(f"{base_cmap_name}_custom", custom_colors)


    # Plot heatmap
    plt.figure(figsize=(width, len(df_filtered) * 0.4))
    plt.grid(False)

    ax = sns.heatmap(
        df_filtered,
        cmap=custom_cmap,
        norm=LogNorm(vmin=vmin, vmax=1000),
        cbar_kws={"label": "Flux (log scale)"},
        mask=df_filtered == 0,  # hide zero-flux cells
        linewidths=0.01,
        linecolor="whitesmoke"
    )

    # Clean up axes and background grid artifacts
    ax.set_facecolor("white")  # fill background of masked cells with white
    sns.despine(left=True, bottom=True)  # remove spines
    ax.tick_params(length=0)  # remove tick marks

    plt.title(f"{type} fluxes across members in community")
    plt.xlabel("Models")
    plt.ylabel(f"reactions")
    #plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, format="svg", bbox_inches="tight")

    plt.show()

    return df_filtered


###
# VISUALISATIONS - BUDGET PLOTS
###

def fba_and_query(model, met_query, medium):
    with model:
        fluxes = get_pfba_fluxes(model, medium)
        if fluxes is None:
            return None, None

    ident = get_identity_model(model)
    if ident == "com":
        solution_frame = fluxes.to_frame(name="fluxes")   # rename column to "fluxes"
        solution_frame.index.name = None
    else:
        solution_frame=fluxes.to_frame()

    budget_mets = []
    for met in model.metabolites.query(met_query):
        budget_mets.append(met)

    print(budget_mets)
    return budget_mets, solution_frame


#Remove reactions with negative flux from old list
def remove_items(test_list, item):
    res = [i for i in test_list if i != item]
    return res


def calc_producer_consumer(model, budget_mets, solution_frame):
    if budget_mets is None:
        return None

    #Defining list of reactions producing and consuming the metabolite
    consumers = []
    producers = []

    #Add reactions to respective list and exclude transport reactions
    for met in budget_mets:
        for reaction in model.reactions:
            if met in reaction.reactants:
                consumers.append(reaction.id)
            elif met in reaction.products:
                producers.append(reaction.id)

    #Get flux values from the simulation for metabolite consuming/producing reactions
    producers_df = solution_frame.loc[producers,:]
    consumers_df = solution_frame.loc[consumers,:]

    #Get values with negative flows: producing reactions with negative flow are consuming and vice-versa
    negative_producers = list(producers_df[producers_df["fluxes"] < 0].index)
    negative_consumers = list(consumers_df[consumers_df["fluxes"] < 0].index)

    #Add reactions to correct list
    consumers.extend(negative_producers)
    producers.extend(negative_consumers)

    for item in negative_producers:
        producers = remove_items(producers, item)

    for item in negative_consumers:
        consumers = remove_items(consumers, item)

    #Get flux values from the simulation for metabolite consuming/producing reactions (correct list)
    producers_df = solution_frame.loc[producers,:]
    #Make all values positive (disregard directionality)
    producers_df["fluxes"] = producers_df["fluxes"].abs()
    #Remove reactions with zero flux
    producers_df = producers_df[producers_df["fluxes"] != 0]

    #Get flux values from the simulation for metabolite consuming/producing reactions (correct list)
    consumers_df = solution_frame.loc[consumers,:]
    #Make all values positive (disregard directionality)
    consumers_df["fluxes"]  = consumers_df["fluxes"].abs()
    #Remove reactions with zero flux
    consumers_df = consumers_df[consumers_df["fluxes"] != 0]

    #Sum the flux values
    print("Sum of consumer fluxes: {}".format(consumers_df ["fluxes"].sum(axis=0)))
    print("Sum of producer fluxes: {}".format(producers_df ["fluxes"].sum(axis=0)))

    producers_df["Status"] = "Producer"
    consumers_df["Status"] = "Consumer"

    frame = [producers_df, consumers_df]
    all_reactions = pd.concat(frame)
    all_reactions["label"] = all_reactions.index

    return all_reactions


def make_budget_plot(all_reactions,budget_mets, save = False):
    if budget_mets is None:
        return None

    #Defining the nº of colors
    number_of_colors = len(all_reactions.index)

    #Getting a list of colors
    random.seed(177)
    color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                 for i in range(number_of_colors)]

    #Getting list of reactions
    reaction_list = list(all_reactions.index)

    #Build color dictionary
    color_dict = {}
    for i in range(len(reaction_list)):
        color_dict[reaction_list[i]] = color[i]


    """
    Plot the pivot table and barplot
    """
    mets = [met.id for met in budget_mets]

    chart = all_reactions.pivot_table(index="Status", columns="label", values="fluxes")
    chart.plot.bar(rot = 0, stacked = True, legend = True, ylabel = "Flux", color = color_dict)
    plt.legend(loc='best', bbox_to_anchor=(1.0, 0.5, 0.5, 0.5), ncol = 2)
    plt.title(f"Budget Plot for {mets}" )
    figsize = [11, 14] #To prevent the cropping of the image

    if save == True:
        plt.savefig('Budget_plot.svg', format='svg', bbox_inches = 'tight', dpi=600) #Line to save the image

    return chart.T

# all together
def budget_plot_allInOne(model, query_term, medium, print_rxns = False):
    budget_mets, solution_frame = fba_and_query(model, query_term, medium)
    if budget_mets is None:
        return "Infeasible or no Growth; no data to plot"
    all_reactions = calc_producer_consumer(model, budget_mets, solution_frame)

    if print_rxns:
        print(all_reactions)

    budget_plot = make_budget_plot(all_reactions, budget_mets)

# alternatively, you can make multiple calls:
# budget_mets, solution_frame = fba_and_query(AA4f, "^nad_c", medium_combined_naveed)
# all_reactions = calc_producer_consumer(AA4f, budget_mets, solution_frame)
# budget_plot = make_budget_plot(all_reactions, budget_mets)
# budget_plot


###
# EXTRAS
###

def pca(df, save_path=None):
    # Standardise the data
    scaler = StandardScaler()
    growth_df_filled = df.fillna(0)
    scaled_growth = scaler.fit_transform(growth_df_filled)

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_growth)

    # Create df for PCA result
    pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"], index=df.index)

    # Plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x="PC1", y="PC2", data=pca_df, s=100)

    # Annotate each point with the community name
    for name, row in pca_df.iterrows():
        if name != "C7_Community":
            plt.text(row["PC1"] + 0.07, row["PC2"]-0.04, name, fontsize=9)
        else:
            plt.text(row["PC1"] + -0.1, row["PC2"]+0.1, name, fontsize=9)

    # Add variance explained info
    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100
    plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
    plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")

    plt.title("PCA")
    plt.grid(True)
    if save_path is not None:
        plt.savefig(save_path, format="svg", bbox_inches="tight")
    plt.show()



#####
# More Community specific stuff
####

# similar to the get_pfba_fluxes function but slighly altered to directly also get growth values;
# also get_pfba_fluxes gets you a series aka one column with all reactions (reactions are exactly named how they are in the syncom,
# so we have rxn_AA1 and rxn_AA2 while with this function we habe rxn and col for AA1 and AA2)

def community_pfba(com_model, medium, frac=1):
    with com_model:
        change_medium(com_model, medium)
        try:
            fluxes = com_model.cooperative_tradeoff(pfba=True, fluxes=True, fraction=frac).fluxes.transpose()
        except (Infeasible, OptimizationError):
            print(f"Model {com_model.id} is infeasible.")
            return None, None

        growth = fluxes[fluxes.index.str.contains("Growth")].transpose()
        growth = growth[~ growth.index.str.contains("medium")]
        growth.index.name = "model"
        growth = growth["Growth"]

        return fluxes, growth


def test_fractions_community(model, medium, drop_out = None, medium_name = None, print_growth=False):
    fractions = np.arange(0, 1.1, 0.1)
    strain_codes = ["Sma", "Bpi", "Cpu", "Elu", "Cin", "Hro", "Ppu"]
    com_title = "community"

    if drop_out is not None:
        if drop_out in strain_codes:
            strain_codes = [s for s in strain_codes if s != drop_out]
            com_title = f"drop-out community (-{drop_out})"

    n_members = len(strain_codes)
    growth_vals_per_member = [[] for _ in range(n_members)]

    for frac in fractions:
        _, growth_vals = community_pfba(model, medium, frac)
        for i in range(n_members):
            if print_growth:
                print(frac, growth_vals)
            growth_vals_per_member[i].append(growth_vals[i])

    # Plot dots & dashed lines for each strain
    for i, code in enumerate(strain_codes):
        plt.plot(fractions, growth_vals_per_member[i],
                 marker='o', linestyle='--', label=code)

    plt.xlabel("Fraction")
    plt.ylabel("Growth value")
    if medium_name is not None:
        plt.title(f"Growth of {com_title} members across fractions\non {medium_name} medium")
    else:
        plt.title(f"Growth of {com_title} members across fractions")
    plt.legend()
    plt.grid(True)
    plt.show()
