"""
This module contains auxiliary functions to be used in Cobrapy model exploration
"""


def set_fixed_flux(r_id, val, model):
    """
    This function receives a reaction ID, a numeric value and a Cobrapy model and fixes the bounds of the reaction
     corresponding to the reaction ID with the value provided
    @param r_id: Valid reaction id of the model being provided
    @param val: numeric value (Ex: 100)
    @param model: Valid Cobrapy model object
    """
    r_obj = model.reactions.get_by_id(r_id)
    r_obj.bounds = (val, val)


def set_bounds(r_id, val_tuple, model):
    """
    This function receives a reaction ID, a tuple with numeric values and a Cobrapy model and fixes the upper and lower
    bound of the reaction corresponding to the reaction ID with the values provided in the tuple
    @param r_id: Valid reaction id of the model being provided
    @param val_tuple: A tuple with numeric value (Ex: (-100,100))
    @param model:  Valid Cobrapy model object
    """
    r_obj = model.reactions.get_by_id(r_id)
    r_obj.bounds = val_tuple


def set_fixed_flux_ratio(r_dict, model):
    """
    This function receives a dictionary containing two reaction IDs and their respective ratios, and implements them as
    a constraint in a valid model
    @param r_dict: Dictionary with the ratio to be established between the two reactions, with the reactions IDs as keys
    and the ratio as values (e.g {'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p':3,'RXN_961_p':1})
    @param model:  Valid Cobrapy model object
    """
    if len(r_dict) == 2:
        r_keys = r_dict.keys()
        r_values = r_dict.values()
        r_id1 = (list(r_keys))[0]
        r_obj1 = model.reactions.get_by_id(r_id1)
        r_v1 = (list(r_values))[0]
        r_id2 = (list(r_keys))[1]
        r_obj2 = model.reactions.get_by_id(r_id2)
        r_v2 = (list(r_values))[1]
        const = model.problem.Constraint(r_v1 * r_obj2.flux_expression - r_v2 * r_obj1.flux_expression, lb = 0, ub = 0)
        model.add_cons_vars(const)
        return const


def metabolite_data(model, save = True, name = 'metabolite_df.csv'):
    """
    This function receives a cobrapy model object and returns a Dataframe with the model metabolite information
    alongside writing it to a .csv file in the current working directory.

    :param
    model: A valid cobrapy model object
    save: A boolean value to wheter the pandas dataframe should be saved to a .csv file. Default is True.
    name: Name of the .csv file to store the dataframe. Default is 'metabolite_df.csv'
    :return:
    A pandas Dataframe with the model's metabolites ID, Name, Compartment and Formula
    A "metabolite_df.csv" file with the information in the pandas Dataframe
    """
    import pandas as pd
    import numpy as np
    import pathlib
    path = str(pathlib.Path().absolute()) + '\{}'.format(name)
    ID = []
    Name = []
    Compartment = []
    Formula = []
    for met in model.metabolites:
        ID.append(met.id)
        Name.append(met.name)
        Compartment.append(met.compartment)
        Formula.append(met.formula)
    df = pd.DataFrame(np.array([ID, Name, Compartment, Formula]))
    df = df.T
    df.columns = ['ID', 'Name', 'Compartment', 'Formula']
    if save:
        df.to_csv(r"{}".format(path), index = False, header = True)
    return df


def reaction_data(model, reaction_id, save = True, name = 'reaction_df.csv'):
    """
    This function receives a reaction ID and returns a Dataframe with the reaction metabolites information alongside
    writing it to a .csv file in the current working directory.

    :param
    model: A valid cobrapy model object;
    reaction_id: An ID for a reaction in the model in string format (Ex : "Bio_opt");
    save: A boolean value to wheter the pandas dataframe should be saved to a .csv file. Default is True
    :return:
    A pandas Dataframe with the reaction's metabolites ID, Name, Compartment, Formula and Stoichiometric coefficient;
    A "reaction_df.csv" file with the information in the pandas Dataframe
    """
    import pandas as pd
    import numpy as np
    ID = []
    Name = []
    Compartment = []
    Formula = []
    Coefficient = []
    import pathlib
    path = str(pathlib.Path().absolute()) + '\{}'.format(name)
    biomass = model.reactions.get_by_id(reaction_id)
    for met in biomass.metabolites:
        ID.append(met.id)
        Name.append(met.name)
        Compartment.append(met.compartment)
        Formula.append(met.formula)
        Coefficient.append(biomass.get_coefficient(met))
    df = pd.DataFrame(np.array([ID, Name, Compartment, Formula, Coefficient]))
    df = df.T
    df.columns = ['ID', 'Name', 'Compartment', 'Formula', 'Coefficient']
    if save:
        df.to_csv(r"{}".format(path), index = False, header = True)
    return df
