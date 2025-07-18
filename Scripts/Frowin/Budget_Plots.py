import matplotlib.pyplot as plt
from cobra.flux_analysis import pfba
from cobra.exceptions import Infeasible
import random
import pandas as pd

###
# FUNCTIONS
###

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


def fba_and_query(model, met_query, medium=None):
    #med = minimal_medium(model, 3, minimize_components=True)
    with model:
        try:
            if medium is not None:
                change_medium(model, medium)
            fluxes = pfba(model)
            if fluxes["Growth"] == 0:
                print("No growth")
        except Infeasible:
            print("Infeasible")
            return None, None
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


def make_budget_plot(all_reactions, save = False):
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


###
# QUERY - MAIN
###
# currently we're doing pFBA
# put the metabolite you want into the function; regex is possible
# you can also specify a certain medium in the function (I added that to the code and is my own function,
# so i dont know yet how versatile it is and if it will work with your code/medium
budget_mets, solution_frame = fba_and_query(AA3, "^gthrd_e") # AA3 is my model

###
# PLOT
###
all_reactions = calc_producer_consumer(AA3, budget_mets, solution_frame) # again put your model here instead of AA3
# always be aware that model is the same, otherwise it isnt working
budget_plot = make_budget_plot(all_reactions)
budget_plot