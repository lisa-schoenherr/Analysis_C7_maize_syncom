# this small script can visualise reaction/metabolite networks based on a list of reaction IDs
# this is only suitable for a subset of reactions (dont visualise your whole model with it lol)

# pip install pyvis
from pyvis.network import Network

def visualize_network(rxn_list, model):
    rxns = set(rxn_list) # specify list of reactions that you want to visualise, e.g. ["rxn1", "rxn2", ...]
    net = Network(height="700px", width="100%", notebook=True, directed=True)

    for rxn in model.reactions: # AA7 is my model, change that
        if rxn.id in rxns:
            # Add reactants → reaction node
            for met in rxn.reactants:
                net.add_node(met.id, label=met.id, color="lightblue", shape="ellipse")
                net.add_node(rxn.id, label=rxn.id, color="red", shape="box")
                net.add_edge(met.id, rxn.id, color='black') # reactant edges
            # Add reaction node → products
            for met in rxn.products:
                net.add_node(met.id, label=met.id, color="lightblue", shape="ellipse")
                net.add_edge(rxn.id, met.id, color='black') # product edges

    # Display in notebook (this file is also exported and then can be opened and examined your browser)
    net.show("network_vis.html") # it will save the html in your current folder where your script is
