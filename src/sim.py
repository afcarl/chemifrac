"""
Simulation of metabolic networks
"""
import numpy as np
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt


def sample_components(dmG, multitons=10, singletons=0):
    """
    Randomly sample components of graph

    Parameters
    ----------
    dmG : nx.Graph
        Network to sample components from
    multitons : int
        Total number of components with more than 1 node
        to sample
    singletons : int
        Total number of singletons to sample

    Returns
    -------
    sampled_graphs : list, nx.Graph
        List of sampled subgraphs
    """
    connected_graphs = list(nx.connected_component_subgraphs(dmG))
    Gsingletons = [connected_graphs[i] for i in range(len(connected_graphs))
                  if len(connected_graphs[i])==1]
    Gmultitons = [connected_graphs[i] for i in range(len(connected_graphs))
                 if len(connected_graphs[i])>1]
    single_idx = np.random.choice(len(Gsingletons), singletons)
    multi_idx = np.random.choice(len(Gmultitons), multitons)
    sampled_graphs = [Gsingletons[i] for i in single_idx]
    sampled_graphs += [Gmultitons[i] for i in multi_idx]
    return sampled_graphs


def generate_sample(pvals, sampling_depth=1000):
    """
    Generates a sample

    Parameters
    ----------
    pvals : pd.Series
       proportions for each node
    sampling_depth : int
       number of draw per node

    Returns
    -------
    pd.Series
       Array of counts and their associated node id
    """
    cnts = np.random.multinomial(sampling_depth, pvals.values)
    return pd.Series(cnts, index=pvals.index)


def merge(graphs):
    """
    Unions a bunch of graphs together

    Parameters
    ----------
    graphs : list nx.Graph

    Returns
    -------
    nx.Graph
    """
    sim_graph = graphs[0]
    for i in range(1, len(graphs)):
        sim_graph = nx.union(graphs[i], sim_graph)
    return sim_graph

def draw_graph(dmG, layout):
    """
    Draws a graph and saves it to a file

    Parameters
    ----------
    dmG : nx.Graph
        Graph to draw
    layout : function
        networkx layout to use
    fname : str
        file to save graph visualization to

    Returns
    -------
    fig : plt.figure
    """
    fig = plt.figure()
    labels = nx.nodes(sim_graph)
    labels = dict( zip(labels, labels) )
    pos=layout(sim_graph)
    nx.draw_networkx_nodes(sim_graph, pos)
    nx.draw_networkx_edges(sim_graph, pos, width=1.0, alpha=0.5)
    nx.draw_networkx_labels(sim_graph, pos, labels)
    edge_labels = dict([((u,v,), str(np.around(d['weight'], decimals=2)))
                        for u,v,d in sim_graph.edges(data=True)])
    nx.draw_networkx_edge_labels(sim_graph, pos, edge_labels=edge_labels)
    return fig
