from ..rig import rig_component, rig_grow, rig
import numpy as np
import networkx as nx
import unittest
from fractions import Fraction

class TestGrowth(unittest.TestCase):

    # def test_growth_equal(self):
    #     G = nx.Graph()
    #     G.add_nodes_from(['i', 'j'])
    #     G.add_edge('i','j')
    #     G['i']['j']['weight'] = 1
    #     G.node['i']['presence'] = [1, 1]
    #     G.node['j']['presence'] = [1, 1]
    #     modG, d = growth(G)

    # def test_growth_one_hole(self):
    #     G = nx.Graph()
    #     G.add_nodes_from(['i', 'j'])
    #     G.add_edge('i','j')
    #     G['i']['j']['weight'] = 1
    #     G.node['i']['presence'] = [1, 0]
    #     G.node['j']['presence'] = [1, 1]
    #     modG, d = growth(G)

    # def test_growth_two_holes1(self):
    #     G = nx.Graph()
    #     G.add_nodes_from(['i', 'j'])
    #     G.add_edge('i','j')
    #     G.add_edge('k','j')
    #     G['i']['j']['weight'] = 1
    #     G['j']['k']['weight'] = 10
    #     G.node['i']['presence'] = [1, 0]
    #     G.node['j']['presence'] = [0, 0]
    #     G.node['k']['presence'] = [1, 1]
    #     modG, d = growth(G)

    # def growth_two_holes(self):
    #     G = nx.Graph()
    #     G.add_nodes_from(['i', 'j'])
    #     G.add_edge('i','j')
    #     G['i']['j']['weight'] = 1
    #     G.node['i']['presence'] = [1, 0]
    #     G.node['j']['presence'] = [0, 1]
    #     modG, d = growth(G)

    # def growth_two_holes2(self):
    #     G = nx.Graph()
    #     G.add_nodes_from(['i', 'j'])
    #     G.add_edge('i','j')
    #     G.add_edge('k','j')
    #     G['i']['j']['weight'] = 1
    #     G['j']['k']['weight'] = 10
    #     G.node['i']['presence'] = [1, 0]
    #     G.node['j']['presence'] = [0, 0]
    #     G.node['k']['presence'] = [0, 1]
    #     modG, d = growth(G)

    # def shortest_path1(self):
    #     G = nx.Graph()
    #     G.add_nodes_from(['i', 'j'])
    #     G.add_edge('i','j')
    #     G['i']['j']['weight'] = 1
    #     G.node['i']['presence'] = [1, 1]
    #     G.node['j']['presence'] = [1, 1]
    #     dm = pd.DataFrame([[0, 1],
    #                        [1, 0]],
    #                       index=['i','j'],
    #                       columns=['i','j'])

    #     d, dest = _shortest_distance(G, dm, 'i', 0)

    # def shortest_path2(self):
    #     G = nx.Graph()
    #     G.add_nodes_from(['i', 'j'])
    #     G.add_edge('i','j')
    #     G['i']['j']['weight'] = 1
    #     G.node['i']['presence'] = [1, 0]
    #     G.node['j']['presence'] = [0, 1]
    #     dm = pd.DataFrame([[0, 1],
    #                        [1, 0]],
    #                       index=['i','j'],
    #                       columns=['i','j'])

    #     d1, dest1 = _shortest_distance(G, dm, 'i', 0)
    #     d2, dest2 = _shortest_distance(G, dm, 'j', 1)

    # def shortest_path3(self):
    #     G = nx.Graph()
    #     G.add_nodes_from(['i', 'j', 'k'])
    #     G.add_edge('i','j')
    #     G.add_edge('j','k')
    #     G['i']['j']['weight'] = 1
    #     G['i']['j']['weight'] = 10
    #     G.node['i']['presence'] = [1, 0]
    #     G.node['j']['presence'] = [0, 0]
    #     G.node['k']['presence'] = [0, 1]
    #     dm = pd.DataFrame([[0, 1, 11],
    #                        [1, 0, 10],
    #                        [11, 10, 0]],
    #                       index=['i','j', 'k'],
    #                       columns=['i','j', 'k'])

    #     d1, dest1 = _shortest_distance(G, dm, 'i', 0)
    #     d2, dest2 = _shortest_distance(G, dm, 'k', 1)

    # def test_filter_edges1(self):
    #     G = nx.Graph()
    #     G.add_nodes_from(['i', 'j', 'k'])
    #     G.add_edge('i','j')
    #     G.add_edge('j','k')
    #     G['i']['j']['weight'] = 1
    #     G['i']['j']['weight'] = 10
    #     edges = G.edges(data='weight')
    #     filtered_edges = _filter_edges(edges, ['i'], ['j'])

    # def test_min_edge_cover(self):
    #     G = nx.Graph()
    #     G.add_nodes_from(['i', 'j', 'k'])
    #     G.add_edge('i','j', {'weight': 1})
    #     G.add_edge('j','k', {'weight': 10})

    def test_rig_grow1(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k'])
        G.add_edge('i','j', {'weight': 1})
        G.add_edge('j','k', {'weight': 10})
        res = rig_grow(G, {'i'}, {'k'})
        self.assertEquals(res, 22)

    def test_rig_grow2(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k'])
        G.add_edge('i','j', {'weight': 1})
        G.add_edge('j','k', {'weight': 10})
        res = rig_grow(G, {'i', 'j'}, {'k'})
        self.assertEquals(res, 21)

    def test_rig_grow3(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k', 'l'])
        G.add_edge('i','j', {'weight': 1})
        G.add_edge('j','k', {'weight': 10})
        G.add_edge('k','l', {'weight': 20})
        res = rig_grow(G, {'i', 'j'}, {'j', 'k', 'l'})
        self.assertEquals(res, Fraction(209, 372))

    def test_rig_grow4(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k', 'l'])
        G.add_edge('i','j', {'weight': 1})
        G.add_edge('j','k', {'weight': 10})
        G.add_edge('k','l', {'weight': 20})
        res = rig_grow(G, {'i'}, {'j', 'k', 'l'})
        self.assertEquals(res, Fraction(521, 372))

    def test_rig_grow5(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k', 'l'])
        G.add_edge('i','j', {'weight': 1})
        G.add_edge('j','k', {'weight': 10})
        G.add_edge('k','l', {'weight': 20})
        res = rig_grow(G, {'i', 'l'}, {'j', 'k'})
        self.assertEquals(res, Fraction(521, 372))


if __name__=='__main__':
    unittest.main()
