from src.rig import rig_component, rig_grow, rig
import numpy as np
import networkx as nx
import unittest
from fractions import Fraction

class TestGrowth(unittest.TestCase):

    def test_rig_grow1(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k'])
        G.add_edge('i','j', {'weight': 1})
        G.add_edge('j','k', {'weight': 10})
        res = rig_grow(G, {'i'}, {'k'})
        self.assertEquals(int(res), 22)

    def test_rig_grow2(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k'])
        G.add_edge('i','j', {'weight': 1})
        G.add_edge('j','k', {'weight': 10})
        res = rig_grow(G, {'i', 'j'}, {'k'})
        self.assertEquals(int(res), 21)

    def test_rig_grow3(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k', 'l'])
        G.add_edge('i','j', {'weight': 1})
        G.add_edge('j','k', {'weight': 10})
        G.add_edge('k','l', {'weight': 20})
        res = rig_grow(G, {'i', 'j'}, {'j', 'k', 'l'})
        self.assertEquals(res, Fraction(209, 6))

    def test_rig_grow4(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k', 'l'])
        G.add_edge('i','j', {'weight': 1})
        G.add_edge('j','k', {'weight': 10})
        G.add_edge('k','l', {'weight': 20})
        res = rig_grow(G, {'i'}, {'j', 'k', 'l'})
        self.assertEquals(res, Fraction(86, 3))

    def test_rig_grow5(self):
        # Tests the scenario where x is a subset of y
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k', 'l'])
        G.add_edge('i','j', {'weight': 1})
        G.add_edge('j','k', {'weight': 10})
        G.add_edge('k','l', {'weight': 20})
        res = rig_grow(G, {'i', 'l'}, {'i', 'j', 'k', 'l'})
        self.assertEquals(res, Fraction(21, 2))

if __name__=='__main__':
    unittest.main()
