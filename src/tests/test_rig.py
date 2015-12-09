from src.rig import rig_component,  rig
import pandas as pd
import numpy as np
import networkx as nx
import unittest
from fractions import Fraction

class TestRig(unittest.TestCase):

    # This tests are broken
    def test_rig(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k'])
        G.add_edge('i','j', {'weight': .9})

        a1 = rig(G,
                 pd.Series({'i':1, 'j':1, 'k':1}),
                 pd.Series({'i':1, 'j':1, 'k':1}))
        a2 = rig(G,
                 pd.Series({'j':Fraction(1, 2),
                            'k':Fraction(1, 2)}),
                 pd.Series({'i':Fraction(1, 3),
                            'j':Fraction(1, 3),
                            'k':Fraction(1, 3)}))
        a3 = rig(G,
                 pd.Series({'j':1,'k':1}),
                 pd.Series({'i':1,'k':1}))
        a4 = rig(G,
                 pd.Series({'i':Fraction(1, 2),
                            'j':Fraction(1, 2)}),
                 pd.Series({'i':Fraction(1, 3),
                            'j':Fraction(1, 3),
                            'k':Fraction(1, 3)}))
        a5 = rig(G,
                 pd.Series({'i':1,'k':1}),
                 pd.Series({'i':1,'j':1}))

        a6 = rig(G,
                 pd.Series({'i':Fraction(1, 3),
                            'j':Fraction(1, 3),
                            'k':Fraction(1, 3)}),
                 pd.Series({'j':1}))
        a7 = rig(G,
                 pd.Series({'i':Fraction(1, 3),
                            'j':Fraction(1, 3),
                            'k':Fraction(1, 3)}),
                 pd.Series({'k':1}))
        a8 = rig(G,
                 pd.Series({'i':Fraction(1, 2),
                            'k':Fraction(1, 2)}),
                 pd.Series({'j':1}))
        a9 = rig(G,
                 pd.Series({'i':Fraction(1, 2),
                            'j':Fraction(1, 2)}),
                 pd.Series({'k':1}))

        self.assertLessEqual(a1, a2)
        self.assertLessEqual(a2, a3)
        self.assertLessEqual(a3, a4)
        self.assertLessEqual(a4, a5)
        self.assertLessEqual(a5, a6)
        self.assertLessEqual(a6, a7)
        self.assertLessEqual(a7, a8)
        self.assertLessEqual(a8, a9)

    def test_rig_component1(self):
        import pickle
        G = pickle.load(open('test.pickle', 'rb'))
        x = pd.read_csv('X', index_col=0, header=None, squeeze=True)
        y = pd.read_csv('Y', index_col=0, header=None, squeeze=True)
        maxW = 491534586766
        d = rig_component(G, x, y, maxW)

    def test_rig_component2(self):
        import pickle
        G = pickle.load(open('test2.pickle', 'rb'))
        x = pd.read_csv('X2', index_col=0, header=None, squeeze=True)
        y = pd.read_csv('Y2', index_col=0, header=None, squeeze=True)
        maxW = 491534586766
        d = rig_component(G, x, y, maxW)

    def test_simplex(self)
        G = pickle.load(open('test_directed.pickle', 'rb'))
        N = list(G)                                # nodes
        D = [G.node[u].get('demand', 0) for u in N]  # node demands
        W, _ = nx.network_simplex(G.to_directed())

        from itertools import chain, repeat
        from numpy import inf
        n = 81
        U, C, D = pickle.load(open('UCD', 'rb'))
        faux_inf = 3 * max(chain([sum(u for u in U if u < inf),
                                  sum(abs(c) for c in C)],
                                  (abs(d) for d in D))) or 1
        C.extend(repeat(faux_inf, n))
        U.extend(repeat(faux_inf, n))



    def test_rig_pairwise(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k'])
        G.add_edge('i','j', {'weight': .9})
        X = np.array([[1, 1, 1],
                      [1, 1, 2]])
        rig_pairwise(G, X)

if __name__=='__main__':
    unittest.main()
