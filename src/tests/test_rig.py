from src.rig import rig_component,  rig
import pandas as pd
import numpy as np
import networkx as nx
import unittest
from fractions import Fraction

class TestRig(unittest.TestCase):

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

    def test_rig_pairwise(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k'])
        G.add_edge('i','j', {'weight': .9})
        X = np.array([[1, 1, 1],
                      [1, 1, 2]])
        rig_pairwise(G, X)

if __name__=='__main__':
    unittest.main()
