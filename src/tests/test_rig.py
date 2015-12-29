from src.rig import rig_component,  rig, rig_pairwise
import pandas as pd
import numpy as np
import networkx as nx
import unittest
from fractions import Fraction
import pickle

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
        G = pickle.load(open('test.pickle', 'rb'))
        x = pd.read_csv('X', index_col=0, header=None, squeeze=True)
        y = pd.read_csv('Y', index_col=0, header=None, squeeze=True)
        maxW = 491534586766
        d = rig_component(G, x, y, maxW)
        self(d, 491534586766)

    def test_rig_component2(self):
        G = pickle.load(open('test2.pickle', 'rb'))
        x = pd.read_csv('X2', index_col=0, header=None, squeeze=True)
        y = pd.read_csv('Y2', index_col=0, header=None, squeeze=True)
        maxW = 491534586766
        d = rig_component(G, x, y, maxW)
        self.assertEquals(2272109295852131273, d)

    def test_rig_component3(self):
        G = pickle.load(open('test3.pickle', 'rb'))
        x = pd.read_csv('X3', index_col=0, header=None, squeeze=True)
        y = pd.read_csv('Y3', index_col=0, header=None, squeeze=True)
        maxW = 491534586766

        d = rig_component(G, x, y, maxW)
        self.assertEquals(d, 750686345168434197)

    def test_rig_pairwise(self):

        from rig import rig_component, rig, rig_pairwise
        import pandas as pd
        import numpy as np
        import networkx as nx
        import unittest
        from fractions import Fraction
        import pickle
        from functools import partial
        from sklearn.metrics.pairwise import pairwise_distances

        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k'])
        G.add_edge('i','j', {'weight': .9})
        X = pd.DataFrame(np.array([[1, 1, 1],
                                   [1, 2, 1]]),
                          columns=['i', 'j', 'k'],
                          index=['s1', 's2'])

        labs = X.columns.values
        rig_func = partial(rig, G=G, labs=labs)
        # Note cannot work with pandas
        dm = pairwise_distances(X.values, metric=rig_func, n_jobs=1)

        dm = rig_pairwise(G, X)

        # This scenario is broken - 2 get normalized to 1.
        x = pd.Series([1, 1, 1], index=['i', 'j', 'k'])
        y = pd.Series([1, 1, 2], index=['i', 'j', 'k'])
        d = rig(x, y, G)


if __name__=='__main__':
    unittest.main()
