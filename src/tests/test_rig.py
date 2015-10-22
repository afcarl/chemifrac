from ..rig import retract, _shortest_distance
import numpy as np
import networkx as nx
import unittest


class TestRetact(unittest.TestCase):
    def setUp(self):
        self.g1 = nx.Graph()
        self.g1.add_nodes_from(['i', 'j', 'k'])
        self.g1.add_edge('i','j')
        self.g1['i']['j']['weight'] = 0.1
        self.g1.node['i']['presence'] = [1, 1]
        self.g1.node['j']['presence'] = [1, 1]
        self.g1.node['k']['presence'] = [1, 1]

    def retract_simple(self):
        pass

    def shortest_path1(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j'])
        G.add_edge('i','j')
        G['i']['j']['weight'] = 1
        G.node['i']['presence'] = [1, 1]
        G.node['j']['presence'] = [1, 1]
        dm = pd.DataFrame([[0, 1], [1, 1]],
                          index=['i','j'],
                          columns=['i','j'])

        d, dest = _shortest_distance(G, dm, 'i', 0)

        self.assertEquals(d, 0)
        self.assertEquals(dest, 'i')


    def shortest_path2(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j'])
        G.add_edge('i','j')
        G['i']['j']['weight'] = 1
        G.node['i']['presence'] = [1, 0]
        G.node['j']['presence'] = [0, 1]
        dm = pd.DataFrame([[0, 1], [1, 1]],
                          index=['i','j'],
                          columns=['i','j'])

        d1, dest1 = _shortest_distance(G, dm, 'i', 0)
        d2, dest2 = _shortest_distance(G, dm, 'j', 1)

        self.assertEquals(d1, 1)
        self.assertEquals(dest1, 'j')
        self.assertEquals(d2, 1)
        self.assertEquals(dest2, 'i')


    def shortest_path3(self):
        G = nx.Graph()
        G.add_nodes_from(['i', 'j', 'k'])
        G.add_edge('i','j')
        G.add_edge('j','k')
        G['i']['j']['weight'] = 1
        G['i']['j']['weight'] = 10
        G.node['i']['presence'] = [1, 0]
        G.node['j']['presence'] = [0, 0]
        G.node['k']['presence'] = [0, 1]
        dm = pd.DataFrame([[0, 1, 11],
                           [1, 0, 10],
                           [11, 10, 0]],
                          index=['i','j', 'k'],
                          columns=['i','j', 'k'])

        d1, dest1 = _shortest_distance(G, dm, 'i', 0)
        d2, dest2 = _shortest_distance(G, dm, 'k', 1)

        self.assertEquals(d1, 11)
        self.assertEquals(dest1, 'k')
        self.assertEquals(d2, 11)
        self.assertEquals(dest2, 'i')


if __name__=='__main__':
    unittest.main()
