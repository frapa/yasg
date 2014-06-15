import numpy
import scipy.spatial
import igraph

class Map(object):
    def __init__(self):
        self.polygons = igraph.Graph()
        self.vertices = igraph.Graph()
