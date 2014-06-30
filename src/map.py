import numpy
import scipy.spatial
import igraph
import itertools

class Map(object):
    def __init__(self):
        self.polygons = igraph.Graph()
        self.vertices = igraph.Graph()
    
    def point_in_polygon(self, poly, x, y):
        """Checks if the point (x, y) is inside poly.
        """
        
        itr = itertools.cycle(self.polygons.vs[poly]["vertices"])
        pv = iv = next(itr)
        pcross = None
        for v in itr:
            edge = numpy.array((
                self.vertices.vs[v]["x"] - self.vertices.vs[pv]["x"],
                self.vertices.vs[v]["y"] - self.vertices.vs[pv]["y"]
            ))
            
            point = numpy.array((
                x - self.vertices.vs[pv]["x"],
                y - self.vertices.vs[pv]["y"]
            ))
            
            cross = numpy.cross(edge, point)
            
            if pcross is None:
                pcross = cross
            elif pcross*cross <= 0:
                return False
            
            pv = v
            if v == iv:
                return True
