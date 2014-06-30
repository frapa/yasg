import itertools

import pyglet

class MapView(object):
    def __init__(self, world):
        self.world = world
        render_method = 0
        
        # Create a pyglet vertex list for each polygon
        self.polygons_vl = []
        print "Polygon count:", len(self.world.polygons.vs) # debug
        for n in range(len(self.world.polygons.vs)):
            p = self.world.polygons.vs[n]
            
            verts = tuple(itertools.chain([p["x"]*800, p["y"]*600],
                [self.world.vertices.vs[p["vertices"][-1]]["x"]*800,
                self.world.vertices.vs[p["vertices"][-1]]["y"]*600],
                *[(self.world.vertices.vs[v]["x"]*800, self.world.vertices.vs[v]["y"]*600)
                    for v in p["vertices"]]))
            
            # Color polygon by height
            if render_method == 0:
                h = p["h"]
                
                if h == 0:
                    c = [0, 0, 255]
                else:
                    comp = int(h*255)
                    c = [comp, 55+int(h*200), comp]
                
                # All vertices of the same color
                colors = c * (2 + len(p["vertices"]))
            # Color vertices by height
            elif render_method == 1:
                colors = tuple(itertools.chain([int((p["h"])*255), 0, 255],
                    [int((self.world.vertices.vs[p["vertices"][-1]]["h"])*255), 0, 255],
                    *[[int((self.world.vertices.vs[v]["h"])*255), 0, 255]
                        for v in p["vertices"]]))
            # CORNER ORDER
            elif render_method == 2:
                l = len(p["vertices"])
                colors = tuple(itertools.chain([0, 0, 255], [0, 0, 0],
                    *[[int(255.0*n/l), 0, 0] for n, v in enumerate(p["vertices"])]))
            # Color by sea/land
            elif render_method == 3:
                h = p["h"]
                
                if h == 0:
                    c = [0, 0, 255]
                else:
                    c = [0, 255, 0]
                
                # All vertices of the same color
                colors = c * (2 + len(p["vertices"]))
            
            # Actually create list
            self.polygons_vl.append(pyglet.graphics.vertex_list(
                len(p["vertices"]) + 2,
                ('v2f', verts),
                ('c3B', colors)
            ))
    
    
    def draw(self):
        for pvl in self.polygons_vl:
            pvl.draw(pyglet.gl.GL_TRIANGLE_FAN)