# -*- encoding: utf-8 -*-

import numpy
import scipy.spatial
import igraph

import map as m
import perlin

def check_region(vertices):
    """Checks if a region is completely contained
    inside the [0, 1]×[0, 1] area."""
    if -1 in vertices:
        return False
    else:
        # We must check that the vertices are inside the [0, 1]×[0, 1] area.
        for (x, y) in vertices:
            # If point is not inside...
            if not (0 < x < 1) or not (0 < y < 1):
                return False
        
        return True

def relaxation(voronoi):
    """Takes a scipy.spatial.Voronoi and returns a new voronoi
    diagram with points relaxed (= make polygons more regular).
    """
    new_points = []
    for n, r in enumerate(voronoi.point_region):
        # Indices of the Voronoi vertices forming each Voronoi region
        vi = voronoi.regions[r]
        v = voronoi.vertices[vi]
        
        # Region is contained in the [0, 1]×[0, 1] area...
        if check_region(v):
            # New point is the average of the region's vertices
            np = [numpy.mean(v[:, 0]), numpy.mean(v[:, 1])]
                
            new_points.append(np)
        else:
            # skip this region
            new_points.append(voronoi.points[n])

    new_points = numpy.array(new_points)
    return scipy.spatial.Voronoi(new_points)

def generate_map(npoly=200, nrelax=2):
    """Generate a map with the given preset and number of polygons

    Arguments:
        npoly - Number of polygons
        nrelax - number of times to run the relaxation procedure

    Returns: map.Map instance
    """

    # 'Centers' of the polygons
    points = numpy.random.rand(npoly, 2)

    # Create Voronoi diagram of the points
    vor = scipy.spatial.Voronoi(points)

    #import matplotlib.pyplot as plt
    #scipy.spatial.voronoi_plot_2d(vor)
    #plt.axes().set_xlim((0, 1))
    #plt.axes().set_ylim((0, 1))

    # Relax the diagram nrelax time, to make polygons more regular
    for i in range(nrelax):
        vor = relaxation(vor)

    # Now, lets create the map object!
    world = m.Map()

    # Add the polygons but only those that are inside the [0,1]x[0,1] area.
    # First we get the regions which are inside (with vertices!)
    regions = []
    verts = set()
    for n in range(len(vor.points)):
        vi = vor.regions[vor.point_region[n]]
        v = vor.vertices[vi]
        
        # If region is inside
        if check_region(v):
            regions.append(n)
            
            verts = verts.union(tuple(vi))
    
    # Second: add vertices of the regions that are insides
    dict_vert_node = {}
    world.vertices.add_vertices(len(verts))
    v_list = list(verts)
    for n, ((x, y), vi) in enumerate(zip(vor.vertices[v_list], v_list)):
        dict_vert_node[vi] = n
        
        world.vertices.vs[n]["x"] = x
        world.vertices.vs[n]["y"] = y
    
    # Add edges between vertices
    for v1, v2 in vor.ridge_vertices:
        if v1 in verts and v2 in verts:
            world.vertices.add_edges((dict_vert_node[v1], dict_vert_node[v2]))
    
    # Then we add the polygons to the map
    world.polygons.add_vertices(len(regions))
    for n, i in enumerate(regions):
        (x, y) = vor.points[i]
        
        world.polygons.vs[n]["x"] = x
        world.polygons.vs[n]["y"] = y
        
        # Get vertices of each polygon and link them
        vi = vor.regions[vor.point_region[i]]
        
        world.polygons.vs[n]["vertices"] = [world.vertices.vs[dict_vert_node[i]]
            for i in vi]

    # Add connections between polygons
    #for p1, p2 in vor.ridge_points:
    #    if p1 != -1 and p2 != -1:
    #        world.polygons.add_edges((p1, p2))

    #scipy.spatial.voronoi_plot_2d(vor)
    #plt.axes().set_xlim((0, 1))
    #plt.axes().set_ylim((0, 1))
    #plt.show()

    return world

generate_map()