# -*- encoding: utf-8 -*-

import time
import random
from math import *

import numpy
import scipy.spatial
import igraph

import map as m
import perlin

from matplotlib import pyplot as plt

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

def gradient(centers, radii, x, y):
    """Given a set of centers and radii, calculate
    the value of the gradient at a given point (x, y).
    """
    
    value = 1.0
    for (cx, cy), r in zip(centers, radii):
        d = (x - cx)**2 + (y - cy)**2
        
        # If we are inside the circle
        if d < r**2:
            value -= 1 - (d / r**1.5)
            if value < 0:
                value = 0

    return value

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

def create_structure(world, vor):
    print "    Searching for polygons and vertices..."
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
    
    print "    Adding vertices..."
    # Second: add vertices of the regions that are insides
    dict_vert_node = {} # vertex index: igraph node
    world.vertices.add_vertices(len(verts))
    v_list = list(verts)
    for n, ((x, y), vi) in enumerate(zip(vor.vertices[v_list], v_list)):
        dict_vert_node[vi] = n
        
        world.vertices.vs[n]["x"] = x
        world.vertices.vs[n]["y"] = y
        world.vertices.vs[n]["polygons"] = []
    
    print "    Adding edges..."
    # Add edges between vertices
    for v1, v2 in vor.ridge_vertices:
        if v1 in verts and v2 in verts:
            world.vertices.add_edges((dict_vert_node[v1], dict_vert_node[v2]))
    
    print "    Adding polygons..."
    # Then we add the polygons to the map
    dict_poly_node = {} # region index: igraph node
    world.polygons.add_vertices(len(regions))
    for n, i in enumerate(regions):
        dict_poly_node[i] = n
        (x, y) = vor.points[i]
        
        world.polygons.vs[n]["x"] = x
        world.polygons.vs[n]["y"] = y
        
        # Get vertices of each polygon and link them
        vi = vor.regions[vor.point_region[i]]
        
        world.polygons.vs[n]["vertices"] = [dict_vert_node[i] for i in vi]
        
        # Do the opposite: link regions close to each vertex
        for i in vi:
            world.vertices.vs[dict_vert_node[i]]["polygons"].append(n)
    
    
    return 0
    print "    Adding connections..."
    # Add connections between polygons
    regions_set = set(regions)
    n = 0 # count
    for p1, p2 in vor.ridge_points:
        if p1 in regions_set and p2 in regions_set:
            # Little bit of mess with the indices!
            world.polygons.add_edges((dict_poly_node[p1], dict_poly_node[p2]))
            
            # Now add links between polygon-polygon connections and
            # vertex-vertex connections. See relation between Voronoi diagram
            # and Delaunay triangulation!
            vs_p1 = set(world.polygons.vs[dict_poly_node[p1]]["vertices"])
            vs_p2 = set(world.polygons.vs[dict_poly_node[p2]]["vertices"])
            # intersection
            v1, v2 = vs_p1 & vs_p2
            
            # Get edge between them
            eid = world.vertices.get_eid(v1, v2)
            # Link them
            world.polygons.es[n]["dual"] = eid
            world.vertices.es[eid]["dual"] = n
            
            n += 1

def set_elevation(world, seed):
    #, gradient_num=3, min_grad_radius=0.1, max_grad_radius=0.2):
    """Set elevation on the vertices and polygons centers of a map.
    This is done randomly using simplex noise to generate islands.
    """
    
    #minimum = max_grad_radius
    #maximum = 1 - max_grad_radius
    
    # Now get the coordinates of gradient_num points,
    # which are the center of the gradients
    #centers = minimum + numpy.random.rand(gradient_num, 2) * (maximum - minimum)
    #radii = min_grad_radius + numpy.random.rand(gradient_num) * max_grad_radius
    
    centers = [(0.5, 0.5)]
    radii = [0.4]
    
    # Create simplex noise object
    simplex = perlin.SimplexNoise()
    
    # Noise scale
    scale = 4
    
    # Set height of vertices
    for v in world.vertices.vs:
        h = (simplex.noise3(v["x"]*scale, v["y"]*scale, seed) * 0.8 + 0.2 - 
            gradient(centers, radii, v["x"], v["y"]))
        
        if 0 <= h <= 1:
            v["h"] = h**2
            v["land"] = True
        elif h < 0:
            v["h"] = 0
            v["land"] = False
    
    # Set height of polygons as the mean of the heights of the vertices
    for p in world.polygons.vs:
        h = numpy.mean(numpy.array([world.vertices.vs[v]["h"]
            for v in p["vertices"]]))
        
        p["h"] = h
        
        if h > 0:
            p["land"] = True
        else:
            p["land"] = False
    
    # Now create subgraphs for pathfinding
    world.land_verts = world.vertices.subgraph(
        world.vertices.vs.select(land_eq=True))
    world.sea_verts = world.vertices.subgraph(
        world.vertices.vs.select(land_eq=False))
    
    world.land_polys = world.polygons.subgraph(
        world.polygons.vs.select(land_eq=True))
    world.sea_polys = world.polygons.subgraph(
        world.polygons.vs.select(land_eq=False))
    
    # Set the downhill direction for every edge in the vertices graph.
    # This is useful to compute the direction of rivers
    for edge in world.vertices.es:
        h1 = world.vertices.vs[edge.source]["h"]
        h2 = world.vertices.vs[edge.target]["h"]
        
        if h1 == h2:
            edge["downhill"] = None
        elif h1 > h2:
            edge["downhill"] = edge.target
        else:
            edge["downhill"] = edge.source


def generate_rivers(world, seed, river_number=None, dryness=125):
    """Generate random rivers based on elevation data
    """
    
    # Compute number of rivers, if not given, using the number of
    # dry vertices and dryness
    if river_number is None:
        vert_num = len(world.land_verts.vs)
        
        river_number = int(vert_num / dryness)
    
    # Randomly choose river sources
    possible_sources = world.land_verts.vs.select(h_gt=0.3)
    sources = set([random.choice(possible_sources)
        for i in range(river_number)])
    
    

def generate_map(npoly=200, nrelax=2, seed=time.time()):
    """Generate a map with the given preset and number of polygons

    Arguments:
        npoly - Number of polygons
        nrelax - number of times to run the relaxation procedure

    Returns: map.Map instance
    """
    
    print "Creating polygons..."
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
        print "Relaxing [{}] ...".format(i+1)
        vor = relaxation(vor)

    # Now, lets create the map object!
    world = m.Map()

    print "Generating graph..."
    # Create graph structure
    create_structure(world, vor)

    #scipy.spatial.voronoi_plot_2d(vor)
    #plt.axes().set_xlim((0, 1))
    #plt.axes().set_ylim((0, 1))
    #plt.show()
    
    print "Creating elevation..."
    # Now set the elevation using perlin noise
    set_elevation(world, seed)
    
    print "Generating rivers..."
    # Generate rivers and other fresh water features
    generate_rivers(world, seed)

    print "Map generation successful!"
    return world