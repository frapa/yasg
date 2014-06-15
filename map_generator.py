# -*- encoding: utf-8 -*-

import numpy
import scipy.spatial
import igraph

import map as m
import perlin

borders = (
    ((0, 0), (0, 1)),
    ((0, 0), (1, 0)),
    ((1, 0), (1, 1)),
    ((0, 1), (1, 1))
)

def segment_intersection(l1, l2):
    (x1, y1), (ax, ay) = l1
    (x2, y2), (bx, by) = l2
    
    vx = ax - x1
    vy = ay - y1
    ux = bx - x2
    uy = by - y2
    
    s = ((y2 - y1) - ((x2 - x1) * vy) / vx) * vx / (ux*vy - uy*vx)
    t = (x2 - x1 + s*ux) / vx
    
    # check if segments collide
    if 0 <= t <= 1 and 0 <= s <= 1:
        return (x1 + vx*t, y1 + vy*t)
    else:
        return None

def relaxation(voronoi):
    """Takes a scipy.spatial.Voronoi and returns a new voronoi
    diagram with points relaxed (= make polygons more regular)
    """
    new_points = []
    for n, p in enumerate(voronoi.point_region):
        r = voronoi.regions[p] # Indices of the Voronoi vertices forming each Voronoi region
        v = voronoi.vertices[r]
        
        # We must check that the vertices are inside the [0, 1]×[0, 1] area
        # In case some of them are out we use the intersection of ridges with
        # the borders of the [0, 1]×[0, 1] box as vertices. To do this we
        # overwrite the `v` variable with new points.
        nv = []
        infinite_ridges_done = False
        for (x, y), i in zip(v, r):
            # If point is inside leave it as is.
            if 0 < x < 1 and 0 < y < 1 and i != -1:
                nv.append(numpy.array([x, y]))
            # If the point is at infinity...
            elif i == -1 and not infinite_ridges_done:
                # We want to get the other point of each semi-infinite ridge
                h = []
                for e in v:
                    # number of ridges with vertex e in region r
                    num_ridges_in_region = 0
                    for rg in voronoi.ridge_vertices:
                        # If the vertex is in ridge and the other vertex of the
                        # ridge is in the region, count the ridge as a region's
                        # ridge.
                        if e in rg and rg[1 - rg.index(e)] in r:
                            num_ridges_in_region += 1
                    
                    if num_ridges_in_region == 1:
                        h.append(e)
                
                # Now we find the regions on the other side of the 
            # Otherwise change it to intersection
            else:
                # Find other vertex that make the ridge
                rgs = [voronoi.vertices[rg]
                    for rg in voronoi.ridge_vertices if i in rg]
                
                # Find intersection with borders
                for b in borders:
                    for rg in rgs:
                        ints = segment_intersection(rg, b)
                        if ints is not None:
                            nv.append(ints)

        v = numpy.array(nv)
        
        # New point is the average of the region's vertices
        np = [numpy.mean(v[:, 0]), numpy.mean(v[:, 1])]
        
        new_points.append(np)

    new_points = numpy.array(new_points)
    return scipy.spatial.Voronoi(new_points)

def generate_map(npoly=100, nrelax=1):
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

    import matplotlib.pyplot as plt
    scipy.spatial.voronoi_plot_2d(vor)
    plt.axes().set_xlim((0, 1))
    plt.axes().set_ylim((0, 1))

    # Relax the diagram nrelax time, to make polygons more regular
    for i in range(nrelax):
        vor = relaxation(vor)

    # Now, lets create the map object!
    world = m.Map()

    # Add vertices
    world.vertices.add_vertices(len(vor.vertices))
    for n, (x, y) in enumerate(vor.vertices):
        world.vertices.vs[n]["x"] = x
        world.vertices.vs[n]["y"] = y

    # Add edges between vertices
    for v1, v2 in vor.ridge_vertices:
        if v1 != -1 and v2 != -1:
            world.vertices.add_edges((v1, v2))

    # Add the polygons
    world.polygons.add_vertices(len(vor.points))
    for n, (x, y) in enumerate(vor.points):
        world.polygons.vs[n]["x"] = x
        world.polygons.vs[n]["y"] = y
        
        # Get vertices of each polygon and link them
        region = vor.point_region[n]
        vs = vor.regions[region]

        #world.polygons[n]["vertices"] = list([world.vertices[i] for i in vs if i >= 0])

    # Add connections between polygons
    for p1, p2 in vor.ridge_points:
        if p1 != -1 and p2 != -1:
            world.polygons.add_edges((p1, p2))

    scipy.spatial.voronoi_plot_2d(vor)
    plt.axes().set_xlim((0, 1))
    plt.axes().set_ylim((0, 1))
    plt.show()

    return world

generate_map()