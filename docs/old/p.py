import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial

points = np.random.rand(100, 2)

vor = scipy.spatial.Voronoi(points)

new_points = []
for r in vor.regions:
    v = vor.vertices[r]
    new_points.append([np.mean(v[:, 0]), np.mean(v[:, 1])])

new_points = np.array(new_points)
new_vor = scipy.spatial.Voronoi(points)

new_points = []
for r in vor.regions:
    v = new_vor.vertices[r]
    new_points.append([np.mean(v[:, 0]), np.mean(v[:, 1])])

new_points = np.array(new_points)
new_vor = scipy.spatial.Voronoi(points)

print(len(vor.points), len(vor.regions), len(vor.point_region), vor.point_region)

# drawings
f1 = plt.figure(figsize=(13, 13), dpi=65)
ax1 = f1.add_subplot(1, 1, 1)

ax1.errorbar(x=points[:,0], y=points[:,1], fmt='o', c='blue')
ax1.errorbar(x=vor.vertices[:,0], y=vor.vertices[:,1], fmt='.', c='red')

for i1, i2 in vor.ridge_vertices:
    ax1.errorbar(x=(vor.vertices[i1,0],vor.vertices[i2,0]), y=(vor.vertices[i1,1],vor.vertices[i2,1]), fmt='-', c='gray')

f2 = plt.figure(figsize=(13, 13), dpi=65)
ax2 = f2.add_subplot(1, 1, 1)

ax2.errorbar(x=new_points[:,0], y=new_points[:,1], fmt='o', c='blue')
ax2.errorbar(x=new_vor.vertices[:,0], y=new_vor.vertices[:,1], fmt='.', c='red')

for i1, i2 in new_vor.ridge_vertices:
    ax2.errorbar(x=(new_vor.vertices[i1,0], new_vor.vertices[i2,0]), y=(new_vor.vertices[i1,1], new_vor.vertices[i2,1]), fmt='-', c='gray')

# mostra grafico
plt.show()
