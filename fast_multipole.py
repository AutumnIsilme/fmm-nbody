import numpy as np
from scipy.special import comb
from matplotlib import pyplot as plt
import cloud
import quadtree

#np.random.seed(0)

bodies = cloud.generate_2d_bodies_uniform_random(500, 10)
bodies_tree = quadtree.create_quadtree(bodies, 5)
strata = quadtree.stratify_quadtree(bodies_tree)
quadtree.populate_list_1(strata)
quadtree.populate_list_2_and_5(strata)
quadtree.populate_list_3_and_4(strata)
leaves = quadtree.leaves(strata)

# Debug: Pick a leaf and draw it's lists
def plot_lists(box, ax):

    if box.colleagues is not None:
        for colleague in box.colleagues:
            quadtree.draw_box_square(colleague, ax, 'red')
            ax.text(colleague.centre[0], colleague.centre[1], 'C', ha='center', va='center')

    #for list_box in box.U:
    #    quadtree.draw_box_square(list_box, ax, 'red')
    #    ax.text(list_box.centre[0], list_box.centre[1], 'U', ha='center', va='center')
    #for list_box in box.V:
    #    quadtree.draw_box_square(list_box, ax, 'darkorange')
    #    ax.text(list_box.centre[0], list_box.centre[1], 'V', ha='center', va='center')
    #for list_box in box.W:
    #    quadtree.draw_box_square(list_box, ax, 'gold')
    #    ax.text(list_box.centre[0], list_box.centre[1], 'W', ha='center', va='center')
    #for list_box in box.X:
    #    quadtree.draw_box_square(list_box, ax, 'greenyellow')
    #    ax.text(list_box.centre[0], list_box.centre[1], 'X', ha='center', va='center')
    #for list_box in box.Y:
    #    quadtree.draw_box_square(list_box, ax, 'lime')
    #    ax.text(list_box.centre[0], list_box.centre[1], 'Y', ha='center', va='center')

    quadtree.graph_box(box, ax)

leaf = np.random.choice(strata[-1])
boxes_to_show_lists = leaf.ancestors[1:] + [leaf]
print(boxes_to_show_lists)

fig = plt.figure()
for box in boxes_to_show_lists:
    ax = fig.add_subplot()
    ax.set_ylim([-1., 1.])
    ax.set_xlim([-1., 1.])
    quadtree.graph_box(bodies_tree, ax, colour='dimgray')

    plot_lists(box, ax)

fig.show()

epsilon = 0.0000001
p = int(np.ceil(-np.log2(epsilon)))
print(f"p: {p}")

# Step 2.1
for leaf in leaves:
    leaf.centre = leaf.root + leaf.extent/2
    for body in leaf.bodies_in_box:
        relative_position = body[:2] - leaf.centre
        z_i = complex(relative_position[0], relative_position[1])
        mass = body[-1]
        leaf.multipole_expansion[0] += mass
        for k in range(1, p+1):
            leaf.multipole_expansion[k] += -mass/k * (z_i**k)

'''
# Debug: Pick a leaf and draw its multipole expansion
fig = plt.figure()
ax = fig.add_subplot()
ax.set_ylim([-1., 1.])
ax.set_xlim([-1., 1.])

random_leaf = leaves[np.random.randint(0, len(leaves))]
quadtree.graph_box(random_leaf, ax)

ax.add_artist(plt.Circle(random_leaf.centre, 1/np.sqrt(2) * random_leaf.extent, fill=False))

delta = (0.0125+0.0125/2) / 4
xs = ys = np.arange(-1., 1., delta)
X, Y = np.meshgrid(xs, ys)
Z = X
for i, x in enumerate(xs):
    for j, y in enumerate(ys):
        z = complex(x, y)
        relative = z - complex(random_leaf.centre[0], random_leaf.centre[1])
        Z[j, i] = (random_leaf.multipole_expansion[0] * np.log(relative)).real
        for k in range(1, p+1):
            Z[j, i] += (random_leaf.multipole_expansion[k] / (relative**k)).real

im = ax.imshow(Z, interpolation='bilinear', cmap="RdYlBu",
               origin='lower', extent=[-1, 1, -1, 1],
               vmax=10, vmin=-10)

fig.show()
'''

# Step 2.2
for stratum in strata[-2::-1]:
    for box in stratum:
        if not box.has_child_boxes:
            continue
        for child in box.child_boxes:
            if child is None:
                continue
            relative_centre = child.centre - box.centre
            z_0 = complex(relative_centre[0], relative_centre[1])
            box.multipole_expansion[0] += child.multipole_expansion[0]
            for l in range(1, p+1):
                for k in range(1, l):
                    box.multipole_expansion[l] += child.multipole_expansion[k] * z_0**(l-k) * comb(l-1, k-1)
                box.multipole_expansion[l] -= child.multipole_expansion[0] / l * z_0**l

'''
# Debug: Pick a parent and draw the force field
fig = plt.figure()
ax = fig.add_subplot()
ax.set_ylim([-1., 1.])
ax.set_xlim([-1., 1.])

parents = list(filter(lambda b: b.has_child_boxes, strata[2]))
random_box = np.random.choice(parents)
quadtree.graph_box(random_box, ax)

ax.add_artist(plt.Circle(random_box.centre, 1/np.sqrt(2) * random_box.extent, fill=False))

delta = (0.0125+0.0125/2) / 4
xs = ys = np.arange(-1., 1., delta)
X, Y = np.meshgrid(xs, ys)
Z = X
for i, x in enumerate(xs):
    for j, y in enumerate(ys):
        z = complex(x, y)
        relative = z - complex(random_box.centre[0], random_box.centre[1])
        z = random_box.multipole_expansion[0]/relative
        for k in range(1, p+1):
            z += random_box.multipole_expansion[k] / (relative**(k+1))
        Z[j, i] = np.sqrt(z.real**2 + z.imag**2)

im = ax.imshow(Z, interpolation='bilinear', cmap="RdYlBu",
               origin='lower', extent=[-1, 1, -1, 1],
               vmax=1000, vmin=0)

fig.show()
'''

# Step 3
