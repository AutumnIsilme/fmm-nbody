import numpy as np
from scipy.special import binom
from matplotlib import pyplot as plt
import numba
import cloud
import quadtree

import time

#np.random.seed(0)

def basic(N):
    start = time.time()
    bodies = cloud.generate_2d_bodies_uniform_random(N, 10)
    box = quadtree.Box([-1, -1], 2., bodies_in_box=bodies)

    box.forces = np.zeros((len(box.bodies_in_box), 2))

    constructed = time.time()
    for i, left in enumerate(box.bodies_in_box[:-1]):
        for j, right in enumerate(box.bodies_in_box[i+1:]):
            left_pos = left[:2]
            right_pos = right[:2]
            left_to_right = right_pos - left_pos
            norm = np.linalg.norm(left_to_right)
            left_to_right /= norm
            mass_product = left[-1] * right[-1]
            force = left_to_right * mass_product / norm**2
            box.forces[i] += force
            box.forces[i + 1 + j] -= force

    done = time.time()
    print(f"Total time: {done - start}")
    print(f"Construction time: {constructed - start}")
    print(f"Processing time: {done - constructed}")


def main(N=500):
    start = time.time()

    bodies = cloud.generate_2d_bodies_uniform_random(N, 10)
    bodies_tree = quadtree.create_quadtree(bodies, 5)
    strata = quadtree.stratify_quadtree(bodies_tree)
    leaves = quadtree.leaves(strata)
    quadtree.populate_list_1(leaves)
    quadtree.populate_list_2(strata)
    quadtree.populate_list_3_and_4(leaves)

    populated = time.time()

    '''
    # Debug: Pick a leaf and draw its lists
    def plot_lists(box, ax):
        if box.has_child_boxes:
            for colleague in box.colleagues:
                if colleague is None:
                    continue
                quadtree.draw_box_square(colleague, ax, 'magenta')
                ax.text(colleague.centre[0], colleague.centre[1], 'C', ha='center', va='center')

        for list_box in box.U:
            quadtree.draw_box_square(list_box, ax, 'red')
            ax.text(list_box.centre[0], list_box.centre[1], 'U', ha='center', va='center')
        for list_box in box.V:
            quadtree.draw_box_square(list_box, ax, 'darkorange')
            ax.text(list_box.centre[0], list_box.centre[1], 'V', ha='center', va='center')
        for list_box in box.W:
            quadtree.draw_box_square(list_box, ax, 'gold')
            ax.text(list_box.centre[0], list_box.centre[1], 'W', ha='center', va='center')
        for list_box in box.X:
            quadtree.draw_box_square(list_box, ax, 'greenyellow')
            ax.text(list_box.centre[0], list_box.centre[1], 'X', ha='center', va='center')
        #for list_box in box.Y:
        #    quadtree.draw_box_square(list_box, ax, 'lime')
        #    ax.text(list_box.centre[0], list_box.centre[1], 'Y', ha='center', va='center')

        quadtree.draw_box_square(box, ax)

    leaf = np.random.choice(strata[-1])
    boxes_to_show_lists = leaf.ancestors[1:] + [leaf]
    print(boxes_to_show_lists)

    fig, axs = plt.subplots(3, 2)
    axs_flat = axs.flatten()
    for i, box in enumerate(boxes_to_show_lists):
        axs_flat[i].set_ylim([-1., 1.])
        axs_flat[i].set_xlim([-1., 1.])
        quadtree.graph_box(bodies_tree, axs_flat[i], colour='dimgray')

        plot_lists(box, axs_flat[i])

    fig.show()

    leaf = np.random.choice(list(filter(lambda b: not b.has_child_boxes, strata[3])))
    boxes_to_show_lists = leaf.ancestors[1:] + [leaf]
    print(boxes_to_show_lists)

    fig, axs = plt.subplots(3, 2)
    axs_flat = axs.flatten()
    for i, box in enumerate(boxes_to_show_lists):
        axs_flat[i].set_ylim([-1., 1.])
        axs_flat[i].set_xlim([-1., 1.])
        quadtree.graph_box(bodies_tree, axs_flat[i], colour='dimgray')

        plot_lists(box, axs_flat[i])

    fig.show()
    '''

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

    step_2_1 = time.time()

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
    for stratum in strata[-2:1:-1]:
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
                        box.multipole_expansion[l] += child.multipole_expansion[k] * z_0**(l-k) * binom(l-1, k-1)
                    box.multipole_expansion[l] -= child.multipole_expansion[0] / l * z_0**l

    step_2_2 = time.time()

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
                z += -k * random_box.multipole_expansion[k] / (relative**(k+1))
            Z[j, i] = np.sqrt(z.real**2 + z.imag**2)

    im = ax.imshow(Z, interpolation='bilinear', cmap="RdYlBu",
                   origin='lower', extent=[-1, 1, -1, 1],
                   vmax=1000, vmin=0)

    fig.show()
    '''

    # Step 3
    for leaf in leaves: # Initialise forces on a leaf only once, not every time it is visited
        leaf.forces = np.zeros((len(leaf.bodies_in_box), 2))

    for leaf in leaves:
        leaf.forces = np.zeros((len(leaf.bodies_in_box), 2))
        for adjacent in leaf.U:
            for i, left in enumerate(leaf.bodies_in_box):
                for j, right in enumerate(adjacent.bodies_in_box):
                    left_pos = left[:2]
                    right_pos = right[:2]
                    left_to_right = right_pos - left_pos
                    norm = np.linalg.norm(left_to_right)
                    left_to_right /= norm
                    mass_product = left[-1] * right[-1]
                    force = left_to_right * mass_product / norm**2
                    leaf.forces[i] += force
                    adjacent.forces[j] -= force
            adjacent.U.remove(leaf)
        for i, left in enumerate(leaf.bodies_in_box[:-1]):
            for j, right in enumerate(leaf.bodies_in_box[i+1:]):
                left_pos = left[:2]
                right_pos = right[:2]
                left_to_right = right_pos - left_pos
                norm = np.linalg.norm(left_to_right)
                left_to_right /= norm
                mass_product = left[-1] * right[-1]
                force = left_to_right * mass_product / norm**2
                leaf.forces[i] += force
                leaf.forces[j + i + 1] -= force # j starts at index 0, not i + 1

    step_3 = time.time()


    # Step 4
    #loops = 0
    for stratum in strata[2:]:
        for box in stratum:
            for b_j in box.V:
                #loops += 1
                relative_centre = b_j.centre - box.centre
                z_0 = complex(relative_centre[0], relative_centre[1])
                #a = b_j.multipole_expansion
                
                recip_z_0 = 1/z_0
                neg_recip_z_0 = -1/z_0
                #recip_z_0_k = 1
                k = np.arange(1, p+1)
                #for k in range(1, p+1):
                    #recip_z_0_k *= recip_z_0
                    #local_expansion += b_j.multipole_expansion[k] * (-1)**k * recip_z_0**k
                box.local_expansion[0] += b_j.multipole_expansion[0] * np.log(-z_0) + np.sum(b_j.multipole_expansion[k] * np.power(neg_recip_z_0, k))
                #recip_z_0_l = 1
                l = np.arange(1, p+1)
                K, L = np.meshgrid(l, k)

                box.local_expansion[1:] += np.power(recip_z_0, l) * (-b_j.multipole_expansion[0] / l + np.sum(b_j.multipole_expansion[K] * np.power(neg_recip_z_0, K) * binom(L+K-1, K-1), axis=1))

                #for l in range(1, p+1):
                    #recip_z_0_l *= recip_z_0
                    #recip_z_0_k = 1
                    #box.local_expansion[l] += -recip_z_0**l * b_j.multipole_expansion[0]/l + recip_z_0**l * np.sum(b_j.multipole_expansion[k] * np.power(neg_recip_z_0, k) * binom(l+k-1, k-1))
                    #box.local_expansion[l] += recip_z_0**l * np.sum(b_j.multipole_expansion[k] * np.power(neg_recip_z_0, k) * binom(l+k-1, k-1))
                    #for k in range(1, p+1):
                        #recip_z_0_k *= recip_z_0
                        #local_expansion += recip_z_0**l * (-1)**k * b_j.multipole_expansion[k] * recip_z_0**k * binom(l+k-1, k-1)

    #print(f"loops: {loops}")

    step_4 = time.time()

    '''
    # Debug: Pick a box and draw its local expansion
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_ylim([-1., 1.])
    ax.set_xlim([-1., 1.])

    random_box = np.random.choice(list(filter(lambda b: len(b.V)>0, sum(strata, start=[]))))
    quadtree.graph_box(random_box, ax)
    for box in random_box.V:
        quadtree.graph_box(box, ax, particle_colour='lime')

    ax.add_artist(plt.Circle(random_box.centre, 1/np.sqrt(2) * random_box.extent, fill=False))

    delta = (0.0125+0.0125/2) / 4
    xs = ys = np.arange(-1., 1., delta)
    X, Y = np.meshgrid(xs, ys)
    Z = X
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            z = complex(x, y)
            relative = z - complex(random_box.centre[0], random_box.centre[1])
            z = 0
            for k in range(1, p+1):
                z += k*random_box.local_expansion[k]*relative**(k-1)
            Z[j, i] = np.sqrt(z.real**2 + z.imag**2)

    im = ax.imshow(Z, interpolation='bilinear', cmap="RdYlBu",
                   origin='lower', extent=[-1, 1, -1, 1],
                   vmax=1000, vmin=-1000)

    fig.show()
    '''

    # Step 5
    for leaf in leaves:
        for i, body in enumerate(leaf.bodies_in_box):
            z = complex(body[0], body[1])
            for box in leaf.W:
                relative = z - complex(box.centre[0], box.centre[1])
                force = box.multipole_expansion[0]/relative
                for k in range(1, p+1):
                    force += -k * box.multipole_expansion[k] / (relative**(k+1))
                leaf.forces[i] += np.array(force.real, force.imag)

    step_5 = time.time()


    # Step 6
    for stratum in strata[2:]:
        for box in stratum:
            if len(box.X) == 0:
                continue
            centre = complex(box.centre[0], box.centre[1])
            for big_leaf in box.X:
                for body in big_leaf.bodies_in_box:
                    z = centre - complex(body[0], body[1])
                    box.local_expansion[0] += body[-1] * np.log(-z)
                    for l in range(1, p+1):
                        box.local_expansion[l] -= body[-1]/(l*z**l)

    step_6 = time.time()


    # Step 7
    for stratum in strata[2:]:
        for box in stratum:
            if not box.has_child_boxes:
                continue
            centre = complex(box.centre[0], box.centre[1])
            for child in box.child_boxes:
                if child is None:
                    continue
                z = centre - complex(child.centre[0], child.centre[1])
                for l in range(0, p+1):
                    for k in range(l, p+1):
                        child.local_expansion[l] += box.local_expansion[k] * (-z)**(k-l) * binom(k, l)

    step_7 = time.time()


    # Step 8
    for leaf in leaves:
        centre = complex(leaf.centre[0], leaf.centre[1])
        for i, body in enumerate(leaf.bodies_in_box):
            z = complex(body[0], body[1]) - centre
            force = 0
            for l in range(0, p+1):
                force += l * leaf.local_expansion[l] * z**(l-1)
            leaf.forces[i] += np.array([force.real, force.imag])

    step_8 = time.time()

    end = time.time()

    print(f"Total time: {end-start}")
    print(f"Step 1: {populated - start}")
    print(f"Step 2.1: {step_2_1 - populated}")
    print(f"Step 2.2: {step_2_2 - step_2_1}")
    print(f"Step 3: {step_3 - step_2_2}")
    print(f"Step 4: {step_4 - step_3}")
    print(f"Step 5: {step_5 - step_4}")
    print(f"Step 6: {step_6 - step_5}")
    print(f"Step 7: {step_7 - step_6}")
    print(f"Step 8: {step_8 - step_7}")

if __name__ == "__main__":
    main()

