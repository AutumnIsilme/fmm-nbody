import numpy as np
import dataclasses
from matplotlib import pyplot as plt
from itertools import product, combinations

@dataclasses.dataclass
class Box:
    root: np.array # The root coordinates of this box
    extent: float  # The size of this box
    parents: list
    has_child_boxes: bool
    child_boxes: list = None # Only exists for stem boxes
    bodies_in_box: list = None # Only exists for leaf boxes

    colleagues: list = dataclasses.field(default_factory=list)
    U: list = dataclasses.field(default_factory=list)
    V: list = dataclasses.field(default_factory=list)
    W: list = dataclasses.field(default_factory=list)
    X: list = dataclasses.field(default_factory=list)
    Y: list = dataclasses.field(default_factory=list)

    def __repr__(self) -> str:
        parent_str = f"{self.parents[-1].root}x{self.parents[-1].extent}" if len(self.parents) != 0 else "None"
        return f"Box[root={self.root}, extend={self.extent}, parent={parent_str}, has_child_boxes={self.has_child_boxes}, len(child_boxes)={len(self.child_boxes) if self.has_child_boxes else 0}, bodies_in_box={self.bodies_in_box}]"

def split(box: Box, bodies: np.array, bodies_per_box: int) -> Box:
    split_size = box.extent / 2.0

    parents = box.parents + [box] if box.parents else [box]
    sub_boxes = []
    sub_boxes.append(Box(box.root, split_size, parents, False))
    sub_boxes.append(Box(box.root + np.array([0., 0., split_size]), split_size, parents, False))
    sub_boxes.append(Box(box.root + np.array([0., split_size, 0.]), split_size, parents, False))
    sub_boxes.append(Box(box.root + np.array([split_size, 0., 0.]), split_size, parents, False))
    sub_boxes.append(Box(box.root + np.array([0., split_size, split_size]), split_size, parents, False))
    sub_boxes.append(Box(box.root + np.array([split_size, 0., split_size]), split_size, parents, False))
    sub_boxes.append(Box(box.root + np.array([split_size, split_size, 0.]), split_size, parents, False))
    sub_boxes.append(Box(box.root + np.array([split_size, split_size, split_size]), split_size, parents, False))

    sub_box_bodies = [[] for _ in range(8)]
    for body in bodies:
        normalised_position = body[:3] - box.root
        if normalised_position[0] < split_size and normalised_position[1] < split_size and normalised_position[2] < split_size:
            sub_box_bodies[0].append(body)
        elif normalised_position[0] < split_size and normalised_position[1] < split_size and normalised_position[2] >= split_size:
            sub_box_bodies[1].append(body)
        elif normalised_position[0] < split_size and normalised_position[1] >= split_size and normalised_position[2] < split_size:
            sub_box_bodies[2].append(body)
        elif normalised_position[0] >= split_size and normalised_position[1] < split_size and normalised_position[2] < split_size:
            sub_box_bodies[3].append(body)
        elif normalised_position[0] < split_size and normalised_position[1] >= split_size and normalised_position[2] >= split_size:
            sub_box_bodies[4].append(body)
        elif normalised_position[0] >= split_size and normalised_position[1] < split_size and normalised_position[2] >= split_size:
            sub_box_bodies[5].append(body)
        elif normalised_position[0] >= split_size and normalised_position[1] >= split_size and normalised_position[2] < split_size:
            sub_box_bodies[6].append(body)
        elif normalised_position[0] >= split_size and normalised_position[1] >= split_size and normalised_position[2] >= split_size:
            sub_box_bodies[7].append(body)
    
    for i, sub_box in enumerate(sub_boxes):
        sub_box.colleagues = sub_boxes[:]
        if len(sub_box_bodies[i]) > bodies_per_box:
            final_sub_boxes.append(split(sub_box, sub_box_bodies[i], bodies_per_box))
        elif len(sub_box_bodies[i]) == 0:
            sub_boxes[i] = None
        else:
            sub_box.bodies_in_box = sub_box_bodies[i]

    box.has_child_boxes = True
    box.child_boxes = sub_boxes
    return box

def create_octree(bodies: np.array, bodies_per_box: int) -> Box:

    box_0_root = np.array([-1., -1., -1.])
    box_0_extent = 2.0
    box_0 = Box(box_0_root, box_0_extent, None, False, None)
    if len(bodies) <= bodies_per_box:
        box_0.bodies_in_box = bodies
    else:
        box_0 = split(box_0, bodies, bodies_per_box)
    return box_0

def stratify_octree(root: Box) -> list:
    stratum = [[root]]
    i = 0
    while i < len(stratum):
        print(i)
        for box in stratum[i]:
            if box.has_child_boxes:
                if len(stratum) == i+1:
                    stratum.append([])
                stratum[i+1].extend(box.child_boxes)
        stratum[i] = filter(lambda b: b is not None, stratum[i])
        i += 1
    return stratum

def colleagify_octree(root: Box):
    box.child_boxes[1]

def graph_box(box, ax):
    if box.has_child_boxes:
        for child in box.child_boxes:
            if child != None:
                graph_box(child, ax)
    else:
        # Borrowed from https://stackoverflow.com/questions/11140163/plotting-a-3d-cube-a-sphere-and-a-vector
        r = [-1, 1]
        for s, e in combinations(np.array(list(product(r, r, r))), 2):
            if np.sum(np.abs(s-e)) == r[1]-r[0]:
                ax.plot3D(*zip(box.extent*(s+1)/2+box.root, box.extent*(e+1)/2+box.root), color="b")
        for body in box.bodies_in_box:
            ax.scatter(body[0], body[1], body[2])

def graph_octree(octree_root_box):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    graph_box(octree_root_box, ax)

    fig.show()