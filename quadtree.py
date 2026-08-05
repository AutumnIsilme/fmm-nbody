import numpy as np
import dataclasses
from matplotlib import pyplot as plt
from itertools import product, combinations

@dataclasses.dataclass
class Box:
    root: np.array # The root coordinates of this box
    extent: float  # The size of this box
    ancestors: list = dataclasses.field(default_factory=list)
    has_child_boxes: bool = False
    child_boxes: list = None # Only exists for stem boxes
    bodies_in_box: list = None # Only exists for leaf boxes

    colleagues: list = dataclasses.field(default_factory=lambda:[None, None, None, None, None, None, None, None][:])
    U: list = dataclasses.field(default_factory=set)
    V: list = dataclasses.field(default_factory=set)
    W: list = dataclasses.field(default_factory=set)
    X: list = dataclasses.field(default_factory=set)
    #Y: list = dataclasses.field(default_factory=set)

    centre: np.array = None
    multipole_expansion: np.array = dataclasses.field(default_factory=lambda:np.zeros(25, complex))
    local_expansion: np.array = dataclasses.field(default_factory=lambda:np.zeros(25, complex))
    forces: list = dataclasses.field(default_factory=list)

    def __repr__(self) -> str:
        parent_str = f"{self.ancestors[-1].root}x{self.ancestors[-1].extent}" if len(self.ancestors) != 0 else "None"
        return f"Box[root={self.root}, extend={self.extent}, parent={parent_str}, has_child_boxes={self.has_child_boxes}, len(child_boxes)={len(self.child_boxes) if self.has_child_boxes else 0}, bodies_in_box={self.bodies_in_box}]"
    
    def __hash__(self):
        return hash((self.root[0], self.root[1], self.extent))
    
    def __eq__(self, other):
        if other is None:
            return False
        return hash((self.root[0], self.root[1], self.extent)) == hash((other.root[0], other.root[1], other.extent))

def split(box: Box, bodies: np.array, bodies_per_box: int) -> Box:
    split_size = box.extent / 2.0

    ancestors = box.ancestors + [box]
    sub_boxes = []
    sub_boxes.append(Box(box.root, split_size, ancestors, centre=box.root+split_size/2))
    sub_boxes.append(Box(box.root+np.array([0., split_size]), split_size, ancestors, centre=box.root+np.array([0., split_size])+split_size/2))
    sub_boxes.append(Box(box.root+np.array([split_size, 0.]), split_size, ancestors, centre=box.root+np.array([split_size, 0.])+split_size/2))
    sub_boxes.append(Box(box.root+np.array([split_size, split_size]), split_size, ancestors, centre=box.root+np.array([split_size, split_size])+split_size/2))

    sub_box_bodies = [[] for _ in range(8)]
    for body in bodies:
        normalised_position = body[:2] - box.root
        if normalised_position[0] < split_size and normalised_position[1] < split_size:
            sub_box_bodies[0].append(body)
        elif normalised_position[0] < split_size and normalised_position[1] >= split_size:
            sub_box_bodies[1].append(body)
        elif normalised_position[0] >= split_size and normalised_position[1] < split_size:
            sub_box_bodies[2].append(body)
        elif normalised_position[0] >= split_size and normalised_position[1] >= split_size:
            sub_box_bodies[3].append(body)
    
    for i, sub_box in enumerate(sub_boxes):
        sub_box.colleagues = sub_boxes[:]
        if len(sub_box_bodies[i]) > bodies_per_box:
            split(sub_box, sub_box_bodies[i], bodies_per_box)
        elif len(sub_box_bodies[i]) == 0:
            sub_boxes[i] = None
        else:
            sub_box.bodies_in_box = sub_box_bodies[i]

    box.has_child_boxes = True
    box.child_boxes = sub_boxes
    return box

def colleagify_quadtree(box: Box):
    if box.child_boxes[0]:
        box.child_boxes[0].colleagues = [
            box.child_boxes[1],
            box.child_boxes[3],
            box.child_boxes[2],
            box.colleagues[4].child_boxes[3] if box.colleagues[4] and box.colleagues[4].has_child_boxes else box.colleagues[4],
            box.colleagues[4].child_boxes[1] if box.colleagues[4] and box.colleagues[4].has_child_boxes else box.colleagues[4],
            box.colleagues[5].child_boxes[3] if box.colleagues[5] and box.colleagues[5].has_child_boxes else box.colleagues[5],
            box.colleagues[6].child_boxes[2] if box.colleagues[6] and box.colleagues[6].has_child_boxes else box.colleagues[6],
            box.colleagues[6].child_boxes[3] if box.colleagues[6] and box.colleagues[6].has_child_boxes else box.colleagues[6]
        ]
    if box.child_boxes[1]:
        box.child_boxes[1].colleagues = [
            box.colleagues[0].child_boxes[0] if box.colleagues[0] and box.colleagues[0].has_child_boxes else box.colleagues[0],
            box.colleagues[0].child_boxes[2] if box.colleagues[0] and box.colleagues[0].has_child_boxes else box.colleagues[0],
            box.child_boxes[3],
            box.child_boxes[2],
            box.child_boxes[0],
            box.colleagues[6].child_boxes[2] if box.colleagues[6] and box.colleagues[6].has_child_boxes else box.colleagues[6],
            box.colleagues[6].child_boxes[3] if box.colleagues[6] and box.colleagues[6].has_child_boxes else box.colleagues[6],
            box.colleagues[7].child_boxes[2] if box.colleagues[7] and box.colleagues[7].has_child_boxes else box.colleagues[7]
        ]
    if box.child_boxes[2]:
        box.child_boxes[2].colleagues = [
            box.child_boxes[3],
            box.colleagues[2].child_boxes[1] if box.colleagues[2] and box.colleagues[2].has_child_boxes else box.colleagues[2],
            box.colleagues[2].child_boxes[0] if box.colleagues[2] and box.colleagues[2].has_child_boxes else box.colleagues[2],
            box.colleagues[3].child_boxes[1] if box.colleagues[3] and box.colleagues[3].has_child_boxes else box.colleagues[3],
            box.colleagues[4].child_boxes[3] if box.colleagues[4] and box.colleagues[4].has_child_boxes else box.colleagues[4],
            box.colleagues[4].child_boxes[1] if box.colleagues[4] and box.colleagues[4].has_child_boxes else box.colleagues[4],
            box.child_boxes[0],
            box.child_boxes[1]
        ]
    if box.child_boxes[3]:
        box.child_boxes[3].colleagues = [
            box.colleagues[0].child_boxes[2] if box.colleagues[0] and box.colleagues[0].has_child_boxes else box.colleagues[0],
            box.colleagues[1].child_boxes[0] if box.colleagues[1] and box.colleagues[1].has_child_boxes else box.colleagues[1],
            box.colleagues[2].child_boxes[1] if box.colleagues[2] and box.colleagues[2].has_child_boxes else box.colleagues[2],
            box.colleagues[2].child_boxes[0] if box.colleagues[2] and box.colleagues[2].has_child_boxes else box.colleagues[2],
            box.child_boxes[2],
            box.child_boxes[0],
            box.child_boxes[1],
            box.colleagues[0].child_boxes[0] if box.colleagues[0] and box.colleagues[0].has_child_boxes else box.colleagues[0]
        ]
    
    if box.child_boxes[0] and box.child_boxes[0].has_child_boxes:
        colleagify_quadtree(box.child_boxes[0])
    if box.child_boxes[1] and box.child_boxes[1].has_child_boxes:
        colleagify_quadtree(box.child_boxes[1])
    if box.child_boxes[2] and box.child_boxes[2].has_child_boxes:
        colleagify_quadtree(box.child_boxes[2])
    if box.child_boxes[3] and box.child_boxes[3].has_child_boxes:
        colleagify_quadtree(box.child_boxes[3])

def create_quadtree(bodies: np.array, bodies_per_box: int) -> Box:
    box_0_root = np.array([-1., -1.])
    box_0_extent = 2.0
    box_0 = Box(box_0_root, box_0_extent, centre=box_0_root+box_0_extent/2)
    if len(bodies) <= bodies_per_box:
        box_0.bodies_in_box = bodies
    else:
        box_0 = split(box_0, bodies, bodies_per_box)
        colleagify_quadtree(box_0)
    
    return box_0

def stratify_quadtree(root: Box) -> list:
    strata = [[root]]
    i = 0
    while i < len(strata):
        for box in strata[i]:
            if box is not None and box.has_child_boxes:
                if len(strata) == i+1:
                    strata.append([])
                strata[i+1].extend(box.child_boxes)
        strata[i] = list(filter(lambda b: b is not None, strata[i]))
        i += 1
    return strata

def populate_list_1(leaves: list):
    for box in leaves:
        box.U.update(filter(lambda b: b is not None, [colleague if colleague is not None and not colleague.has_child_boxes else None for colleague in box.colleagues]))
        for cobox in box.U:
            cobox.U.update([box])

def populate_list_2(strata: list):
    i = 0
    while i < len(strata) - 1:
        #stratum = set(strata[i])
        for box in strata[i]:
            if box.has_child_boxes:
                colleague_children = sum([colleague.child_boxes if colleague is not None and colleague.has_child_boxes else [None] for colleague in box.colleagues], [])
                #set_colleagues = set(box.colleagues)
                #well_separated = stratum - set_colleagues
                #well_separated.remove(box)
                for child in box.child_boxes:
                    if child is not None:
                        child.V.update(filter(lambda b: b is not None and b not in child.U and b not in child.colleagues, colleague_children))
                        #child.Y.update(well_separated)
        i += 1

def populate_list_3_and_4(leaves: list):
    #for stratum in strata:
    #    for box in stratum:
    #        if len(box.ancestors) != 0:
    #            for colleague in box.ancestors[-1].colleagues:
    #                if colleague and not colleague.has_child_boxes and colleague not in box.U:
    #                    box.X.update([colleague])
    #                    colleague.W.update([box])
    def W(X, box, qualifying: set):
        l = []
        for i, child in enumerate(box.child_boxes):
            if child is None:
                continue
            if i in qualifying:
                l.append(child)
                child.X.update([X])
            elif child.has_child_boxes:
                l += W(X, child, qualifying)
        return l
    
    for box in leaves:
        if box.colleagues[0] is not None and box.colleagues[0].has_child_boxes:
            box.W.update(W(box, box.colleagues[0], set([1, 3])))
        if box.colleagues[1] is not None and box.colleagues[1].has_child_boxes:
            box.W.update(W(box, box.colleagues[1], set([1, 2, 3])))
        if box.colleagues[2] is not None and box.colleagues[2].has_child_boxes:
            box.W.update(W(box, box.colleagues[2], set([2, 3])))
        if box.colleagues[3] is not None and box.colleagues[3].has_child_boxes:
            box.W.update(W(box, box.colleagues[3], set([0, 2, 3])))
        if box.colleagues[4] is not None and box.colleagues[4].has_child_boxes:
            box.W.update(W(box, box.colleagues[4], set([0, 2])))
        if box.colleagues[5] is not None and box.colleagues[5].has_child_boxes:
            box.W.update(W(box, box.colleagues[5], set([0, 1, 2])))
        if box.colleagues[6] is not None and box.colleagues[6].has_child_boxes:
            box.W.update(W(box, box.colleagues[6], set([0, 1])))
        if box.colleagues[7] is not None and box.colleagues[7].has_child_boxes:
            box.W.update(W(box, box.colleagues[7], set([0, 1, 3])))

def leaves(strata: list):
    result_leaves = []
    for stratum in strata:
        for box in stratum:
            if not box.has_child_boxes:
                result_leaves.append(box)
    return result_leaves

def draw_box_square(box, ax, colour='b'):
    if box.has_child_boxes:
        for child in box.child_boxes:
            if child != None:
                draw_box_square(child, ax, colour)
    e = np.array([0, 0])
    f = np.array([1, 1])
    g = np.array([0, 1])
    h = np.array([1, 0])
    #print(box.extent, box.root, box.extent*e+box.root, box.extent*f+box.root)
    ax.plot([box.root[0], box.root[0]+box.extent], [box.root[1], box.root[1]], color=colour)
    ax.plot([box.root[0], box.root[0]], [box.root[1], box.root[1]+box.extent], color=colour)
    ax.plot([box.root[0]+box.extent, box.root[0]+box.extent], [box.root[1]+box.extent, box.root[1]], color=colour)
    ax.plot([box.root[0]+box.extent, box.root[0]], [box.root[1]+box.extent, box.root[1]+box.extent], color=colour)

def graph_box(box, ax, colour='b', particle_colour='dimgray'):
    if box.has_child_boxes:
        for child in box.child_boxes:
            if child != None:
                graph_box(child, ax, colour, particle_colour)
    else:
        # Borrowed from https://stackoverflow.com/questions/11140163/plotting-a-3d-cube-a-sphere-and-a-vector
        e = np.array([0, 0])
        f = np.array([1, 1])
        g = np.array([0, 1])
        h = np.array([1, 0])
        #print(box.extent, box.root, box.extent*e+box.root, box.extent*f+box.root)
        ax.plot([box.root[0], box.root[0]+box.extent], [box.root[1], box.root[1]], color=colour)
        ax.plot([box.root[0], box.root[0]], [box.root[1], box.root[1]+box.extent], color=colour)
        ax.plot([box.root[0]+box.extent, box.root[0]+box.extent], [box.root[1]+box.extent, box.root[1]], color=colour)
        ax.plot([box.root[0]+box.extent, box.root[0]], [box.root[1]+box.extent, box.root[1]+box.extent], color=colour)
        for body in box.bodies_in_box:
            ax.plot(body[0], body[1], marker=",", color=particle_colour)

def graph_quadtree(quadtree_root_box):
    fig = plt.figure()
    ax = fig.add_subplot()
    
    graph_box(quadtree_root_box, ax)

    fig.show()
