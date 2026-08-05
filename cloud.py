import numpy as np

def lines(bodies: np.array):
    for body in bodies:
        yield f"{body[0]}, {body[1]}, {body[2]}, {body[3]}, {body[4]}\n"

def write_points_and_masses(bodies: np.array, filename: str):
    with open(filename, "w") as file:
        file.writelines(lines(bodies))

"""
Body format:
    [x, y, z, vx, vy, vz, m]
"""

def generate_points_and_masses_uniform_random_distribution(n_points: int, max_mass: float) -> np.array:
    data = []
    for i in range(n_points):
        data.append([np.random.random() * 2 - 1, np.random.random() * 2 - 1, np.random.random() * 2 - 1, 0., 0., 0., np.random.random() * max_mass])
    return np.array(data)


def generate_2d_bodies_uniform_random(n_points: int, max_mass: float) -> np.array:
    data = []
    for i in range(n_points):
        data.append([np.random.random() * 2 - 1, np.random.random() * 2 - 1, 0., 0., np.random.random() * max_mass])
    return np.array(data)

