import matplotlib.pyplot as plt
import numpy as np
import sys

def read_data(path, filename):
    data = np.loadtxt(path + filename)
    pos_x = data[:, 0]
    pos_y = data[:, 1]
    pos_z = data[:, 2]

    v_x = data[:, 3]
    v_y = data[:, 4]
    v_z = data[:, 5]

    a_x = data[:, 6]
    a_y = data[:, 7]
    a_z = data[:, 8]

    potential = data[:, 9]
    kinetic = data[:, 10]
    total_energy = data[:, 11]
    return potential, kinetic, total_energy


path = "../src/"
filename = sys.argv[1]

potential, kinetic, total_energy = read_data(path, filename)
fig = plt.figure(1, dpi=500, figsize=(7.4, 4.8), facecolor="white")
ax = fig.add_subplot(111)

ax.plot(potential, color="orange", label="potential", linestyle=":")
ax.plot(kinetic, color="blue", label="kinetic", linestyle="-.")
ax.plot(total_energy, label="total energy", color="green")

ax.set_xlabel("time")
ax.set_ylabel("energy")

plt.legend()
plt.savefig(path + filename + "_energy.png")