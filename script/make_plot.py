import matplotlib.pyplot as plt
import numpy as np
import sys

def read_data(path, filename):
    data = np.loadtxt(path + filename)
    pos_x = data[:, 0]
    # pos_x_sample = []
    # for i, j in pos_x:
    #     if i % 100 == 0:
    #         pos_x_sample.append(j)
    return pos_x


path = "../src/"
filename = sys.argv[1]

pos_x = read_data(path, filename)
fig = plt.figure(1, dpi=500, figsize=(7.4, 4.8), facecolor="white")
ax = fig.add_subplot(111)

ax.plot(pos_x, color="orange")

plt.savefig(path + filename + "_pos_x.png")