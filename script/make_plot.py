import matplotlib.pyplot as plt
import numpy as np
import sys, time


def energy_plot(path, input_filename):
    data = np.loadtxt(path + input_filename)

    time = data[:, 0]
    potential = data[:, 1]
    kinetic = data[:, 2]
    total_energy = data[:, 3]

    fig = plt.figure(1, dpi=500, figsize=(7.4, 4.8), facecolor="white")
    ax = fig.add_subplot(111)

    ax.plot(time, potential, color="orange", label="potential", linestyle=":")
    ax.plot(time, kinetic, color="blue", label="kinetic", linestyle="-.")
    ax.plot(time, total_energy, label="total energy", color="green")

    ax.set_xlabel("time/s")
    ax.set_ylabel("energy/J")

    plt.legend()
    plt.savefig(path + input_filename[:-4] + ".png")
    plt.clf()


def particle_plot(path, input_filename):
    data = np.loadtxt(path + input_filename)
    
    time = data[:, 0]

    pos_x = data[:, 1]
    pos_y = data[:, 2]
    pos_z = data[:, 3]
    
    v_x = data[:, 4]
    v_y = data[:, 5]
    v_z = data[:, 6]
    
    a_x = data[:, 7]
    a_y = data[:, 8]
    a_z = data[:, 9]

    fig = plt.figure(1, dpi=500, figsize=(22.2, 4.8), facecolor="white")

    plt1 = plt.subplot(131)
    plt1.plot(time, pos_x, color="orange", label="position_x", linestyle=":")
    plt1.plot(time, pos_y, color="blue", label="position_y", linestyle="-.")
    plt1.plot(time, pos_z, label="position_z", color="green")
    plt1.set_xlabel("time/s")
    plt1.set_ylabel("position/m")
    plt1.legend()

    plt2 = plt.subplot(132)
    plt2.plot(time, v_x, color="orange", label="v_x", linestyle=":")
    plt2.plot(time, v_y, color="blue", label="v_y", linestyle="-.")
    plt2.plot(time, v_z, label="v_z", color="green")
    plt2.set_xlabel("time/s")
    plt2.set_ylabel("velocity/(m/s)")
    plt2.legend()

    plt3 = plt.subplot(133)
    plt3.plot(time, a_x, color="orange", label="a_x", linestyle=":")
    plt3.plot(time, a_y, color="blue", label="a_y", linestyle="-.")
    plt3.plot(time, a_z, label="a_z", color="green")
    plt3.set_xlabel("time/s")
    plt3.set_ylabel("acceleration/(m/s2)")
    plt3.legend()

    plt.savefig(path + input_filename[:-4] + ".png")
    plt.clf()


def rdf_plot(path, input_filename):
    data = np.loadtxt(path + input_filename)

    r = []
    for i in data[:, 0]:
        if i <= 2.5:
            r.append(i)
    g = data[:len(r), 1]

    fig = plt.figure(1, dpi=500, figsize=(7.4, 4.8), facecolor="white")
    ax = fig.add_subplot(111)

    ax.plot(r, g, color="orange", label="radial distribution function")

    ax.set_xlabel("distance")
    ax.set_ylabel("g(r) (Normalized)")

    plt.legend()
    plt.savefig(path + input_filename[:-4] + ".png")
    plt.clf()


if __name__ == "__main__":
    path = __file__[:-12] + "../output/"
    print("[MD Script]", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), "Plotting energy.csv")
    energy_plot(path, "energy.csv")
    print("[MD Script]", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), "Saved to energy.png")
    print("[MD Script]", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), "Plotting particle.csv")
    particle_plot(path, "particle.csv")
    print("[MD Script]", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), "Saved to particle.png")
    print("[MD Script]", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), "Plotting rdf.csv")
    rdf_plot(path, "rdf.csv")
    print("[MD Script]", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), "Saved to rdf.png")