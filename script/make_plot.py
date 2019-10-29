import matplotlib.pyplot as plt
import numpy as np
import sys, getopt


def energy_plot(path, input_filename):
    data = np.loadtxt(path + input_filename)

    potential = data[:, 0]
    kinetic = data[:, 1]
    total_energy = data[:, 2]

    fig = plt.figure(1, dpi=500, figsize=(7.4, 4.8), facecolor="white")
    ax = fig.add_subplot(111)

    ax.plot(potential, color="orange", label="potential", linestyle=":")
    ax.plot(kinetic, color="blue", label="kinetic", linestyle="-.")
    ax.plot(total_energy, label="total energy", color="green")

    ax.set_xlabel("time")
    ax.set_ylabel("energy")

    plt.legend()
    plt.savefig(path + input_filename[:-4] + "_energy.png")


def particle_plot(path, input_filename):
    data = np.loadtxt(path + input_filename)
    pos_x = data[:, 0]
    pos_y = data[:, 1]
    pos_z = data[:, 2]
    
    v_x = data[:, 3]
    v_y = data[:, 4]
    v_z = data[:, 5]
    
    a_x = data[:, 6]
    a_y = data[:, 7]
    a_z = data[:, 8]

    fig = plt.figure(1, dpi=500, figsize=(22.2, 4.8), facecolor="white")

    plt1 = plt.subplot(131)
    plt1.plot(pos_x, color="orange", label="position_x", linestyle=":")
    plt1.plot(pos_y, color="blue", label="position_y", linestyle="-.")
    plt1.plot(pos_z, label="position_z", color="green")
    plt1.set_xlabel("time")
    plt1.set_ylabel("position")
    plt1.legend()

    plt2 = plt.subplot(132)
    plt2.plot(v_x, color="orange", label="v_x", linestyle=":")
    plt2.plot(v_y, color="blue", label="v_y", linestyle="-.")
    plt2.plot(v_z, label="v_z", color="green")
    plt2.set_xlabel("time")
    plt2.set_ylabel("velocity")
    plt2.legend()

    plt3 = plt.subplot(133)
    plt3.plot(a_x, color="orange", label="a_x", linestyle=":")
    plt3.plot(a_y, color="blue", label="a_y", linestyle="-.")
    plt3.plot(a_z, label="a_z", color="green")
    plt3.set_xlabel("time")
    plt3.set_ylabel("acceleration")
    plt3.legend()

    plt.savefig(path + input_filename[:-4] + "_movement.png")


def rdf_plot(path, input_filename):
    data = np.loadtxt(path + input_filename)

    r = data[:, 0]
    g = data[:, 1]

    fig = plt.figure(1, dpi=500, figsize=(7.4, 4.8), facecolor="white")
    ax = fig.add_subplot(111)

    ax.plot(r, g, color="orange", label="radial distribution function")

    ax.set_xlabel("distance")
    ax.set_ylabel("g(r) (Normalized)")

    plt.legend()
    plt.savefig(path + input_filename[:-4] + "_rdf.png")

if __name__ == "__main__":

    path = "../src/"

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:erp")
    except getopt.GetoptError:
        print('Error command line')
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-i":
            input_filename = arg
        elif opt == "-e":
            energy_plot(path, input_filename)
        elif opt == "-p":
            particle_plot(path, input_filename)
        elif opt == "-r":
            rdf_plot(path, input_filename)
