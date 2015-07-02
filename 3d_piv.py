#!/usr/bin/env python

###
#
# This script creates a 3D vector map using data from
# the micro-PIV experiments. It averages over each height
# and combines this into a single graph.
#
# Author: Callum Kift
#
# Directory setup:
# main_dir > different_height_directories + height_file.txt> txt files for given height
#
###

import os
import collections
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plot, 'projection=3d'


def sqr(a):
    """
    Returns the square of the input.
    """
    return a * a


def read_hf():
    """
    Reads the height file and returns a list of the heights
    """
    hf = main_dir + "height_file.txt"
    height_list = []
    with open(hf, 'r') as f:
        for line in f:
                line = line.strip()
                column = line.split()
                if len(column) == 1:
                    height_list.append(float(column[0]))
                else:
                    print "Error: height file has wrong format!"
                    return

    return np.array(height_list)

def get_subdirectories(a_dir):
    """
    The input is a directory. Retrieves the sub-directories of the folder. Returns a list of
    subdirectories.
    """
    return [a_dir + name + "/" for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]


def get_files(a_dir):
    """
    The input is a directory. Retrieves all txt files from the directory. Returns a list of
    files
    """
    gf = []
    for file in os.listdir(a_dir):
        if file.endswith(".txt"):
            gf.append(a_dir + "/" + str(file))
    if len(gf) != 0:
        return gf
    else:
        print "Error: Cannot find TXT files in subdirectory!\n\t (%s)" % a_dir


def get_data(eh, file_list):
    """
    The input is a list of files. Reads all files from the list, averages the velocities for
    each position and returns the position and height with their corresponding averaged
    x-/y-/abs-velocities in an array.
    """
    x_pos = []
    y_pos = []
    x_vel = []
    y_vel = []
    z_vel = []
    unique_y = []

    # reading data
    for file in file_list:
        with open(file, 'r') as f:
            f.readline()  # Ignores first line
            for line in f:
                line = line.strip()
                column = line.split()
                if len(column) == 4:
                    if file == file_list[0]:
                        # Only takes position data from first file as the same in each file
                        x_pos.append(float(column[0]))
                        y_pos.append(float(column[1]))

                        x_vel.append(float(column[2]))
                        y_vel.append(float(column[3]))
                        z_vel.append(0.0)

                        if float(column[1]) not in unique_y:
                            unique_y.append(float(column[1]))
                    else:
                        x_vel.append(float(column[2]))
                        y_vel.append(float(column[3]))
                        z_vel.append(0.0)
                else:
                    print "Error: TXT file is not correct!"

    if eh == exp_h_list[-1]:
        print "All data read."

    uy = len(unique_y)

    # checks list lengths to ensure matching and then averages the velocities for all files
    # and then returns an array with position and average velocities
    if len(x_pos) == len(y_pos):
        pos_count = len(x_pos)
        if len(x_vel) == len(y_vel) and len(x_vel) == len(z_vel):
            vel_count = len(x_vel)
            nof = vel_count / pos_count  # equals number of files for each height
            ax_vel, ay_vel, az_vel = avg_data_each_h(nof, pos_count, x_vel, y_vel, z_vel)

            if make_sg:
                subgrid_array = sub_grid(uy, x_pos, y_pos, eh, ax_vel, ay_vel, az_vel)
                return subgrid_array
            else:
                z_pos = [eh] * len(x_pos)
                return zip(x_pos, y_pos, z_pos, ax_vel, ay_vel, az_vel)
        else:
            print "Error: different number of velocities!"
    else:
        print "Error: not all x-positions have a corresponding y-position!"


def avg_data_each_h(nof, lof, x_vel, y_vel, z_vel):
    """
    Averages the x-/y-/abs-velocities for each position at a certain height. Works by knowing
    how long each data set is and adding the corresponding lines together, before dividing by
    the number of data sets.
    """
    sx_vel = []
    sy_vel = []
    sz_vel = []

    for i in range(lof):
        sx = 0
        sy = 0
        sa = 0
        for j in range(nof):
            sx += x_vel[i + (j * nof)]
            sy += y_vel[i + (j * nof)]
            sa += z_vel[i + (j * nof)]

        sx_vel.append(sx)
        sy_vel.append(sy)
        sz_vel.append(sa)

    # checks lengths match and then averages them and returns the average velocities
    if len(sx_vel) == len(sy_vel) and len(sx_vel) == len(sz_vel):
        if len(sx_vel) == lof:
            ax_vel = np.array(sx_vel) / nof
            ay_vel = np.array(sy_vel) / nof
            az_vel = np.array(sz_vel) / nof

            return ax_vel, ay_vel, az_vel
        else:
            print "Error: summed velocity array is the wrong length!"
    else:
        print "Error: summed velocity data not matching!"


def sub_grid(unique_y, xpos, ypos, zpos, axvel, ayvel, azvel):
    """
    Reduces data by subdividing our roi into n*n grids, with each grid containing the average
    of the n*n velocities and positions.
    """

    n = sgs
    ssgh_array = []

    i = 0
    while i + n + ((n - 1) * unique_y) < len(xpos):
        # Makes sure that subgrid can be fromed
        sxp = 0
        syp = 0
        szp = 0
        sxv = 0
        syv = 0
        szv = 0

        for j in range(n):
            for k in range(n):
                sxp += xpos[i + j + (k * unique_y)]
                syp += ypos[i + j + (k * unique_y)]
                szp += zpos
                sxv += axvel[i + j + (k * unique_y)]
                syv += ayvel[i + j + (k * unique_y)]
                szv += azvel[i + j + (k * unique_y)]
        ssgh_array.append([sxp, syp, szp, sxv, syv, szv])

        if (i + n) < len(xpos):
            i += n
        else:
            pl = unique_y - (i % unique_y)
            i += pl + ((n - 1) * unique_y)

    return np.array(ssgh_array) / sqr(n)


def dict_to_array(dict_array):
    """
    Puts the dictionary-filled-array into plottable array.
    """
    plottable_array = []
    for k in dict_array:
        for i in range(len(dict_array[k])):
            plottable_array.append(dict_array[k][i])

    return np.array(plottable_array)


def plot_3d_vector(pa):
    """
    Plots a 3d vector graph
    """
    # Changeable variables
    al = 0.01  # arrow length
    rgba = (0.3, 0.3, 0.3, 0.8)  # rgba for panels
    lw = 1.5  # changes thickness of arrow

    X, Y, Z, U, V, W = zip(*pa)
    A = np.sqrt(np.power(X, 2) + np.power(Y, 2))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    q = ax.quiver(X[::peo3], Y[::peo3], Z[::peo3], U[::peo3], V[::peo3], W[::peo3], A,
                  length=al, lw=lw)
    q.set_array(np.random.rand(10))
    plt.colorbar(q)
    ax.w_xaxis.set_pane_color(rgba)
    ax.w_yaxis.set_pane_color(rgba)
    ax.w_zaxis.set_pane_color(rgba)
    ax.set_zlabel("Height")
    ax.set_title(r"$\mu$-PIV vector plot, %s, %s" % (shark_species, sample_area))

    plt.show()


def vector_plots_2d(dicti):
    """
    Creates single height (2D) arrays and sends them to be plotted.
    The plots include vector graphs for each height and an averaged velocity comparison for all
    heights.
    """
    mean_xs = []
    mean_ys = []
    err_xs = []
    err_ys = []

    hcount = 0
    for k in dicti:
        pa2d = []
        for i in range(len(dicti[k])):
            if dicti[k][i][2] == exp_h_list[hcount]:
                pa2d.append(
                    [dicti[k][i][0], dicti[k][i][1], dicti[k][i][3],
                     dicti[k][i][4]])
        pa2d = np.array(pa2d)
        mxv, myv, errxv, erryv = plot_2d_vector(exp_h_list[hcount], pa2d)
        mean_xs.append(mxv)
        mean_ys.append(myv)
        err_xs.append(errxv)
        err_ys.append(erryv)
        hcount += 1

    plot_2d_mean_roi(mean_xs, mean_ys, err_xs, err_ys)


def plot_2d_vector(eh, pa):
    """
    Plots a 2D vector graph
    """

    X, Y, U, V = zip(*pa)

    if plot2D:
        A = np.sqrt(np.power(X, 2.0) + np.power(Y, 2.0))
        fig = plt.quiver(X[::peo2], Y[::peo2], U[::peo2], V[::peo2], A)
        plt.colorbar(fig)
        plt.title(r"$\mu$-PIV vector plot at height %.3f, %s, %s" % (eh, shark_species, sample_area))
        plt.xlabel(
            r"Average velocity: (%.3f $\bar{x}$ + %.3f $\bar{y}$) $ms^{-1}$" % (np.mean(U), np.mean(V)))
        plt.show()

    if sem_bar:
        return np.mean(U), np.mean(V), np.std(U) / np.size(U), np.std(V) / np.size(V)
    elif sd_bar:
        return np.mean(U), np.mean(V), np.std(U), np.std(V)


def plot_2d_mean_roi(mxa, mya, errx, erry):
    """
    Plots the mean velocity for each height, with error bars.
    :param mxa: array containing mean x-velocity for each height
    :param mya: array containing mean y-velocity for each height
    :param errx: array containing error on mxa
    :param erry: array containing error on mxa
    :return: n/a
    """
    plt.errorbar(mxa, exp_h_list, xerr=errx, marker='o', label=r'<$v_x$>')
    plt.errorbar(mya, exp_h_list, xerr=erry, marker='o', label=r'<$v_y$>')
    plt.plot([0.0, 0.0], [np.amin(exp_h_list), np.amax(exp_h_list)], 'k-')
    plt.legend(loc=2)

    if sem_bar:
        plt.xlabel(r"Velocity ($ms^{-1}$), $\pm$sem")
    elif sd_bar:
        plt.xlabel(r"Velocity ($ms^{-1}$), $\pm$sd")

    plt.ylabel("Height")
    plt.title("Averaged velocity over roi")
    plt.show()


if __name__ == '__main__':

    main_dir = "/home/callumkift/Documents/sharks_dtu/micro_piv/20150701_x10_bonnet_back/"
    # main dir where all the subdirs with the data are

    shark_species = ""
    sample_area = ""
    exp_h_list = read_hf()  # vertical heights of PIV
    lehl = len(exp_h_list)

    make_sg = False # True -> makes subgrids
    sgs = 3

    plot3D = False # True -> plots 3D vector
    peo3 = 10  # plots every nth vector for the 3D plot

    plot2D = False # True -> plots 2D vector plots for all heights
    peo2 = 3  # plots every nth vector for the 2D plot

    # At least one must be True
    sem_bar = False # plots standard error on mean bars on 2d_mean_roi graph
    sd_bar = True # plot standard deviation bars on 2d_mean_roi graph
    # If both true, sem will be plotted

    if sem_bar or sd_bar:
        if lehl != 0:

            sub_dirs = get_subdirectories(main_dir)

            # Store files of each subdir in a dictionary
            height_file_dict = {}
            for dir in sub_dirs:
                height_file_dict["height{0}".format(dir)] = get_files(dir)
            height_file_dict = collections.OrderedDict(
                sorted(height_file_dict.items()))  # sorts dictionary by subdir



            # Stores an array filled with position and velocities for each height
            if len(height_file_dict) == lehl:
                print "\nreading and manipulating data ..."
                h_pos_vel_dict = {}
                hcount = 0
                for k in height_file_dict:

                    # Avoids bug when more than 10 height measurements
                    if hcount < 10:
                        fhc = "0" + str(hcount)
                    else:
                        fhc = str(hcount)

                    h_pos_vel_dict["height{0}".format(fhc)] = get_data(exp_h_list[hcount],
                                                                          height_file_dict[k])
                    hcount += 1
                h_pos_vel_dict = collections.OrderedDict(sorted(h_pos_vel_dict.items()))

                if plot3D:
                    pa = dict_to_array(h_pos_vel_dict)
                    plot_3d_vector(pa)

                vector_plots_2d(h_pos_vel_dict)

            else:
                print "\nError: height list does not match number of subdirectories containing files!"
        else:
            print "\nError: experimental height measurements not given!"
    else:
        print "\nError: please choose which error bars to plot"