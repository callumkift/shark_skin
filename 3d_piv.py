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
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plot, 'projection=3d'


def sqr(a):
    """
    Returns the square of the input.
    :param a: integer or float
    :return: param a squared
    """
    return a * a


def read_hf():
    """
    Reads the height file
    :return: an array containing the heights of the experiment
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
    Retrieves the sub-directories of the given directory.
    :param a_dir: path to directory
    :return: A list of the subdirectories.
    """
    return [a_dir + name + "/" for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]


def get_files(a_dir):
    """
    Creates a list of all txt files in the given directory
    :param a_dir: Path to a directory
    :return: a list of all txt files in the directory
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
    Reads the data in the list of files, averages them and has the possibility of making
    subgrid values.
    :param eh: Experimental height
    :param file_list: List of files corresponding to the given height
    :return: An array of positions and velocities
    """
    x_pos = []
    y_pos = []
    x_vel = []
    y_vel = []
    z_vel = []
    unique_x = []
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

                        if float(column[0]) not in unique_y:
                            unique_x.append(float(column[0]))
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

    ux = len(unique_x)
    xmid = np.median(unique_x)
    ymid = np.median(unique_y)

    # checks list lengths to ensure matching and then averages the velocities for all files
    # and then returns an array with position and average velocities
    if len(x_pos) == len(y_pos):
        pos_count = len(x_pos)
        if len(x_vel) == len(y_vel) and len(x_vel) == len(z_vel):
            vel_count = len(x_vel)
            nof = vel_count / pos_count  # equals number of files for each height
            ax_vel, ay_vel, az_vel = avg_data_each_h(nof, pos_count, x_vel, y_vel, z_vel)

            if make_sg:
                subgrid_array = sub_grid(ux, x_pos, y_pos, eh, ax_vel, ay_vel, az_vel)
                return subgrid_array
            else:
                z_pos = [eh] * len(x_pos)
                return xmid, ymid, zip(x_pos, y_pos, z_pos, ax_vel, ay_vel, az_vel)
        else:
            print "Error: different number of velocities!"
    else:
        print "Error: not all x-positions have a corresponding y-position!"


def avg_data_each_h(nof, lof, x_vel, y_vel, z_vel):
    """
    Averages the components of the velocities
    :param nof: Number of files for each height
    :param lof: Length of the data set in each file
    :param x_vel: x_velocity
    :param y_vel: y_velocity
    :param z_vel: z_velocity
    :return: Three arrays: averaged x-/y-/z- velocities
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


def sub_grid(unique_x, xpos, ypos, zpos, axvel, ayvel, azvel):
    """
    Creates an n*n subgrid
    :param unique_x: Number of unique x positions
    :param xpos: Array of x_positions
    :param ypos: Array of y_positions
    :param zpos: Array of z_positions
    :param axvel: Array of the average x_velocity corresponding to the positions
    :param ayvel: Array of the average y_velocity corresponding to the positions
    :param azvel: Array of the average z_velocity corresponding to the positions
    :return: An array containing the average parameter value for each subgrid
    """

    n = sgs
    ssgh_array = []

    i = 0
    while i + n + ((n - 1) * unique_x) < len(xpos):
        # Makes sure that subgrid can be fromed
        sxp = 0
        syp = 0
        szp = 0
        sxv = 0
        syv = 0
        szv = 0

        for j in range(n):
            for k in range(n):
                sxp += xpos[i + j + (k * unique_x)]
                syp += ypos[i + j + (k * unique_x)]
                szp += zpos
                sxv += axvel[i + j + (k * unique_x)]
                syv += ayvel[i + j + (k * unique_x)]
                szv += azvel[i + j + (k * unique_x)]
        ssgh_array.append([sxp, syp, szp, sxv, syv, szv])

        if (i + n) < len(xpos):
            i += n
        else:
            pl = unique_x - (i % unique_x)
            i += pl + ((n - 1) * unique_x)

    return np.array(ssgh_array) / sqr(n)


def dict_to_array(dict_array):
    """
    Transforms a dictionary to an array
    :param dict_array: Dictionary
    :return: Array
    """
    plottable_array = []
    for k in dict_array:
        for i in range(len(dict_array[k])):
            plottable_array.append(dict_array[k][i])

    return np.array(plottable_array)


def plot_3d_vector(pa):
    """
    Plots 3D vector graph
    :param pa: Array containing [x_position, y_pos, z_pos x_velocity, y_vel, z_vel]
    :return: n/a
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


def plot_2d_vector(dicti):
    """
    Produces all the 2D plots from the given dictionary. These plots are the 2D vector
    plots and the mean velocity for each height.
    :param dicti: Dictionary; for each height there is a corresponding array of
                    [x_position, y_pos, z_pos x_velocity, y_vel, z_vel]
    :return: n/a
    """

    hcount = 0
    for k in dicti:
        pa2d = []
        for i in range(len(dicti[k])):
            if dicti[k][i][2] == exp_h_list[hcount]:
                pa2d.append(
                    [dicti[k][i][0], dicti[k][i][1], dicti[k][i][3],
                     dicti[k][i][4]])

        X, Y, U, V = zip(*pa2d)
        A = np.sqrt(np.power(X, 2.0) + np.power(Y, 2.0))
        fig = plt.quiver(X[::peo2], Y[::peo2], U[::peo2], V[::peo2], A)
        plt.colorbar(fig)
        plt.title(
            r"$\mu$-PIV vector plot at height %.3f, %s, %s" % (exp_h_list[hcount], shark_species,
                                                               sample_area))
        plt.xlabel(
            r"Average velocity: (%.3f $\bar{x}$ + %.3f $\bar{y}$) $ms^{-1}$" % (
            np.mean(U), np.mean(V)))
        plt.show()

        hcount += 1


def plane_plots(xpv, ypv, dicti):
    """
    Makes planar vector plots of the xz- and yz-planes
    :param xpv: x-position for yz-plane
    :param ypv: y-position for xz-plane
    :param dicti: Dictionary; for each height there is a corresponding array of
                    [x_position, y_pos, z_pos x_velocity, y_vel, z_vel]
    :return: n/a
    """
    xz = []
    yz = []

    for k in dicti:
        for i in range(len(dicti[k])):
            if dicti[k][i][1] == ypv:
                xz.append([dicti[k][i][0], dicti[k][i][2], dicti[k][i][3], dicti[k][i][5]])
            if dicti[k][i][0] == xpv:
                yz.append([dicti[k][i][1], dicti[k][i][2], dicti[k][i][4], dicti[k][i][5]])

    xzx, xzz, xzxv, xzzv = zip(*xz)
    yzy, yzz, yzyv, yzzv = zip(*yz)

    axv, exv= mean_vel(xzxv, xzz)
    ayv, eyv = mean_vel(yzyv, yzz)

    # Changes height values to distance from dermal denticles.
    dfdd = np.array(abs(exp_h_list) - np.amin(abs(exp_h_list)))
    xzz = np.array(abs(np.array(xzz)) - np.amin(abs(np.array(xzz))))
    yzz = np.array(abs(np.array(yzz)) - np.amin(abs(np.array(yzz))))

    f, axarr = plt.subplots(2,2, sharey=True)

    fig = axarr[0,0].quiver(xzx, xzz, xzxv, xzzv, yzyv)
    # Colour uses yzyv so that it matches the other plot
    axarr[0,0].set_xlabel("x")
    axarr[0,0].set_ylabel(r"Height from dd ($mm$)")

    axarr[0,1].errorbar(100*axv, dfdd, xerr=100*exv, marker='o', color='g')
    axarr[0,1].plot([0.0, 0.0], [np.amin(dfdd), np.amax(dfdd)], 'k--')
    axarr[0,1].set_xlabel(r"x-velocity ($\times 10^{-2}ms^{-1}$)")

    axarr[1,0].quiver(yzy, yzz, yzyv, yzzv, yzyv)
    axarr[1,0].set_xlabel("y")
    axarr[1,0].set_ylabel(r"Height from dd ($mm$)")

    axarr[1,1].errorbar(100*ayv, dfdd, xerr=100*eyv, marker='o', color='b')
    axarr[1,1].plot([0.0, 0.0], [np.amin(dfdd), np.amax(dfdd)], 'k--')
    axarr[1,1].set_xlabel(r"y-velocity ($\times 10^{-2}ms^{-1}$)")

    cax, kw = mpl.colorbar.make_axes([ax for ax in axarr.flat])
    f.colorbar(fig, cax=cax, **kw)
    cbar = mpl.colorbar.ColorbarBase(cax, norm=mpl.colors.Normalize(vmin=-0.1, vmax=0.1))
    cbar.set_clim(-0.1,0.1)
    f.suptitle(r"$\mu$-PIV for the %s (%s). Flow direction: %s" %(shark_species, sample_area,
                                                                  flow_direction))
    plt.show()


def mean_vel(vel_array, z_array):

    ava = []
    eva = []

    for i in range(len(exp_h_list)):
        va = []
        for j in range(len(z_array)):
            if z_array[j] == exp_h_list[i]:
                va.append(vel_array[j])

        ava.append(np.mean(va))
        eva.append(np.std(va))

    if sem_bar:
        return np.array(ava), np.array(eva) / len(eva)
    elif sd_bar:
        return np.array(ava), np.array(eva)



if __name__ == '__main__':

    main_dir = "/home/callumkift/Documents/sharks_dtu/micro_piv/20150709_x10_bonnet_c1_s1/"
    # main dir where all the subdirs with the data are

    shark_species = "Bonnethead"
    sample_area = "C1, Side 1"
    flow_direction = "-ve y."
    exp_h_list = read_hf()  # vertical heights of PIV
    lehl = len(exp_h_list)

    make_sg = False  # True -> makes subgrids
    sgs = 3

    plot3D = False  # True -> plots 3D vector
    peo3 = 10  # plots every nth vector for the 3D plot

    plot2D = False  # True -> plots 2D vector plots for all heights
    peo2 = 2  # plots every nth vector for the 2D plot

    # At least one must be True
    sem_bar = True  # plots standard error on mean bars on 2d_mean_roi graph
    sd_bar = False  # plot standard deviation bars on 2d_mean_roi graph
    # If both true, sem will be plotted

    midx = 0
    midy = 0

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

                    midx, midy, h_pos_vel_dict["height{0}".format(fhc)] = get_data(
                        exp_h_list[
                                                                                      hcount],
                                                                       height_file_dict[k])
                    hcount += 1
                h_pos_vel_dict = collections.OrderedDict(sorted(h_pos_vel_dict.items()))

                if plot3D:
                    pa = dict_to_array(h_pos_vel_dict)
                    plot_3d_vector(pa)

                plane_plots(midx, midy, h_pos_vel_dict)

                if plot2D:
                    plot_2d_vector(h_pos_vel_dict)
            else:
                print "\nError: height list does not match number of subdirectories containing files!"
        else:
            print "\nError: experimental height measurements not given!"
    else:
        print "\nError: please choose which error bars to plot"
