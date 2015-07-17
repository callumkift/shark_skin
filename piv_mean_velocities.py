#!/usr/bin/env python

###
#
# This script plots the averaged velocities from the micro-PIV experiments.
#
# Author: Callum Kift
#
# Directory setup:
#
#
###

import os
import collections
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


def get_files(a_dir):
    """
    Creates a list of all txt files in the given directory
    :param a_dir: Path to a directory
    :return: a list of all txt files in the directory
    """
    gf = []
    for file in os.listdir(a_dir):
        if file.endswith(".txt"):
            gf.append(a_dir + str(file))
    if len(gf) != 0:
        return gf
    else:
        print "Error: Cannot find TXT files in subdirectory!\n\t (%s)" % a_dir


def get_data(v_file):
    """
    Retrieves data from the file sent. Returns an array containing heights and averaged velocities.
    """

    height_values = []
    xveloc_values = []
    yveloc_values = []

    f = open(v_file, "r")
    f.readline()  # Ignores first line
    for line in f:
        line = line.strip()
        column = line.split()
        if len(column) == 3:
            height_values.append(float(column[0]))
            xveloc_values.append(float(column[1]))
            yveloc_values.append(float(column[2]))
        else:
            print "Error: File is not in the correct format."
    f.close()

    return [height_values, xveloc_values, yveloc_values]


def plot_xy():
    """
    Plots 2 graphs, showing the averaged x-/y-velocities for each sample.
    """

    f, axarr = plt.subplots(2, sharex=True)
    for key, items in velocity_dict.iteritems():
        xv = 100 * np.array(items[1])
        yv = 100 * np.array(items[2])
        h = items[0]
        axarr[0].plot(xv, h, "*-", label=str(key))
        axarr[1].plot(yv, h, "*-")

    axarr[0].set_xlabel(r"x-velocity ($\times 10^{-2}ms^{-1}$)")
    axarr[0].set_ylabel(r"Height from dd ($mm$)")

    axarr[1].set_xlabel(r"y-velocity ($\times 10^{-2}ms^{-1}$)")
    axarr[1].set_ylabel(r"Height from dd ($mm$)")

    f.suptitle("Averaged velocities")
    axarr[0].legend(loc="upper left", borderaxespad=0.0)
    plt.show()
    return


def plot_raw():
    width_shrink = 0.85
    leg_pos = (1.15, 0.5)
    ax = plt.subplot(111)

    for key, items in velocity_dict.iteritems():
        ax.plot(-100 * np.array(items[2]), items[0], "o-", label=
        string.replace(key[4:], "_x10", ""))

    ax.set_xlabel(r"Flow-velocity ($\times 10^{-2}ms^{-1}$)")
    ax.set_ylabel(r"Height")
    plt.suptitle("Average velocity in direction of flow. Raw data.")

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * width_shrink, box.height])

    fontP = FontProperties()
    fontP.set_size('small')
    ax.legend(prop=fontP, loc=9, bbox_to_anchor=leg_pos)
    plt.show()

    return


def plot_flow():
    """
    Plots the average y-velocity for each sample.
    """
    max_vels = []
    width_shrink = 0.85
    leg_pos = (1.15, 0.5)

    ax = plt.subplot(111)

    for key, items in velocity_dict.iteritems():
        ch = 0.5
        h = np.amax(items[0])
        vel_factor = h / ch

        ax.plot(-100 * vel_factor * np.array(items[2]), np.array(items[0]) * (0.5 / np.amax(items[
                                                                                              0])),
                "o-", label=
        string.replace(key[4:], "_x10", ""))
        max_vels.append(np.amax(-100 * vel_factor * np.array(items[2])))

    max_value = np.amax(max_vels)

    ax.plot((0, max_value), (0.5, 0.5), 'k--')
    ax.set_xlabel(r"Normalised flow-velocity ($\times 10^{-2}ms^{-1}$)")
    ax.set_ylabel(r"Normalised height")
    plt.suptitle("Average velocity in direction of flow. Normalised height, normalised velocity.")

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * width_shrink, box.height])

    fontP = FontProperties()
    fontP.set_size('small')
    ax.legend(prop=fontP, loc=9, bbox_to_anchor=leg_pos)
    plt.show()
    return


if __name__ == '__main__':

    main_dir = "/home/callumkift/Documents/sharks_dtu/micro_piv/mean_vels/"

    # Choose which plots to show
    plotXY = False  # Plots both x- and y-velocities in subplots
    plotFLOW = True  # Plots y-velocities (normally flow direction)

    if not os.path.exists(main_dir):
        print "Error: main directory does not exist!"
    else:

        velocity_files = get_files(main_dir)
        velocity_files.sort()

        velocity_dict = {}

        for file in velocity_files:
            velocity_dict["{}".format(file.rsplit("/", 1)[1][:-4])] = get_data(file)
            # makes the key = filename without .txt, item = array of positions and velocities.

        velocity_dict = collections.OrderedDict(sorted(velocity_dict.items()))
        # Sorts the dictionary

        plot_raw()

        if plotXY:
            plot_xy()
        if plotFLOW:
            plot_flow()
