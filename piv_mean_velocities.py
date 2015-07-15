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
import numpy as np
import matplotlib.pyplot as plt


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
        xv = 100*np.array(items[1])
        yv = 100*np.array(items[2])
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


def plot_flow():
    """
    Plots the average y-velocity for each sample.
    """
    for key, items in velocity_dict.iteritems():
        plt.plot(-100*np.array(items[2]), items[0]/np.amax(items[0]), "o-", label = str(key))

    plt.plot((0,10), (0.5,0.5), 'k--')
    plt.xlabel(r"Flow-velocity ($\times 10^{-2}ms^{-1}$)")
    plt.ylabel(r"Height from dd")
    plt.title("Average velocity in direction of flow. Normalised height.")
    plt.legend(loc=1)
    plt.show()
    return



if __name__ == '__main__':

    main_dir = "/home/callumkift/Documents/sharks_dtu/micro_piv/mean_vels/"

    # Choose which plots to show
    plotXY = False # Plots both x- and y-velocities in subplots
    plotFLOW = True # Plots y-velocities (normally flow direction)

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

        if plotXY:
            plot_xy()
        if plotFLOW:
            plot_flow()
