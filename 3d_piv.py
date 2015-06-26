#!/usr/bin/env python

###
#
# This script creates a 3D vector map using data from
# the micro-PIV experiments. It averages over each hight
# and combines this into a single graph
#
###

import numpy as np
import os
import collections
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mplc

def sqr(a):
	"""
	Returns the square of the input.
	"""
	return a*a

def get_subdirectories(a_dir):
	"""
	The input is a direcotry. Retrieves the sub-directories of the folder. Returns a list of
	subdirectories.
	"""
	return [a_dir+name+"/" for name in os.listdir(a_dir)
			if os.path.isdir(os.path.join(a_dir, name))]

def get_files(a_dir):
	"""
	The input is a directory. Retrieves all txt files from the directory. Returns a list of
	files
	"""
	gf = []
	for file in os.listdir(dir):
		if file.endswith(".txt"):
			gf.append(dir + "/" + str(file))
	if len(gf) != 0:
		return gf
	else:
		print "Error: Cannot find TXT files in subdirectories!"

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
	unique_x = []
	unique_y = []

	# reading data
	for file in file_list:
		with open (file, 'r') as f:
			header1 = f.readline()
			for line in f:
				line = line.strip()
				column = line.split()
				if len(column) == 4:
					if file == file_list[0]:
						# Only takes position data from first file as the same in each file
						x_pos.append(float(column[0]))
						if float(column[0]) not in unique_x:
							unique_x.append(float(column[0]))
						y_pos.append(float(column[1]))
						if float(column[1]) not in unique_y:
							unique_y.append(float(column[1]))
						x_vel.append(float(column[2]))
						y_vel.append(float(column[3]))
						# z_vel.append(math.sqrt(sqr(float(column[2])) + sqr(float(column[3]))))
						z_vel.append(0.0)
					else:
						x_vel.append(float(column[2]))
						y_vel.append(float(column[3]))
						# z_vel.append(math.sqrt(sqr(float(column[2])) + sqr(float(column[3]))))
						z_vel.append(0.0)
				else:
					print "Error: TXT file is not correct!"

	ux = len(unique_x)
	uy = len(unique_y)

	# checks list lengths to ensure matching and then averages the velocities for all files
	# and then returns an array with position and average velocities
	if len(x_pos) == len(y_pos):
		pos_count = len(x_pos)
		if len(x_vel) == len(y_vel) and len(x_vel) == len(z_vel):
			vel_count = len(x_vel)
			nof = vel_count/pos_count # equals number of files
			ax_vel, ay_vel, az_vel = avg_data_each_h(nof, pos_count, x_vel, y_vel, z_vel)

			if len(ax_vel) == len(x_pos):
				subgrid_array = sub_grid(ux, uy, x_pos, y_pos, eh, ax_vel, ay_vel, az_vel)
				return subgrid_array
			else:
				print "Error: averaged velocities do not match with position data!"

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
		for j in range (nof):
			sx += x_vel[i + (j*nof)]
			sy += y_vel[i + (j*nof)]
			sa += z_vel[i + (j*nof)]

		sx_vel.append(sx)
		sy_vel.append(sy)
		sz_vel.append(sa)

	# checks lengths match and then averages them and returns the average velocities
	if len(sx_vel) == len(sy_vel) and len(sx_vel) == len(sz_vel):
		ax_vel = np.array(sx_vel)/nof
		ay_vel = np.array(sy_vel)/nof
		az_vel = np.array(sz_vel)/nof

		if len(ax_vel) == len(ay_vel) and len(ax_vel) == len(az_vel):
			return ax_vel, ay_vel, az_vel
		else:
			print "Error: averaged velocity data not matching!"
	else:
		print "Error: summed velocity data not matching!"

def sub_grid(unique_x, unique_y, xpos, ypos, zpos, axvel, ayvel, azvel):
	"""
	Reduces data by subdividing our roi into n*n grids, with each grid containing the average
	of the n*n velocities and positions.
	"""
	print unique_x, unique_y, len(xpos), len(axvel)

	n = sgs
	ssgh_array = []

	i = 0
	while i < (len(xpos)):
		sxp = 0
		syp = 0
		szp = 0
		sxv = 0
		syv = 0
		szv = 0
		if (i + n) < len(xpos) and (i + n + (n-1)*unique_y) < len(xpos):
			for j in range(n):
				for k in range(n):
					sxp += xpos[i + j + (k*unique_y)]/sqr(n)
					syp += ypos[i + j + (k*unique_y)]/sqr(n)
					szp += zpos/sqr(n)
					sxv += axvel[i + j + (k*unique_y)]/sqr(n)
					syv += ayvel[i + j + (k*unique_y)]/sqr(n)
					szv += azvel[i + j + (k*unique_y)]/sqr(n)
			ssgh_array.append([sxp, syp, szp, sxv, syv, szv])
		
		if (i+n) < len(xpos): 
				i += n
		else:
			pl = unique_y - (i % unique_y)
			i += pl + ((n-1)*unique_y)

	print np.array(ssgh_array).shape
	return np.array(ssgh_array)

def arrays_to_plot(dict_array):
	"""
	Puts the dictionary-filled-array into plottable array.
	"""
	plottable_array = []
	for k in dict_array:
		for i in range(len(dict_array[k])):
			plottable_array.append(dict_array[k][i])

	return np.array(plottable_array)

def plot_vectors(pa):
	"""
	Plots a 3d vector graph
	"""
	X, Y, Z, U, V, W = zip(*pa)
	A = np.sqrt(np.power(X,2) + np.power(Y,2))

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	q = ax.quiver(X, Y, Z, U, V, W, A, length=0.005, arrow_length_ratio=0.9)
	q.set_array(np.random.rand(100))
	plt.colorbar(q)
	ax.w_xaxis.set_pane_color((0.3,0.3,0.3,0.99))
	ax.w_yaxis.set_pane_color((0.3,0.3,0.3,0.99))
	ax.w_zaxis.set_pane_color((0.3,0.3,0.3,0.99))

	plt.show()

if __name__ == '__main__':

	exp_h_list = [10.0, 15.0, 30.0, 40.0, 45.0, 60.0, 70.0, 80.0] # vertical heights of PIV
	lehl = len(exp_h_list)
	sgs = 3 # sub grid size

	if lehl != 0:

		main_dir = "/home/callumkift/Documents/sharks_dtu/micro_piv/20150609_x20_bdft_os/" #main dir where all the subdirs with the data are
		sub_dirs = get_subdirectories(main_dir)
		
		# Store files of each subdir in a dictionary
		height_file_dict = {}
		for dir in sub_dirs:
			height_file_dict["height{0}".format(dir)]=get_files(dir)
		height_file_dict = collections.OrderedDict(sorted(height_file_dict.items())) #sorts dictionary by subdir

		# Stores an array filled with position and velocities for each height
		if len(height_file_dict) == lehl:
			h_pos_vel_dict = {}
			hcount = 1
			for k in height_file_dict:
				h_pos_vel_dict["height{0}".format(hcount)]=get_data(exp_h_list[hcount-1], height_file_dict[k])
				hcount += 1
			h_pos_vel_dict = collections.OrderedDict(sorted(h_pos_vel_dict.items())) 

		else:
			print "Error: height list does not match number of subdirectories!"


		# Turn the dictionary into plottable arrayslea
		pa = arrays_to_plot(h_pos_vel_dict)

		plot_vectors(pa)

	else:
		print "Error: experimental height measurements not given!"

