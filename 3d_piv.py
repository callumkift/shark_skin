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
	a_vel = []

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
						y_pos.append(float(column[1]))
						x_vel.append(float(column[2]))
						y_vel.append(float(column[3]))
						a_vel.append(math.sqrt(sqr(float(column[2])) + sqr(float(column[3]))))
					else:
						x_vel.append(float(column[2]))
						y_vel.append(float(column[3]))
						a_vel.append(math.sqrt(sqr(float(column[2])) + sqr(float(column[3]))))
				else:
					print "Error: TXT file is not correct!"

	# checks list lengths to ensure matching and then averages the velocities for all files
	# and then returns an array with position and average velocities
	if len(x_pos) == len(y_pos):
		pos_count = len(x_pos)
		if len(x_vel) == len(y_vel) and len(x_vel) == len(a_vel):
			vel_count = len(x_vel)
			nof = vel_count/pos_count # equals number of files
			ax_vel, ay_vel, aa_vel = avg_data_each_h(nof, pos_count, x_vel, y_vel, a_vel)

			if len(ax_vel) == len(x_pos):
				height_array = []
				for i in range(len(x_pos)):
					row = []
					row.append(x_pos[i])
					row.append(y_pos[i])
					row.append(eh)
					row.append(ax_vel[i])
					row.append(ay_vel[i])
					row.append(aa_vel[i])
					height_array.append(np.array(row))

				height_array = np.array(height_array)
				return height_array
			else:
				print "Error: averaged velocities do not match with position data!"

		else:
			print "Error: different number of velocities!"
	else:
		print "Error: not all x-positions have a corresponding y-position!"

	



def avg_data_each_h(nof, lof, x_vel, y_vel, a_vel):
	"""
	Averages the x-/y-/abs-velocities for each position at a certain height.
	"""
	sx_vel = []
	sy_vel = []
	sa_vel = []

	for i in range(lof):
		sx = 0
		sy = 0
		sa = 0
		for j in range (nof):
			sx += x_vel[i + (j*nof)]
			sy += y_vel[i + (j*nof)]
			sa += a_vel[i + (j*nof)]

		sx_vel.append(sx)
		sy_vel.append(sy)
		sa_vel.append(sa)

	# checks lengths match and then averages them and returns the average velocities
	if len(sx_vel) == len(sy_vel) and len(sx_vel) == len(sa_vel):
		ax_vel = np.array(sx_vel)/nof
		ay_vel = np.array(sy_vel)/nof
		aa_vel = np.array(sa_vel)/nof

		if len(ax_vel) == len(ay_vel) and len(ax_vel) == len(aa_vel):
			return ax_vel, ay_vel, aa_vel
		else:
			print "Error: averaged velocity data not matching!"
	else:
		print "Error: summed velocity data not matching!"


if __name__ == '__main__':

	exp_h_list = []

	if len(exp_h_list) != 0:

		main_dir = "/home/callumkift/Documents/sharks_dtu/micro_piv/20150609_x20_bdft_os/" #main dir where all the subdirs with the data are
		sub_dirs = get_subdirectories(main_dir)
		
		# Store files of each subdir in a dictionary
		height_file_dict = {}
		for dir in sub_dirs:
			height_file_dict["height{0}".format(dir)]=get_files(dir)
		height_file_dict = collections.OrderedDict(sorted(height_file_dict.items())) #sorts dictionary by subdir

		if len(height_file_dict) == len(exp_h_list):
			# Stores an array filled with position and velocities for each height
			h_pos_vel_dict = {}
			hcount = 1
			for k in height_file_dict:
				h_pos_vel_dict["height{0}".format(hcount)]=get_data(exp_h_list[hcount-1], height_file_dict[k])
				hcount += 1
			h_pos_vel_dict = collections.OrderedDict(sorted(h_pos_vel_dict.items())) 

			for k in h_pos_vel_dict:
				print k, h_pos_vel_dict[k][0]
		else:
			print "Error: height list does not match number of subdirectories!"

	else:
		print "Error: experimental height measurements not given!"

