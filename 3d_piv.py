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
	return gf

def get_data(file_list):
	"""
	The input is a list of files. 
	"""
	x_pos = []
	y_pos = []
	x_vel = []
	y_vel = []
	abs_vel = []

	file_dict = {}

	for i in range(len(file_list)):
		file_dict["f{0}".format(i)]=file_list[i]

	file_dict = collections.OrderedDict(sorted(file_dict.items()))
	for k in file_dict:
		print k, file_dict[k]

	# for file in file_list:
	# 	with open( file, 'r' ) as f : 
	# 		# Read and ignore header lines:
	# 		header1 = f.readline()
	# 		for line in f:
	# 			line = line.strip() #splits the file into columns
	# 			columns = line.split() #splits the line into columns at every space in the line
	# 			if (len(columns) == 4) : 
	# 				x_vel.append(float(columns[2]))
	# 				y_vel.append(float(columns[3])) 
	# 				abs_vel.append(sqrt(sqr(float(columns[2]))+sqr(float(columns[3]))))
	# 			else:
	# 				print "txt file is not correct."

if __name__ == '__main__':

	main_dir = "/home/callumkift/Documents/sharks_dtu/micro_piv/20150609_x20_bdft_os/" #main dir where all the subdirs with the data are
	sub_dirs = get_subdirectories(main_dir)
	
	# Store files of each subdir in a dictionary
	height_file_dict = {}
	for dir in sub_dirs:
		height_file_dict["height{0}".format(dir)]=get_files(dir)

	height_file_dict = collections.OrderedDict(sorted(height_file_dict.items())) #sorts dictionary by subdir

	for k in height_file_dict:
		print k, height_file_dict[k]
		raw_input("")
		get_data(height_file_dict[k])
		raw_input("")

