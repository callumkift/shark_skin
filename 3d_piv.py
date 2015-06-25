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
from math import sqrt
import matplotlib.pyplot as plt


def get_subdirectories(a_dir):
	"""
	The input is a direcotry. Retrieves the sub-directories of the folder. Returns a list of
	subdirectories.
	"""
	return [a_dir+name+"/" for name in os.listdir(a_dir)
			if os.path.isdir(os.path.join(a_dir, name))]



if __name__ == '__main__':

	main_dir = "/home/callumkift/Documents/sharks_dtu/micro_piv/20150609_x20_bdft_os/"
	sub_dirs = get_subdirectories(main_dir)
