import numpy as np
import matplotlib.pyplot as plt
import sys, getopt

argv = sys.argv[1:]

geoname = 'geometry'
solname = 'solution'
path = './'
try:
  	opts, args = getopt.getopt(argv,"p:g:s:",["path=","gfile=","sfile="])
except getopt.GetoptError:
  	print('Error parsing')
  	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print('test.py -p <path> -g <geoname> -o <solname>')
		sys.exit()
	elif opt in ("-p", "--path"):
		path = arg
	elif opt in ("-g", "--gfile"):
		geoname = arg
	elif opt in ("-s", "--sfile"):
		solname = arg
print('Input file is "'+path+geoname+'"')
print('Output file is "'+path+solname+'"')

mat = np.genfromtxt(path+solname,delimiter=",")
print(mat)

plt.figure()
plt.scatter(mat[:,0],mat[:,1])
plt.show()
