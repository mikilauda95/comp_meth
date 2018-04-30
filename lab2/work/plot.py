import numpy as np
import matplotlib.pyplot as plt
import sys

f = open(sys.argv[1], "r")
meas_list=f.readlines()
meas_list=[x.split() for x in meas_list]
scalar=[float(x[1]) for x in meas_list if x[0]=="scalar"]
vector=[float(x[1]) for x in meas_list if x[0]=="vector"]
vector256=[float(x[1]) for x in meas_list if x[0]=="vector256"]
print(len(scalar), len(vector))
print(scalar[-10:-1])
x=np.arange(1,len(scalar)+1)

colors = ["r", "b", "g"]
plt.plot(x, scalar, colors[0], label = "Scalar" )
plt.plot(x,vector, colors[1], label = "SSE4" )
plt.plot(x, vector256, colors[2], label = "AVX2")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
plt.savefig("plots/"+sys.argv[2])
plt.show()
