import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np


paral = ["scalar", "128", "256"]
filename = sys.argv[1]
mode = sys.argv[2]

for par in paral:
    lab2_launch = ["./lab2", "min", "data.txt", mode, par]
    for i in range(50):
        p=subprocess.Popen(
            lab2_launch,
            stdout=subprocess.PIPE,bufsize=1, universal_newlines=True
        )
        meas_list = p.stdout.read()
        meas_list = meas_list.split("\n")[:-1]
        meas_list=[x.split() for x in meas_list]
        if par == "scalar":
            scalar=[float(x[1]) for x in meas_list if x[0]=="scalar"]
        if par == "128":
            vector=[float(x[1]) for x in meas_list if x[0]=="vector"]
        if par == "256":
            vector256=[float(x[1]) for x in meas_list if x[0]=="vector256"]

        if i == 0:
            if par == "scalar":
                scalar_minimum = scalar
            if par == "128":
                vector_minimum = vector
            if par == "256":
                vector256_minimum = vector256
        else:
            if par == "scalar":
                scalar_minimum = np.minimum(scalar_minimum, scalar)
            if par == "128":
                vector_minimum = np.minimum(vector_minimum, vector)
            if par == "256":
                vector256_minimum = np.minimum(vector256_minimum, vector256)



vector_minimum[0] = vector_minimum[1]
vector256_minimum[0] = vector256_minimum[1]
x=np.arange(len(scalar_minimum))
colors = ["r", "b", "g"]
plt.plot(x, scalar_minimum, colors[0], label = "Scalar" )
plt.plot(x,vector_minimum, colors[1], label = "SSE4" )
plt.plot(x, vector256_minimum, colors[2], label = "AVX2")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
plt.savefig("plots/"+filename)
plt.show()
