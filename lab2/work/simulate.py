import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np


filename = sys.argv[1]
mode = sys.argv[2]
lab2_launch = ["./lab2", "min", "data.txt", mode]

for i in range(50):
    p=subprocess.Popen(
        lab2_launch,
        stdout=subprocess.PIPE,bufsize=1, universal_newlines=True
    )
    meas_list = p.stdout.read()
    meas_list = meas_list.split("\n")[:-1]
    meas_list=[x.split() for x in meas_list]
    scalar=[float(x[1]) for x in meas_list if x[0]=="scalar"]
    vector=[float(x[1]) for x in meas_list if x[0]=="vector"]
    vector256=[float(x[1]) for x in meas_list if x[0]=="vector256"]

    if i == 0:
        scalar_minimum = scalar
        vector_minimum = vector
        vector256_minimum = vector256
    else:
        scalar_minimum = np.minimum(scalar_minimum, scalar)
        vector_minimum = np.minimum(vector_minimum, vector)
        vector256_minimum = np.minimum(vector256_minimum, vector256)



vector_minimum[0] = vector_minimum[1]
x=np.arange(len(scalar_minimum))
colors = ["r", "b", "g"]
plt.plot(x, scalar_minimum, colors[0], label = "Scalar" )
plt.plot(x,vector_minimum, colors[1], label = "SSE4" )
plt.plot(x, vector256_minimum, colors[2], label = "AVX2")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
plt.savefig("plots/"+filename)
plt.show()
