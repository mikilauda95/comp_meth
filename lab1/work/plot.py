import sys
import matplotlib.pyplot as plt
import subprocess
import numpy

#argv[]= "num_point_iter" "test" "configuration" "title"

num_point_iter = int(sys.argv[1])
test= int(sys.argv[2])
configuration= int(sys.argv[3])
title=sys.argv[4]

N=1
# print(SNR)
# set the configuration for the resolution
num_points=[64, 256, 1024, 4096]
# num_points=[64, 256, 1024]
# num_points=[64, 256]
num_points=num_points[:num_point_iter]
SNR=numpy.zeros((4, 80))
for (indexnumpoints,num_point) in enumerate(num_points):
    for count in range(N):
        p=subprocess.Popen(["./fft", str(num_point), str(test), str(configuration)], stdout=subprocess.PIPE,bufsize=1, universal_newlines=True )

        f=str(p.stdout.read()).split("\n")
        f=f[1:]
        f=[x.split(", ") for x in f]
        data=[(float(x[0]), float(x[1])) for x in f[:-1]]
        db=numpy.array([x[0] for x in data])
        SNR[indexnumpoints]+=numpy.array([x[1] for x in data])

SNR/= N

if num_point_iter == 1:
    plt.plot(db, SNR[0],"r")
elif num_point_iter == 2:
    plt.plot(db, SNR[0],"r", db, SNR[1], "b")
elif num_point_iter == 3:
    plt.plot(db, SNR[0],"r", db, SNR[1], "b", db, SNR[2], "g")
elif num_point_iter == 4:
    plt.plot(db, SNR[0],"r", db, SNR[1], "b", db, SNR[2], "g", SNR[3], "y")

plt.xlabel("POWER (dB)")
plt.ylabel("Distortion (dB)")
plt.title(title)

plt.show()
