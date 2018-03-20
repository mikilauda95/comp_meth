import sys
import matplotlib.pyplot as plt
import subprocess
import numpy

N=1
# print(SNR)
# set the configuration for the resolution
configuration=0
num_points=[64, 256, 1024, 4096]
SNR=numpy.zeros((4, 40))
for test in range(1):
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

plt.plot(db, SNR[0],"r", db, SNR[1], "b", db, SNR[2], "g", SNR[3], "g")
plt.xlabel("POWER (dB)")
plt.ylabel("Distortion (dB)")
plt.title("TEST0")
    # if test==0:
        # plt.subplot(221)
        # plt.plot(db,SNR[0], 1, 1)
        # plt.plot(db,SNR[1], 1, 1)
        # plt.xlabel("POWER (dB)")
        # plt.ylabel("Distortion (dB)")
        # plt.title("TEST0")
        # pass
    # elif test==1:
        # plt.subplot(222)
        # plt.plot(db,SNR, 1, 1)
        # plt.xlabel("POWER (dB)")
        # plt.ylabel("Distortion (dB)")
        # plt.title("TEST1")
    # elif test==2:
        # print("test 2")
        # plt.subplot(223)
        # plt.plot(db,SNR, 1, 1)
        # plt.xlabel("POWER (dB)")
        # plt.ylabel("Distortion (dB)")
        # plt.title("TEST2")

plt.show()
