import sys

import matplotlib.pyplot as plt
import subprocess
import numpy

# Load the points

filename=sys.argv[1]

data=numpy.load(filename)

SNR=data[0:len(data):2]
db=data[1:len(data):2]
print(len(data))
print(SNR)
print(db)

# Compute the cut off ranges

title = sys.argv[1]
num_points = [64, 256, 1024, 4096]

cut1=numpy.zeros((4, 2))
cut2=numpy.zeros((4, 2))

for indexnumpoints in range(4):
    flag=0
    for x in zip(SNR[indexnumpoints], db[indexnumpoints]):
        if flag==0:
            if x[0] > 50 and x[1] < -20:
                cut1[indexnumpoints][0]=x[0]
                cut1[indexnumpoints][1]=x[1]
                flag=1
        elif flag==1:
            if x[0] < 50 and x[1] > -20:
                cut2[indexnumpoints][0]=x[0]
                cut2[indexnumpoints][1]=x[1]
                flag=0

labels=["N = "+ str(x) for x in num_points]
print(labels)
plt.subplot(111)

#for plotting the comparison on the two rounding
if sys.argv[1] == "rounding":
    if ind_sim==0:
        plt.plot(db, SNR[0],"r", label=labels[0])
    else:
        plt.plot(db, SNR[0],"b", label=labels[0])
#Plot what you need
plt.plot(db[0], SNR[0],"r", label=labels[0])
plt.plot(db[1], SNR[1],"b", label=labels[1])
plt.plot(db[2], SNR[2],"g", label=labels[2])
plt.plot(db[3], SNR[3],"y", label=labels[3])


th=[50]*len(db[0])
plt.plot(db[0], th, "orange", label = "Threshold")


# plot the cutoff points
# plt.scatter(cut1[0][1], cut1[0][0],color="r")

# if cut2[0][1] and cut2[0][0]:
    # plt.scatter(cut2[0][1], cut2[0][0],color="r")

# plt.scatter(cut1[1][1], cut1[1][0],color="b")

# if cut2[1][1] and cut2[1][0]:
    # plt.scatter(cut2[1][1], cut2[1][0],color="b")

# plt.scatter(cut1[2][1], cut1[2][0],color="g")

# if cut2[2][1] and cut2[2][0]:
    # plt.scatter(cut2[2][1], cut2[2][0],color="b")

# plt.scatter(cut1[3][1], cut1[3][0],color="y")

# if cut2[3][1] and cut2[3][0]:
    # plt.scatter(cut2[3][1], cut2[3][0],color="b")

#Write the labels, title, and axes
plt.xlabel("Input Power (dB)")
plt.ylabel("SNR (dB)")
# plt.title(title)
# plt.legend(bbox_to_anchor=(0.3, 0.2))
# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)

print(cut1, cut2)

plt.show()
# plt.savefig(filename+".png")
# plt.close()

#Write the filenames and save
# if sys.argv[1]=="rounding":
    # filename=(d_sign[test]+"_"+d_repr[configuration]+"round_compare").replace(" ", "_")
# else:
#for plotting the two configurations with round and floor
# if sys.argv[1]=="rounding":
    # plt.savefig("plots/"+filename+".png")
    # plt.close()
    # plt.show()
