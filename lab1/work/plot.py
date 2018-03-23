import sys

import matplotlib.pyplot as plt
import subprocess
import numpy

#argv[]= "num_point_iter" "test" "configuration" "mc iterations" "title"

if sys.argv[1]=="all":
    #all the possible configurations
    simulation_list=[
        [[64, 256, 1024, 4096], 0, 0, 0, 10],
        [[64, 256, 1024, 4096], 0, 1, 0, 10],
        [[64, 256, 1024, 4096], 2, 0, 0, 10],
        [[64, 256, 1024, 4096], 2, 1, 0, 10],
        [[64, 256, 1024, 4096], 3, 0, 0, 10],
        [[64, 256, 1024, 4096], 3, 1, 0, 10]
    ]
elif sys.argv[1]=="rounding":
    simulation_list=[
        [[1024], 0, 0, 0, 10],
        [[1024], 0, 0, 1, 10]
    ]
else:
    simulation_list=[
        [[int(sys.argv[1])], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])]
    ]

d_sign={
    0 : "Cosine",
    1 : "QPSK",
    2 : "16QAM",
    3 : "White Gaussian Noise"
}

d_repr={
    0 : "Q15",
    1 : "Q24"
}


d_round={
    0 : "Twiddle Floored",
    1 : "Twiddle Rounded"
}

for ind_sim, simulation in enumerate(simulation_list):
    [num_points, test, configuration, rounding, N] = simulation
    print ( simulation )

    title=d_sign[test]+" - "+d_repr[configuration]
    if sys.argv[1]=="rounding":
        title += " - " + d_round[rounding]

    SNR=numpy.zeros((4, 40))
    for (indexnumpoints,num_point) in enumerate(num_points):
        print("running num_points = ", num_point)
        for count in range(N):
            p=subprocess.Popen(
                               ["./fft", str(num_point), str(test), str(configuration), str(rounding)],
                               stdout=subprocess.PIPE,bufsize=1, universal_newlines=True
                              )
            f=str(p.stdout.read()).split("\n")
            f=f[1:]
            f=[x.split(", ") for x in f]
            data=[(float(x[0]), float(x[1])) for x in f[:-1]]
            db=numpy.array([x[0] for x in data])
            SNR[indexnumpoints]+=numpy.array([x[1] for x in data])

    SNR/= N

    labels=["N = "+ str(x) for x in num_points]
    print(labels)
    #for plotting the comparison on the two rounding
    if sys.argv[1] == "rounding":
        if ind_sim==0:
            plt.plot(db, SNR[0],"r", label=labels[0])
        else:
            plt.plot(db, SNR[0],"b", label=labels[0])
    #Plot what you need
    elif len(num_points) == 1:
        plt.plot(db, SNR[0],"r", label=labels[0])
    elif len(num_points) == 2:
        plt.plot(db, SNR[0],"r", label=labels[0])
        plt.plot(db, SNR[1],"b", label=labels[1])
    elif len(num_points) == 3:
        plt.plot(db, SNR[0],"r", label=labels[0])
        plt.plot(db, SNR[1],"b", label=labels[1])
        plt.plot(db, SNR[2],"g", label=labels[2])
    elif len(num_points) == 4:
        plt.plot(db, SNR[0],"r", label=labels[0])
        plt.plot(db, SNR[1],"b", label=labels[1])
        plt.plot(db, SNR[2],"g", label=labels[2])
        plt.plot(db, SNR[3],"y", label=labels[3])

    #Write the labels, title, and axes
    plt.xlabel("Input Power (dB)")
    plt.ylabel("SNR (dB)")
    plt.title(title)
    plt.legend()

    #Write the filenames and save
    filename=d_sign[test].replace(" ", "_")+"_"+d_repr[configuration].replace(" ", "_")
    if sys.argv[1]=="rounding":
        filename=d_sign[test].replace(" ", "_")+"_"+d_repr[configuration].replace(" ", "_")+"round_compare"
    else:
        plt.savefig("plots/"+filename+".png")
        plt.close()
#for plotting the two configurations with round and floor
if sys.argv[1]=="rounding":
    plt.savefig("plots/"+filename+".png")
    plt.close()
    # plt.show()
