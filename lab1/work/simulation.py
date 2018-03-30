import sys

import matplotlib.pyplot as plt
import subprocess
import numpy

#argv[]= "num_point_iter" "test" "configuration" "mc iterations" "title"

fi=open("ranges", 'w')

if sys.argv[1]=="all":
    #all the possible configurations
    simulation_list=[
        [[64, 256, 1024, 4096], 0, 0, 0, 1],
        [[64, 256, 1024, 4096], 0, 1, 0, 1],
        [[64, 256, 1024, 4096], 2, 0, 0, 10],
        [[64, 256, 1024, 4096], 2, 1, 0, 10],
        [[64, 256, 1024, 4096], 3, 0, 0, 10],
        [[64, 256, 1024, 4096], 3, 1, 0, 10]
    ]
elif sys.argv[1]=="test":
    #Testing values
    simulation_list=[
        [[64, 256], 0, 0, 0, 1],
        [[64, 256], 0, 1, 0, 1],
        [[64, 256], 2, 0, 0, 1],
        [[64, 256], 2, 1, 0, 1],
        [[64, 256], 3, 0, 0, 1],
        [[64, 256], 3, 1, 0, 1]
    ]
elif sys.argv[1]=="rounding":
    simulation_list=[
        [[1024], 0, 0, 0, 1],
        [[1024], 0, 0, 1, 1]
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

plot_points=200

for ind_sim, simulation in enumerate(simulation_list):
    [num_points, test, configuration, rounding, N] = simulation
    print ( simulation )

    title=d_sign[test]+" - "+d_repr[configuration]
    filename=(d_sign[test]+"_"+d_repr[configuration]).replace(" ", "_")
    points_file=filename + "_points"
    base="plots/"
    fpoints = open(base+points_file, "w")
    if sys.argv[1]=="rounding":
        title += " - " + d_round[rounding]

    SNR = numpy.empty((0,plot_points), float)
    for (indexnumpoints,num_point) in enumerate(num_points):
        print(num_point)
        SNR_tmp=numpy.zeros((1,plot_points))
        db = numpy.empty((1,plot_points), float)
        for count in range(N):
            print("I am at the iteration ", str(N))
            p=subprocess.Popen(
                               ["./fft", str(num_point), str(test), str(configuration), str(rounding)],
                               stdout=subprocess.PIPE,bufsize=1, universal_newlines=True
                              )
            f=str(p.stdout.read()).split("\n")
            f=f[1:]
            f=[x.split(", ") for x in f]
            data=[(float(x[0]), float(x[1])) for x in f[:-1]]
            db=numpy.array([[x[0] for x in data]])
            SNR_tmp+=numpy.array([x[1] for x in data])
        SNR_tmp /= N
        SNR_tmp=numpy.append(SNR_tmp, db, axis=0)
        SNR=numpy.append(SNR, SNR_tmp, axis=0)
    numpy.save(fpoints, SNR)
