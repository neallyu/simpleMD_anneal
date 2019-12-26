import matplotlib.pyplot as plt
import numpy as np
import sys, subprocess, time, math


def trendPlot(path, input_filename):
    data = np.loadtxt(path + input_filename)
    potentialValue = data[:, 0]
    gradientSquaredNorm = data[:, 1]

    fig = plt.figure(1, dpi=500, figsize=(9, 4.8), facecolor="white")
    ax = fig.add_subplot(111)
    ax.plot(potentialValue, color="orange", label="Potential")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Value")
    plt.yticks([-44.3270, -44.3265, -44.3260, -44.3255, -44.3250], ["-44.3270", "-44.3265", "-44.3260", "-44.3255", "-44.3250"])
    ax.axis([0, len(data[:, 0]), -44.327, -44.325])
    plt.legend()
    plt.savefig(path + "Potential_" + input_filename[:-4] + ".png")
    plt.clf()

    fig = plt.figure(1, dpi=500, figsize=(7.4, 4.8), facecolor="white")
    ax = fig.add_subplot(111)
    ax.plot(gradientSquaredNorm, color="green", label="Gradient Squared Norm", linestyle=":")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Value")
    plt.legend()
    plt.savefig(path + "Gradient Squared Norm_" + input_filename[:-4] + ".png")
    plt.clf()


def mk_plot(path, filename, plot_function):
    print("[MD Script]", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), "Plotting", filename)
    plot_function(path, filename)
    print("[MD Script]", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), "Saved to", filename[:-4] + ".png")


if __name__ == "__main__":
    # path = __file__[:-18] + "../output3/CG1/"
    # executableCmd = __file__[:-18] + "../src/CG_optimization " + path + "coordinates_rawdata.txt 10 13 " + path
    # subprocess.call(executableCmd, shell=True)
    path = __file__[:-18] + sys.argv[1] + "/"
    mk_plot(path, "trend.csv", trendPlot)