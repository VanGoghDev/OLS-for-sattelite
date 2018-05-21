import numpy as np
import DormanPrins_lab2 as dp
import matplotlib.pyplot as plt
import matplotlib as mpl
import Satellite
import PrincipalComponentsAnalysis
from ISZ import ISZ


def build_plot(m: ISZ):
    rows = m.result.shape[0]
    axis1 = np.zeros(rows)
    axis2 = np.zeros(rows)
    for i in range(rows):
        axis1[i] = m.result[i][1]
        axis2[i] = m.result[i][2]
    mpl.pyplot.scatter(axis1, axis2, marker='.', linewidths=1, label=r'$\Check x$')
    plt.grid(True)
    mpl.pyplot.show()


def build_plot2(m: ISZ):
    rows = m.ElevationAzimut.shape[0]
    axis2 = np.zeros(rows)
    t = m.count
    axis1 = np.empty(rows)
    for i in range(rows):
        axis1[i] = i
    for i in range(rows):
        axis2[i] = m.ElevationAzimut[i]
    mpl.pyplot.scatter(axis1, axis2, marker='.', linewidths=1, label=r'$\Check x$')
    plt.grid(True)
    mpl.pyplot.show()


def main():
    dorm__prins = dp.TDP()
    dorm__prins.geps = 1e-8
    Sputnik = Satellite.Satellite()
    obj = ISZ(0, 17000, 1, 6, Sputnik.x0, Sputnik)
    dorm__prins.run(obj)
    res = obj.OpornResult
    count = obj.count
    pca = PrincipalComponentsAnalysis.PCA(res, count, 6, obj)
    pca.countH()
    build_plot2(obj)


main()