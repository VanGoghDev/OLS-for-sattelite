from ISZ import *
from Satellite import *
from DormanPrins_lab2 import TDP
from numpy import linalg


class PCA:
    def __init__(self, result):

        """
        :param result: выборка элементов. В данном случае, result модели в интервале времени когда ИСЗ был замечен НИПом
        """
        self.dormanPrins = TDP()  # интегратор
        self.dormanPrins.geps = 1e-8

        self.deltaXYZ = 100  # разброс для координат
        self.deltaV = 1  # разброс для скоростей
        self.count = result.shape[0]
        self.n = result.shape[1] - 1

        self.result = result[0]

        for i in range(3):
            self.result[i + 1] += self.deltaXYZ
            self.result[i + 4] += self.deltaV

        resultForIsz = np.zeros(6)
        for i in range(6):
            resultForIsz[i] = self.result[i + 1]

        self.Sputnik = Satellite()
        self.ISZ = ISZ(0, self.count, 1, self.n, resultForIsz, self.Sputnik)  # опорная траектория
        self.dormanPrins.run(self.ISZ)  # интегрирование опорной траектории

        self.count = self.ISZ.count

        self.resultChangePlus = np.zeros((6, self.n))
        self.resultChangeMinus = np.zeros((6, self.n))

        for j in range(6):
            for k in range(6):
                self.resultChangePlus[j][k] = self.ISZ.result[0][k + 1]
                self.resultChangeMinus[j][k] = self.ISZ.result[0][k + 1]

        for k in range(3):
            self.resultChangePlus[k][k] += self.deltaXYZ
            self.resultChangePlus[k + 3][k + 3] += self.deltaV
            self.resultChangeMinus[k][k] -= self.deltaXYZ
            self.resultChangeMinus[k + 3][k + 3] -= self.deltaV

        self.modelPlus = ISZ(0, self.count, 1, self.n, self.resultChangePlus[0], self.Sputnik)
        self.modelMinus = ISZ(0, self.count, 1, self.n, self.resultChangeMinus[0], self.Sputnik)

        # self.count += self.count

        self.Hm = np.zeros((self.count + self.count, self.n))
        self.D = np.zeros((self.count, self.count))
        self.iscElAz = np.zeros((self.count, 1))  # матрица искаженных значений элевации и азимута

        self.countD()  # считаем матрицу D

    def countD(self):
        disp = 3.3
        deltaD = random.normalvariate(0, disp)
        # row = self.ISZ.Azimut.shape[0]
        # for i in range(row):
        # self.iscAzimut[i] = self.ISZ.Azimut[i] + deltaD
        # self.iscElevation[i] = self.ISZ.Elevation[i] + deltaD
        row = self.iscElAz.shape[0]
        for i in range(row):
            self.iscElAz[i] = self.ISZ.ElevationAzimut[i] + deltaD
        for i in range(self.count):
            for j in range(self.count):
                if i != j:
                    self.D[i][j] = 0
                else:
                    self.D[i][j] = disp

    def countH(self):
        column = 0
        while column < self.n:  # до нужного количства диф. уравнений (или столбцов)
            self.modelPlus = ISZ(0, self.count, 1, self.n, self.resultChangePlus[column], self.Sputnik)
            self.modelMinus = ISZ(0, self.count, 1, self.n, self.resultChangeMinus[column], self.Sputnik)
            self.dormanPrins.run(self.modelPlus)
            self.dormanPrins.run(self.modelMinus)
            
            for row in range(self.count):  # количество измерений

                dxtecdx = np.zeros(self.n)
                dfidxtec = np.zeros(self.n)
                dfidxtec1 = np.zeros(self.n)
                dfidxtec2 = np.zeros(self.n)
                for i in range(3):
                    dxtecdx[i] = (self.modelPlus.result[row][i+1] - self.modelMinus.result[row][i+1]) /\
                                 (2 * self.deltaXYZ)
                    dxtecdx[i+3] = (self.modelPlus.result[row][i+4] - self.modelMinus.result[row][i+4]) /\
                                   (2 * self.deltaV)

                    dfidxtec[i] = self.ISZ.d[row][i]

                    dfidxtec1[i] = self.ISZ.ElevationAzimut[row*2]
                    dfidxtec2[i] = self.ISZ.ElevationAzimut[row*2+1]

                self.Hm[row*2][column] = dxtecdx[0] * dfidxtec1[0] + dxtecdx[1] * dfidxtec1[1] + dxtecdx[2] * dfidxtec1[2]
                self.Hm[row*2+1][column] = dxtecdx[0] * dfidxtec2[0] + dxtecdx[1] * dfidxtec2[1] + dxtecdx[2] * dfidxtec2[2]
            column += 1
        self.countK(0)

    # поиск корреляционной матрицы ошибок оценивания
    def countK(self, k):
        D_1 = linalg.inv(self.D)  # обратная матрица D
        Ht = self.Hm.transpose()  # транспонированная матрица H
        HD_1 = np.dot(Ht, D_1)  # произведение H трансп на обратную матрицу D
        K = np.dot(HD_1, self.Hm)
        K = linalg.inv(K)  # находим матрицу K
        KHt = np.dot(K, Ht)
        KHtD_1 = np.dot(KHt, D_1)
        dif = np.zeros((self.ISZ.Elevation.shape[0], 1))
        for i in range(self.ISZ.Elevation.shape[0]):
            dif[i] = self.iscElAz[i] - self.ISZ.ElevationAzimut[i]
        deltaX = np.dot(KHtD_1, dif)
        X1 = np.zeros((self.n, 1))
        for i in range(self.n):
            X1[i] = self.ISZ.result[k][i] + deltaX[i]
        a = 3
