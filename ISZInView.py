from ISZ import ISZ
import model as model
import numpy as np
import math
import random
import Satellite
from numpy.linalg import norm


class ISZInView(ISZ):
    import model as model
    import numpy as np
    import math
    import random
    import Satellite
    from numpy.linalg import norm
    mu = 398600.436  # km/s^2
    re = 6371  # Earth radius, km
    omega = 7.292115e-5  # angular velocity of earth rotation рад/c
    # initializing starting conditions in constructor
    # we have a system of 6 equations for Vx, Vy, Vz and Vx/dx, Vy/dy, Vz/dz
    def __init__(self, t0, t, h, n, x0, startingConditions: Satellite.Satellite):  # n = 7
        super().__init__(t0, t, h, n, x0, startingConditions)
        self.m = 50  # mass of satellite kg
        self.matrix = np.array([[0, 1.225, -0.2639 * 10e-8, 0.7825e-4],
                                [20000, 0.891e-1, 0.4407e-9, 0.16375e-3],
                                [60000, 2.578e-4, -0.2560e-8, 0.5905e-4],
                                [100000, 4.061e-7, 0.1469e-8, 0.1787e-3],
                                [150000, 2.130e-9, 0.8004e-10, 0.3734e-4],
                                [300000, 4.764e-11, 0.7111e-11, 0.1547e-4],
                                [600000, 8.726e-13, 0.1831e-11, 0.9280e-5],
                                [900000, 6.367e-14, 0, 0.9540e-5]])
        self.omega = np.array([0, 0, 7.292115e-5])
        self.v_alpha = np.zeros(3)
        # x0 = startingConditions.x0
        self.x0 = np.zeros(n)
        for i in range(6):
            self.x0[i] = x0[i]
        # self.x0[6] = 0
        self.count = 0  # счетчик того, сколько раз НИП засечет спутник
        self.OpornResult = np.zeros((0, self.n + 1))  # матрица результатов для опорной траектории
        self.ElevationAzimut = np.zeros((0, 1))  # матрица, которая содержит в себе элевацию и азимут
        self.Elevation = np.zeros((0, 1))
        self.Azimut = np.zeros((0, 1))
        self.iscElAz = np.zeros((0, 1))
        self.e = startingConditions.e
        self.angles = np.zeros((0, 1))  # массив углов, которые измеряются пока ИСЗ в зоне видимости НИПа
        self.d = np.zeros((0, 3))
        self.h = 0
        self.rn = np.zeros((0, 3))
        self.dt = np.zeros((0, 3))
    @staticmethod
    def getsize(a):
        size = 0
        for k in range(len(a)):
            size += 1
        return size
    @staticmethod
    def return_temp(i):
        if i < 900000:
            temp = 0
        elif i < 60000:
            temp = 1
        elif i < 100000:
            temp = 2
        elif i < 150000:
            temp = 3
        elif i < 300000:
            temp = 4
        elif i < 600000:
            temp = 5
        elif i < 900000:
            temp = 6
        elif i >= 20000:
            temp = 7
        else:
            temp = 'Error'
        return temp
    @staticmethod
    def Geographical(a):
        result = np.zeros(3)
        result[0] = math.atan(a[1] / a[0])
        result[1] = math.atan(a[2] / (math.sqrt(pow(a[0], 2) + pow(a[1], 2))))
        result[2] = math.sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2))
        return result
    def std_density(self, height):  # ro(h)
        height *= 1000
        i = self.return_temp(height)
        result = self.matrix[i][1] * math.pow(math.e,
                                              (self.matrix[i][2] * math.pow((height - self.matrix[i][0]), 2) -
                                               self.matrix[i][3] * (height - self.matrix[i][0])))
        result /= 1000
        return result
    def aerodynamic_force(self, a):  # R
        module_x = math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])  # length of a vector
        height = module_x - ISZ.re
        a1 = np.array([a[0], a[1], a[2]])
        c = np.zeros(3)
        for i in range(3):
            c[i] = 0
            for j in range(3):
                c[i] += self.omega[i] * a1[j]
        # c = np.dot(self.omega, a1)
        for i in range(3):
            self.v_alpha[i] = a[i + 3] - c[i]
        module_v_alpha = math.sqrt(self.v_alpha[0] * self.v_alpha[0] + self.v_alpha[1] * self.v_alpha[1]
                                   + self.v_alpha[2] * self.v_alpha[2])
        CxS = 1.4  # m^2
        result = -CxS * 0.5 * self.std_density(height) * module_v_alpha
        result = np.dot(result, self.v_alpha)
        result = result / self.m
        return result
    def add_result(self, a, t):
        k = t / self.sampling_increment
        k = int(math.floor(k))
        self.result[k][0] = t
        for i in range(1, self.getsize(self.x0) + 1):
            self.result[k][i] = a[i - 1]
        rows = self.getsize(self.x0)
        for i in range(1, self.getsize(self.x0) + 1):
            self.result[rows][i] = a[i - 1]
        d = np.zeros(3)
        rs = np.zeros(3)
        alphaZ = 75
        # Параметры НИПа
        phi = math.radians(18)
        lamda = math.radians(173)
        Sg = ISZ.omega * t  # + Sg0 которое равно нулю
        NIPLongitude = Sg + lamda  # долгота НИПа
        ISZLongitude = self.Geographical(a)[0]
        rs[0] = math.cos(phi) * math.cos(Sg + lamda)
        rs[1] = math.cos(phi) * math.sin(Sg + lamda)
        rs[2] = math.sin(phi)
        self.rn = np.row_stack((self.rn, rs))
        for i in range(3):
            d[i] = a[i] - rs[i] * (ISZ.re + self.h)
        dRs = np.dot(d, rs)  # скалярное произведение d на rs
        moduleD = norm(d, ord=2)  # модуль вектора d
        moduleRs = norm(rs, ord=2)  # модуль вектора rs
        alpha = dRs / (moduleD * moduleRs)
        alpha = math.acos(alpha)
        alpha = math.degrees(alpha)
        self.d = np.row_stack((self.d, d))
        mod_d = math.sqrt(pow(d[0], 2) + pow(d[1], 2) + pow(d[2], 2))
        mod_rs = math.sqrt(pow(rs[0], 2) + pow(rs[1], 2) + pow(rs[2], 2))
        tempElevation = 90 - math.acos((rs[0] * d[0] + rs[1] * d[1] + rs[2] * d[2]) / (mod_rs * mod_d))
        arr = np.array([[-math.sin(phi) * math.cos(Sg), -math.sin(phi) * math.sin(Sg), math.cos(phi)],
                        [math.cos(phi) * math.cos(Sg), math.cos(phi) * math.sin(Sg), math.sin(phi)],
                        [-math.sin(phi), math.cos(phi), 0]])
        dt = np.dot(arr, d)
        self.dt = np.row_stack((self.dt, dt))
        # north = np.array([1, 0, 0])
        # tempAzimut = np.dot(dt, north) / 1000
        mod_dt = math.sqrt(pow(dt[0], 2) + pow(dt[1], 2) + pow(dt[2], 2))
        if d[2] < 0:
            tempAzimut = -math.acos(dt[0] / mod_dt)
        else:
            tempAzimut = math.acos(dt[0] / mod_dt)
        tempAzimut = math.degrees(tempAzimut)
        tg = ISZLongitude - NIPLongitude
        delta = self.Geographical(a)[1]
        self.OpornResult = np.row_stack((self.OpornResult, self.result[t]))
        # tempElevation = math.asin(math.sin(delta) * math.sin(phi) + math.cos(delta)
        # * math.cos(phi) * math.cos(tg))
        # tempAzimut = math.asin((math.cos(delta) * math.sin(tg)) / math.cos(self.e))
        self.ElevationAzimut = np.row_stack((self.ElevationAzimut, tempElevation))
        self.ElevationAzimut = np.row_stack((self.ElevationAzimut, tempAzimut))
        disp = 3.3  # по-моему слишком сильное отклонение происходит
        deltaD = random.normalvariate(0, disp)
        iscEl = tempElevation + deltaD
        self.iscElAz = np.row_stack((self.iscElAz, iscEl))
        iscAz = tempAzimut + deltaD
        self.iscElAz = np.row_stack((self.iscElAz, iscAz))
        self.Elevation = np.row_stack((self.Elevation, tempElevation))
        self.Azimut = np.row_stack((self.Azimut, tempAzimut))
        self.angles = np.row_stack((self.angles, alpha))
        self.count += 1
    # The core of Dorman-Prins method
    # Here we are solving the system of equations
    def get_right(self, a, t):
        module_x = math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])  # length of a vector
        result = np.empty(self.getsize(self.x0))
        # T = 0.1428571229
        # k = 1.511857892
        result[0] = a[3]  # Vx = dx/dt
        result[1] = a[4]  # Vy = dy/dt
        result[2] = a[5]  # Vz = dz/dt
        result[3] = -ISZ.mu * a[0] / math.pow(module_x, 3) + self.aerodynamic_force(a)[0]
        result[4] = -ISZ.mu * a[1] / math.pow(module_x, 3) + self.aerodynamic_force(a)[1]
        result[5] = -ISZ.mu * a[2] / math.pow(module_x, 3) + self.aerodynamic_force(a)[2]
        # result[6] = 0  # 1 / T * (k * self.wn.getvalue(t) - a[6])
        return result