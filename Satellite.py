import numpy as np
import math


class Satellite:
    def __init__(self) -> None:
        super().__init__()
        # self.Omega = np.array([math.radians(0), math.radians(360)]) долгота восходящего узла
        # self.i = np.array([math.radians(0), math.radians(90)]) наклонение орбиты
        # self.omega = np.array([math.radians(-90), math.radians(90)]) широта перицентра
        # self.a = np.array([6500, 50000]) большая полуось орбиты
        # self.e = np.array([0, 0.5]) эксцентриситет
        # self.v = np.array([math.radians(0), math.radians(360)]) истинная аномалия
        self.i = math.radians(42)
        re = 6371  # Earth radius, km
        h_pi = 21000
        h_alpha = 970
        r_pi = h_pi + re
        r_alpha = h_alpha + re
        self.a = (r_pi + r_alpha) / 2  # km
        self.e = (r_alpha - r_pi) / (r_alpha + r_pi)  # km

        self.Omega = 0
        self.omega = 2
        self.v = 0

        start = self.isk_count()
        self.x0 = np.zeros(6)
        for i in range(6):
            self.x0[i] = start[i]

        # вычисление параметров орбиты в инерциальной прямоугольной геоцентрической системе координат
    def isk_count(self):
        mu = 398600.436  # km^3

        omega = self.Omega
        u = self.omega + self.v
        p = self.a * (1 - self.e * self.e)
        r = p / (1 + self.e * math.cos(self.v))

        vr = math.sqrt(math.fabs(mu / p)) * self.e * math.sin(self.v)
        vn = math.sqrt(math.fabs(mu / p)) * (1 + self.e * math.cos(self.v))

        x = r * (math.cos(u) * math.cos(omega) - math.sin(u) * math.sin(omega) * math.cos(self.i))
        y = r * (math.cos(u) * math.sin(omega) + math.sin(u) * math.cos(omega) * math.cos(self.i))
        z = r * math.cos(u) * math.sin(self.i)

        vx = x / r * vr + (-math.sin(u) * math.cos(omega) - math.cos(u) * math.sin(omega) * math.cos(self.i)) * vn
        vy = y / r * vr + (-math.sin(u) * math.sin(omega) + math.cos(u) * math.cos(omega) * math.cos(self.i))
        vz = z / r * vr + math.cos(u) * math.sin(self.i) * vn
        result = np.array([x, y, z, vx, vy, vz])
        return result
