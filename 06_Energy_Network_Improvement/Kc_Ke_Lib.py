"""core library about Kc & Ke factors for different types of HEX"""
import numpy as np
from scipy import interpolate


class Kc_Ke_Lib:
    """for rectangular passages of the HEX"""

    def Kc_rect(self, delta, Re):
        x = np.loadtxt('rectangular_passage_Kc.csv', dtype=float, delimiter=',', comments='Kc')
        Kc = np.zeros(shape=Re.shape)
        for i in range(Kc.shape[0]):
            if Re[i] < 2000:
                Kc[i] = -0.4252 * delta * delta + 0.0326 * delta + 1.1798
            elif Re[i] > 1e4:
                Kc[i] = -0.419 * delta * delta + 0.0208 * delta + 0.4033
            else:
                delta_set = x[:, 0]
                Re_set = [2e3, 3e3, 5e3, 1e4]
                Kc_set_T = np.array([(x[:, 2]), (x[:, 3]), (x[:, 4]), (x[:, 5])])
                f = interpolate.interp2d(delta_set, Re_set, Kc_set_T, kind='linear')
                Kc[i] = f(delta, Re[i])
        return Kc

    def Ke_rect(self, delta, Re):
        x = np.loadtxt('rectangular_passage_Ke.csv', dtype=float, delimiter=',', comments='Ke')
        Ke = np.zeros(shape=Re.shape)
        for i in range(Ke.shape[0]):
            if Re[i] < 2000:
                Ke[i] = 0.996 * delta * delta - 2.7687 * delta + 1.0016
            elif Re[i] > 1e4:
                Ke[i] = 0.9832 * delta * delta - 1.9823 * delta + 1
            else:
                delta_set = x[:, 0]
                Re_set = [2e3, 3e3, 5e3, 1e4]
                Ke_set_T = np.array([(x[:, 2]), (x[:, 3]), (x[:, 4]), (x[:, 5])])
                f = interpolate.interp2d(delta_set, Re_set, Ke_set_T, kind='linear')
                Ke[i] = f(delta, Re[i])
        return Ke


class Kc_Ke_Lib_double_values:
    """for rectangular passages of the HEX"""

    def Kc_rect(self, delta, Re):
        x = np.loadtxt('rectangular_passage_Kc.csv', dtype=float, delimiter=',', comments='Kc')
        Kc = np.zeros(shape=Re.shape)
        for i in range(Kc.shape[0]):
            if Re[i] < 2000:
                Kc[i] = -0.4252 * delta[i] * delta[i] + 0.0326 * delta[i] + 1.1798
            elif Re[i] > 1e4:
                Kc[i] = -0.419 * delta[i] * delta[i] + 0.0208 * delta[i] + 0.4033
            else:
                delta_set = x[:, 0]
                Re_set = [2e3, 3e3, 5e3, 1e4]
                Kc_set_T = np.array([(x[:, 2]), (x[:, 3]), (x[:, 4]), (x[:, 5])])
                f = interpolate.interp2d(delta_set, Re_set, Kc_set_T, kind='linear')
                Kc[i] = f(delta[i], Re[i])
        return Kc

    def Ke_rect(self, delta, Re):
        x = np.loadtxt('rectangular_passage_Ke.csv', dtype=float, delimiter=',', comments='Ke')
        Ke = np.zeros(shape=Re.shape)
        for i in range(Ke.shape[0]):
            if Re[i] < 2000:
                Ke[i] = 0.996 * delta[i] * delta[i] - 2.7687 * delta[i] + 1.0016
            elif Re[i] > 1e4:
                Ke[i] = 0.9832 * delta[i] * delta[i] - 1.9823 * delta[i] + 1
            else:
                delta_set = x[:, 0]
                Re_set = [2e3, 3e3, 5e3, 1e4]
                Ke_set_T = np.array([(x[:, 2]), (x[:, 3]), (x[:, 4]), (x[:, 5])])
                f = interpolate.interp2d(delta_set, Re_set, Ke_set_T, kind='linear')
                Ke[i] = f(delta[i], Re[i])
        return Ke


class Kc_Ke_Lib_single_value:
    """for rectangular passages of the HEX"""

    def Kc_rect(self, delta, Re):
        x = np.loadtxt('rectangular_passage_Kc.csv', dtype=float, delimiter=',', comments='Kc')
        if Re < 2000:
            Kc = -0.4252 * delta * delta + 0.0326 * delta + 1.1798
        elif Re > 1e4:
            Kc = -0.419 * delta * delta + 0.0208 * delta + 0.4033
        else:
            delta_set = x[:, 0]
            Re_set = [2e3, 3e3, 5e3, 1e4]
            Kc_set_T = np.array([(x[:, 2]), (x[:, 3]), (x[:, 4]), (x[:, 5])])
            f = interpolate.interp2d(delta_set, Re_set, Kc_set_T, kind='linear')
            Kc = f(delta, Re)
        return Kc

    def Ke_rect(self, delta, Re):
        x = np.loadtxt('rectangular_passage_Ke.csv', dtype=float, delimiter=',', comments='Ke')
        if Re < 2000:
            Ke = 0.996 * delta * delta - 2.7687 * delta + 1.0016
        elif Re > 1e4:
            Ke = 0.9832 * delta * delta - 1.9823 * delta + 1
        else:
            delta_set = x[:, 0]
            Re_set = [2e3, 3e3, 5e3, 1e4]
            Ke_set_T = np.array([(x[:, 2]), (x[:, 3]), (x[:, 4]), (x[:, 5])])
            f = interpolate.interp2d(delta_set, Re_set, Ke_set_T, kind='linear')
            Ke = f(delta, Re)
        return Ke
