import csv
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import filedialog
import numpy as np
from sympy import symbols
from sympy import sqrt, Matrix
import sympy as sy
from numpy import linalg
from FitData import FitData


class Phase3(FitData):
    def __init__(self):
        super().__init__()

    def initialise_coefHQ(self):  # On prend Nh=3
        self.getCSV()
        x_pts, y_pts, N = self.x, self.y, 3
        """ Initialise les parametres de l'HQ 

        Un jeu de paramètres est déterminé à partir d'un polygone entourant le nuage de points.
        NB : ne fonctionne correctement que pour N >= 3

        Parameters
        ----------
        x_pts, y_pts: tables de float, coordonnees du nuage de points a fitter
        N: entier, nombre de termes de l'HQ à fitter

        Returns
        -------
        param_ini : table de float de dimension (N,4), paramètres pour initialiser le fit de l'HQ

        """
        # 1) Calcul du centre de gravite et des distances max et min
        xG, yG = np.mean(x_pts), np.mean(y_pts)
        dMax = 2 * np.sqrt(np.max((x_pts - xG) ** 2 + (y_pts - yG) ** 2))
        pad = 0.1 * dMax  # Marge pour s'assurer que les droites incluent tous les points

        i, piN = np.arange(N), np.pi / N
        R = (dMax / 2 + pad) / np.cos(piN)  # Rayon du cercle circonscrit au polygone
        D = 1.5 * (dMax + pad)  # Distance entre les droites d'une paire
        # Definition des N sommets d'un polygone
        ptsPoly = np.zeros((2, N + 1))
        ptsPoly[0, :N] = xG + R * np.cos((2 * i + 1) * piN)
        ptsPoly[1, :N] = yG + R * np.sin((2 * i + 1) * piN)
        ptsPoly[0, -1], ptsPoly[1, -1] = ptsPoly[0, 0], ptsPoly[1, 0]
        # Definition de N points sur les lignes paralleles
        ptsPar = np.zeros((2, N))
        ptsPar[0, :] = xG - (D - R) * np.cos(2 * (i + 1) * piN)
        ptsPar[1, :] = yG - (D - R) * np.sin(2 * (i + 1) * piN)

        # Calcul des coefficients
        Det = ptsPoly[0, i] * ptsPoly[1, i + 1] - ptsPoly[0, i] * ptsPar[1, i] + \
              ptsPar[0, i] * ptsPoly[1, i] - ptsPar[0, i] * ptsPoly[1, i + 1] + \
              ptsPoly[0, i + 1] * ptsPar[1, i] - ptsPoly[0, i + 1] * ptsPoly[1, i]

        param_ini = np.zeros((N, 4))
        param_ini[i, 0] = 2 * (ptsPoly[1, i + 1] - ptsPoly[1, i]) / Det
        param_ini[i, 1] = 2 * (ptsPoly[0, i] - ptsPoly[0, i + 1]) / Det
        param_ini[i, 2] = (ptsPoly[0, i + 1] * ptsPar[1, i] - ptsPoly[0, i] * ptsPoly[1, i + 1] + \
                           ptsPar[0, i] * ptsPoly[1, i] - ptsPoly[0, i] * ptsPar[1, i] + \
                           ptsPoly[0, i + 1] * ptsPoly[1, i] - ptsPar[0, i] * ptsPoly[1, i + 1]) / Det
        param_ini[i, 3] = 4

        return param_ini

    def Fio(self, A, B, C, i):
        s = 0
        Nh = len(A)
        for k in range(Nh):
            s += (A[k] * self.x[i] + B[k] * self.y[i] + C[k]) ** 4
        return s ** 0.25

    def Eof1(self, A, B, C):
        N = len(self.x)
        Nh = len(A)
        s = 0
        u = 0
        for i in range(N):
            s = self.Fio(A, B, C, i)
            u += (1 - s) ** 2
        return 0.5 * u

    def d(self,A, B, C, i):
        n = 1 - self.Fio(A, B, C, i)
        s = self.Fio(A, B, C, i) ** (-3)
        Nh = len(A)
        s1, s2 = 0, 0
        for k in range(Nh):
            s1 += A[k] * (A[k] * self.x[i] + B[k] * self.y[i] + C[k]) ** 3
            s2 += B[k] * (A[k] * self.x[i] + B[k] * self.y[i] + C[k]) ** 3
        d1 = s * s1
        d2 = s * s2
        d = [d1, d2]
        d = np.linalg.norm(d)
        return n / d

    def Eof2(self):
        N = len(self.x)
        s = 0
        for i in range(N):
            s += self.d(i) ** 2
        return 0.5 * s

    # Définie Smin et Smax
    def couronneNuage(self):
        xc, yc = np.mean(self.x), np.mean(self.y)  # coordonnées du centre du nuage
        smin, smax = np.inf, 0  # rayons minimal et maximal

        for u, v in zip(self.x, self.y):
            dist = np.sqrt((u - xc) ** 2 + (v - yc) ** 2)
            smin = min(smin, dist)
            smax = max(smax, dist)

        return smin, smax

    def Pk(self, A, B, C):
        k1, k2, s = 10, 10, 0
        Nh = len(A)
        b = []
        Smin, Smax = self.couronneNuage()
        for k in range(Nh):
            u1 = (2 / (k1 * Smax)) ** 2
            u2 = (2 / (k2 * Smin)) ** 2
            m1 = u1 - (A[k] ** 2 + B[k] ** 2)
            m2 = (A[k] ** 2 + B[k] ** 2) - u2
            s += max(0, m1) ** 2 + max(0, m2) ** 2
        return s

        # définition de circle

    def circle(self, smin, smax):
        xc, yc = np.mean(self.x), np.mean(self.y)  # coordonnées du centre du nuage
        theta = np.linspace(0, 2 * np.pi, 100)
        x1 = smin * np.cos(theta)
        y1 = smin * np.sin(theta)
        x2 = smax * np.cos(theta)
        y2 = smax * np.sin(theta)

        x1 += xc
        y1 += yc
        # plt.figure()
        plt.plot(x1, y1, label='Smin')
        x2 += xc
        y2 += yc
        plt.plot(x2, y2, 'r', label='Smax')
        plt.legend()

    def Eof(self, A, B, C):
        v = 10 ** 8
        s = self.Eof2() + v * self.Pk(self, A, B, C)
        return s

    def J(self, A, B, C,i):
        i = 0
        s = self.Fio(A, B, C, i) ** (-3)
        Nh, N = len(A), len(self.x)
        s1, s2 = 0, 0
        w = [0] * N
        for i in range(N):
            for k in range(Nh):
                s1 += A[k] * (A[k] * self.x[i] + B[k] * self.y[i] + C[k]) ** 3
                s2 += B[k] * (A[k] * self.x[i] + B[k] * self.y[i] + C[k]) ** 3
            d1 = s * s1
            d2 = s * s2
            d = [d1, d2]
            d = np.linalg.norm(d) ** 2
            # w.append(1/d)
            w[i] = float(1 / d)
        p = 0
        for i in range(N):
            p += w[i] * (1 - self.Fio(A, B, C, i)) ** 2
        p = p * 0.5
        for i in range(N):
            p += 10 ** 8 * self.Pk(A, B, C)
        return p

    def derivee(self, f):
        a1, b1, c1, a2, b2, c2, a3, b3, c3, x, y = symbols('a1 b1 c1 a2 b2 c2 a3 b3 c3 x y')
        f = sy.sympify(f)
        print("J(x)= {}".format(f))

        # 1 ere dérivée:
        d1fa1 = sy.diff(f, a1)
        d1fb1 = sy.diff(f, b1)
        d1fc1 = sy.diff(f, c1)
        d1fa2 = sy.diff(f, a2)
        d1fb2 = sy.diff(f, b2)
        d1fc2 = sy.diff(f, c2)
        d1fa3 = sy.diff(f, a3)
        d1fb3 = sy.diff(f, b3)
        d1fc3 = sy.diff(f, c3)

        d1L = [d1fa1, d1fb1, d1fc1, d1fa2, d1fb2, d1fc2, d1fa3, d1fb3, d1fc3]
        print('J\'(x,y)=: {}'.format(d1L))

        # convertir en fonction
        d1fa1 = sy.lambdify(a1, d1fa1)
        d1fb1 = sy.lambdify(b1, d1fb1)
        d1fc1 = sy.lambdify(c1, d1fc1)
        d1fa2 = sy.lambdify(a2, d1fa2)
        d1fb2 = sy.lambdify(b2, d1fb2)
        d1fc2 = sy.lambdify(c2, d1fc2)
        d1fa3 = sy.lambdify(a3, d1fa3)
        d1fb3 = sy.lambdify(b3, d1fb3)
        d1fc3 = sy.lambdify(c3, d1fc3)
        d1f = [d1fa1, d1fb1, d1fc1, d1fa2, d1fb2, d1fc2, d1fa3, d1fb3, d1fc3]


        print("**fin de calcul de derivées 1 ordre et 2eme ordre ***")
        return d1f, d1L

    # Smin,Smax=self.couronneNuage()

    # u1=(2/(10*Smax))**2
    # u2=(2/(10*Smin))**2
    def gradJ(self, A, B, C, i):
        s = self.Fio(A, B, C, i) ** (-3)
        Nh, N = len(A), len(self.x)
        s1, s2 = 0, 0
        w = [0] * N
        for i in range(N):
            for k in range(Nh):
                s1 += A[k] * (A[k] * self.x[i] + B[k] * self.y[i] + C[k]) ** 3
                s2 += B[k] * (A[k] * self.x[i] + B[k] * self.y[i] + C[k]) ** 3
            d1 = s * s1
            d2 = s * s2
            d = [d1, d2]
            d = np.linalg.norm(d) ** 2
            # w.append(1/d)
            w[i] = float(1 / d)
        Smin, Smax = self.couronneNuage()
        u1 = (2 / (10 * Smax)) ** 2
        u2 = (2 / (10 * Smin)) ** 2
        pa1 = A[0] * (-max(0, u1 - (A[0] ** 2 + B[0] ** 2)) + max(0, -u2 + (A[0] ** 2 + B[0] ** 2)))
        pa2 = A[1] * (-max(0, u1 - (A[1] ** 2 + B[1] ** 2)) + max(0, -u2 + (A[1] ** 2 + B[1] ** 2)))
        pa3 = A[2] * (-max(0, u1 - (A[2] ** 2 + B[2] ** 2)) + max(0, -u2 + (A[2] ** 2 + B[2] ** 2)))
        pb1 = B[0] * (-max(0, u1 - (A[0] ** 2 + B[0] ** 2)) + max(0, -u2 + (A[0] ** 2 + B[0] ** 2)))
        pb2 = B[1] * (-max(0, u1 - (A[1] ** 2 + B[1] ** 2)) + max(0, -u2 + (A[1] ** 2 + B[1] ** 2)))
        pb3 = B[2] * (-max(0, u1 - (A[2] ** 2 + B[2] ** 2)) + max(0, -u2 + (A[2] ** 2 + B[2] ** 2)))
        pc1 = 0
        pc2 = 0
        pc3 = 0
        p = [[] for i in range(9)]

        p[0] = pa1
        p[1] = pb1
        p[2] = pc1
        p[3] = pa2
        p[4] = pb2
        p[5] = pc2
        p[6] = pa3
        p[7] = pb3
        p[8] = pc3
        p = np.asarray(p).T

        s = 0
        for i in range(N):
            d1fa1 = -8 * self.x[i] * (A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 3 * (
                        -(A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 4 + (
                            A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 4 + (
                                    A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 4 + 1)
            d1fb1 = -8 * self.y[i] * (A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 3 * (
                        -(A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 4 + (
                            A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 4 + (
                                    A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 4 + 1)
            d1fc1 = -8 * (A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 3 * (
                        -(A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 4 + (
                            A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 4 + (
                                    A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 4 + 1)
            d1fa2 = 8 * self.x[i] * (A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 3 * (
                        -(A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 4 + (
                            A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 4 + (
                                    A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 4 + 1)
            d1fb2 = 8 * self.y[i] * (A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 3 * (
                    -(A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 4 + (
                            A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 4 + (
                            A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 4 + 1)
            d1fc2 = 8 * (A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 3 * (
                        -(A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 4 + (
                            A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 4 + (
                                    A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 4 + 1)
            d1fa3 = 8 * self.x[i] * (A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 3 * (
                    -(A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 4 + (
                    A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 4 + (
                                    A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 4 + 1)
            d1fb3 = 8 * self.y[i] * (A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 3 * (
                    -(A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 4 + (
                    A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 4 + (
                            A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 4 + 1)
            d1fc3 = 8 * (A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 3 * (
                    -(A[0] * self.x[i] + B[0] * self.y[i] + C[0]) ** 4 + (
                    A[1] * self.x[i] + B[1] * self.y[i] + C[1]) ** 4 + (
                            A[2] * self.x[i] + B[2] * self.y[i] + C[2]) ** 4 + 1)
            J = [[] for i in range(9)]

            J[0] = d1fa1
            J[1] = d1fb1
            J[2] = d1fc1
            J[3] = d1fa2
            J[4] = d1fb2
            J[5] = d1fc2
            J[6] = d1fa3
            J[7] = d1fb3
            J[8] = d1fc3
            J = np.asarray(J).T
            # J=np.array([J[0],J[1],J[2],J[3],J[4],J[5],J[6],J[7],J[8]])
            # J=float(J)
            s += w[i] * J
        s = 0.5 * s + (10 ** 8) * p
        return s

    def Hessienne(self, A, B, C):
        ss, i = 0, 0
        s = self.Fio(A, B, C, i) ** (-3)
        Nh, N = len(A), len(self.x)
        s1, s2 = 0, 0
        w = [0] * N
        for i in range(N):
            for k in range(Nh):
                s1 += A[k] * (A[k] * self.x[i] + B[k] * self.y[i] + C[k]) ** 3
                s2 += B[k] * (A[k] * self.x[i] + B[k] * self.y[i] + C[k]) ** 3
            d1 = s * s1
            d2 = s * s2
            d = [d1, d2]
            d = np.linalg.norm(d) ** 2
            # w.append(1/d)
            w[i] = float(1 / d)

        for i in range(N):
            x = self.x[i]
            y = self.y[i]
            a1, a2, a3 = A[0], A[1], A[2]
            b1, b2, b3 = B[0], B[1], B[2]
            c1, c2, c3 = C[0], C[1], C[2]
            # H1=np.asarray([[[32*x**2*(a1*x + b1*y + c1)**6 - 24*x**2*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*x*y*(a1*x + b1*y + c1)**6 - 24*x*y*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*x*(a1*x + b1*y + c1)**6 - 24*x*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), -32*x**2*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*x*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*x*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*x**2*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*x*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*x*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3]], [[32*x*y*(a1*x + b1*y + c1)**6 - 24*x*y*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*y**2*(a1*x + b1*y + c1)**6 - 24*y**2*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*y*(a1*x + b1*y + c1)**6 - 24*y*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), -32*x*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*y**2*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*x*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*y**2*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3]], [[32*x*(a1*x + b1*y + c1)**6 - 24*x*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*y*(a1*x + b1*y + c1)**6 - 24*y*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*(a1*x + b1*y + c1)**6 - 24*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), -32*x*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*x*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3]], [[-32*x**2*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*x*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*x*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, 32*x**2*(a2*x + b2*y + c2)**6 + 24*x**2*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*x*y*(a2*x + b2*y + c2)**6 + 24*x*y*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*x*(a2*x + b2*y + c2)**6 + 24*x*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*x**2*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*x*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*x*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3]], [[-32*x*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*y**2*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, 32*x*y*(a2*x + b2*y + c2)**6 + 24*x*y*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*y**2*(a2*x + b2*y + c2)**6 + 24*y**2*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*y*(a2*x + b2*y + c2)**6 + 24*y*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*x*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*y**2*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3]], [[-32*x*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, -32*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3, 32*x*(a2*x + b2*y + c2)**6 + 24*x*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*y*(a2*x + b2*y + c2)**6 + 24*y*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*(a2*x + b2*y + c2)**6 + 24*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*x*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3]], [[-32*x**2*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*x*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*x*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, 32*x**2*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*x*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*x*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*x**2*(a3*x + b3*y + c3)**6 + 24*x**2*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*x*y*(a3*x + b3*y + c3)**6 + 24*x*y*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*x*(a3*x + b3*y + c3)**6 + 24*x*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1)]], [[-32*x*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*y**2*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, 32*x*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*y**2*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*x*y*(a3*x + b3*y + c3)**6 + 24*x*y*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*y**2*(a3*x + b3*y + c3)**6 + 24*y**2*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*y*(a3*x + b3*y + c3)**6 + 24*y*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1)]], [[-32*x*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, -32*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3, 32*x*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3, 32*x*(a3*x + b3*y + c3)**6 + 24*x*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*y*(a3*x + b3*y + c3)**6 + 24*y*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1), 32*(a3*x + b3*y + c3)**6 + 24*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1)]]])
            # H1=[[[32*self.x[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**6 - 24*self.x[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**6 - 24*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**6 - 24*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), -32*self.x[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.x[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3]], [[32*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**6 - 24*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.y[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**6 - 24*self.y[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**6 - 24*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), -32*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.y[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*self.y[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3]], [[32*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**6 - 24*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**6 - 24*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**6 - 24*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), -32*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3]], [[-32*self.x[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, 32*self.x[i]**2*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**6 + 24*self.x[i]**2*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.x[i]*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**6 + 24*self.x[i]*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.x[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**6 + 24*self.x[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.x[i]**2*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.x[i]*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.x[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3]], [[-32*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.y[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, 32*self.x[i]*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**6 + 24*self.x[i]*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.y[i]**2*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**6 + 24*self.y[i]**2*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**6 + 24*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.x[i]*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.y[i]**2*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3]], [[-32*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, -32*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3, 32*self.x[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**6 + 24*self.x[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**6 + 24*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**6 + 24*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.x[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3]], [[-32*self.x[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.x[i]**2*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.x[i]*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.x[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.x[i]**2*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**6 + 24*self.x[i]**2*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.x[i]*self.y[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**6 + 24*self.x[i]*self.y[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.x[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**6 + 24*self.x[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1)]], [[-32*self.x[i]*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*self.y[i]**2*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.x[i]*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.y[i]**2*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.x[i]*self.y[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**6 + 24*self.x[i]*self.y[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.y[i]**2*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**6 + 24*self.y[i]**2*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.y[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**6 + 24*self.y[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1)]], [[-32*self.x[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*self.y[i]*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, -32*(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.x[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.y[i]*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*(A[1]*self.x[i] + B[1]*self.y[i] + C[1])**3*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**3, 32*self.x[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**6 + 24*self.x[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*self.y[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**6 + 24*self.y[i]*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1), 32*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**6 + 24*(A[2]*self.x[i] + B[2]*self.y[i] + C[2])**2*(-(A[0]*self.x[i] + B[0]*self.y[i] + C[0])**4 + (A[1]*self.x[i] + B[1]*self.y[i] + C[1])**4 + (A[2]*self.x[i] + B[2]*self.y[i] + C[2])**4 + 1)]]]
            # H1=[[[32*x**2*(a1*x + b1*y + c1)**6, - 24*x**2*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1)  ,32*x*y*(a1*x + b1*y + c1)**6 - 24*x*y*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1)  ,32*x*(a1*x + b1*y + c1)**6 - 24*x*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1)  ,-32*x**2*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*x*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3,  -32*x*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3,  -32*x**2*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , -32*x*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3  ,-32*x*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3] [32*x*y*(a1*x + b1*y + c1)**6 - 24*x*y*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*y**2*(a1*x + b1*y + c1)**6 - 24*y**2*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*y*(a1*x + b1y + c1)**6 - 24*y*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , -32*x*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*y**2*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*x*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 ,-32*y**2*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , -32*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3] [32*x*(a1*x + b1*y + c1)**6 - 24*x*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*y*(a1*x + b1*y + c1)**6 - 24*y*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*(a1*x + b1*y + c1)**6 - 24*(a1*x + b1*y + c1)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , -32*x*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*x*(a1*x + b1*y + c1)**3*(a3*x + b3*y  3)**3 , -32*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , -32*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3] [-32*x**2*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*x*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*x*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , 32*x**2*(a2*x + b2*y + c2)**6 + 24*x**2*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*x*y*(a2*x + b2*y + c2)**6 + 24*x*y*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*x*(a2*x + b2*y + c2)**6 + 24*x*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*x**2*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*x*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*x*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3] [-32*x*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*y**2*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , 32*x*y*(a2*x + b2*y + c2)**6 + 2*x*y*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*y**2*(a2*x + b2*y + c2)**6 + 24*y**2*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*y*(a2*x + b2*y + c2)**6 + 24*y*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*x*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*y**2*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3] [-32*x*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*y*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , -32*(a1*x + b1*y + c1)**3*(a2*x + b2*y + c2)**3 , 32*x*(a2*x + b2*y + c2)**6 + 24*x*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1),  32*y*(a2*x + b2*y + c2)**6 + 24*y*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*(a2*x + b2*y + c2)**6 + 24*(a2*x + b2*y + c2)**2*(-(a1*x + b1*y + c)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*x*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3] [-32*x**2*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , -32*x*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , -32*x*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , 32*x**2*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3  ,32*x*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*x*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3  ,32*x**2*(a3*x + b3*y + c3)**6 + 24*x**2*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*x*y*(a3*x + b3*y + c3)**6 + 24*x*y*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*x*(a3*x + b3*y + c3)**6 + 24*x*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1)] [-32*x*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , -32*y**2*(a1*x + b1*y + c1)**3*(a3*x  3*y + c3)**3 , -32*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , 32*x*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*y**2*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*x*y*(a3*x + b3*y + c3)**6 + 24*x*y*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*y**2*(a3*x + b3*y + c3)**6 + 24*y**2*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*y*(a3*x + b3*y + c3)**6 + 24*y*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1)] [-32*x*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , -32*y*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , -32*(a1*x + b1*y + c1)**3*(a3*x + b3*y + c3)**3 , 32*x*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3,  32*y*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*(a2*x + b2*y + c2)**3*(a3*x + b3*y + c3)**3 , 32*x*(a3*x + b3*y + c3)**6 + 24*x*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*y*(a3*x + b3*y + c3)**6 + 24*y*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1) , 32*(a3*x + b3*y + c3)**6 + 24*(a3*x + b3*y + c3)**2*(-(a1*x + b1*y + c1)**4 + (a2*x + b2*y + c2)**4 + (a3*x + b3*y + c3)**4 + 1)]]]
            H1 = [[32 * x ** 2 * (a1 * x + b1 * y + c1) ** 6 - 24 * x ** 2 * (a1 * x + b1 * y + c1) ** 2 * (
                        -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * x * y * (a1 * x + b1 * y + c1) ** 6 - 24 * x * y * (a1 * x + b1 * y + c1) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * x * (a1 * x + b1 * y + c1) ** 6 - 24 * x * (a1 * x + b1 * y + c1) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   -32 * x ** 2 * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   -32 * x * y * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   -32 * x * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   -32 * x ** 2 * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   -32 * x * y * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   -32 * x * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3], [
                      32 * x * y * (a1 * x + b1 * y + c1) ** 6 - 24 * x * y * (a1 * x + b1 * y + c1) ** 2 * (
                                  -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                      a3 * x + b3 * y + c3) ** 4 + 1),
                      32 * y ** 2 * (a1 * x + b1 * y + c1) ** 6 - 24 * y ** 2 * (a1 * x + b1 * y + c1) ** 2 * (
                                  -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                      a3 * x + b3 * y + c3) ** 4 + 1),
                      32 * y * (a1 * x + b1 * y + c1) ** 6 - 24 * y * (a1 * x + b1 * y + c1) ** 2 * (
                                  -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                      a3 * x + b3 * y + c3) ** 4 + 1),
                      -32 * x * y * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                      -32 * y ** 2 * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                      -32 * y * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                      -32 * x * y * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                      -32 * y ** 2 * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                      -32 * y * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3], [
                      32 * x * (a1 * x + b1 * y + c1) ** 6 - 24 * x * (a1 * x + b1 * y + c1) ** 2 * (
                                  -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                      a3 * x + b3 * y + c3) ** 4 + 1),
                      32 * y * (a1 * x + b1 * y + c1) ** 6 - 24 * y * (a1 * x + b1 * y + c1) ** 2 * (
                                  -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                      a3 * x + b3 * y + c3) ** 4 + 1),
                      32 * (a1 * x + b1 * y + c1) ** 6 - 24 * (a1 * x + b1 * y + c1) ** 2 * (
                                  -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                      a3 * x + b3 * y + c3) ** 4 + 1),
                      -32 * x * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                      -32 * y * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                      -32 * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                      -32 * x * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                      -32 * y * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                      -32 * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3],
                  [-32 * x ** 2 * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   -32 * x * y * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   -32 * x * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   32 * x ** 2 * (a2 * x + b2 * y + c2) ** 6 + 24 * x ** 2 * (a2 * x + b2 * y + c2) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * x * y * (a2 * x + b2 * y + c2) ** 6 + 24 * x * y * (a2 * x + b2 * y + c2) ** 2 * (
                           -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * x * (a2 * x + b2 * y + c2) ** 6 + 24 * x * (a2 * x + b2 * y + c2) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * x ** 2 * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * x * y * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * x * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3],
                  [-32 * x * y * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   -32 * y ** 2 * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   -32 * y * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   32 * x * y * (a2 * x + b2 * y + c2) ** 6 + 24 * x * y * (a2 * x + b2 * y + c2) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * y ** 2 * (a2 * x + b2 * y + c2) ** 6 + 24 * y ** 2 * (a2 * x + b2 * y + c2) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * y * (a2 * x + b2 * y + c2) ** 6 + 24 * y * (a2 * x + b2 * y + c2) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * x * y * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * y ** 2 * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * y * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3],
                  [-32 * x * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   -32 * y * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   -32 * (a1 * x + b1 * y + c1) ** 3 * (a2 * x + b2 * y + c2) ** 3,
                   32 * x * (a2 * x + b2 * y + c2) ** 6 + 24 * x * (a2 * x + b2 * y + c2) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * y * (a2 * x + b2 * y + c2) ** 6 + 24 * y * (a2 * x + b2 * y + c2) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * (a2 * x + b2 * y + c2) ** 6 + 24 * (a2 * x + b2 * y + c2) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * x * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * y * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3],
                  [-32 * x ** 2 * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   -32 * x * y * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   -32 * x * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * x ** 2 * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * x * y * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * x * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * x ** 2 * (a3 * x + b3 * y + c3) ** 6 + 24 * x ** 2 * (a3 * x + b3 * y + c3) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * x * y * (a3 * x + b3 * y + c3) ** 6 + 24 * x * y * (a3 * x + b3 * y + c3) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * x * (a3 * x + b3 * y + c3) ** 6 + 24 * x * (a3 * x + b3 * y + c3) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1)],
                  [-32 * x * y * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   -32 * y ** 2 * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   -32 * y * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * x * y * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * y ** 2 * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * y * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * x * y * (a3 * x + b3 * y + c3) ** 6 + 24 * x * y * (a3 * x + b3 * y + c3) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * y ** 2 * (a3 * x + b3 * y + c3) ** 6 + 24 * y ** 2 * (a3 * x + b3 * y + c3) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * y * (a3 * x + b3 * y + c3) ** 6 + 24 * y * (a3 * x + b3 * y + c3) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1)],
                  [-32 * x * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   -32 * y * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   -32 * (a1 * x + b1 * y + c1) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * x * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * y * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * (a2 * x + b2 * y + c2) ** 3 * (a3 * x + b3 * y + c3) ** 3,
                   32 * x * (a3 * x + b3 * y + c3) ** 6 + 24 * x * (a3 * x + b3 * y + c3) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * y * (a3 * x + b3 * y + c3) ** 6 + 24 * y * (a3 * x + b3 * y + c3) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1),
                   32 * (a3 * x + b3 * y + c3) ** 6 + 24 * (a3 * x + b3 * y + c3) ** 2 * (
                               -(a1 * x + b1 * y + c1) ** 4 + (a2 * x + b2 * y + c2) ** 4 + (
                                   a3 * x + b3 * y + c3) ** 4 + 1)]]
            ss += w[i] * np.asarray(H1)
        ss = 0.5 * ss
        # Hessienne de la pénalité
        Smin, Smax = self.couronneNuage()
        u1 = (2 / (10 * Smax)) ** 2
        u2 = (2 / (10 * Smin)) ** 2
        pa1a1 = 4 * (-max(0, u1 - (A[0] ** 2 + B[0] ** 2)) + max(0, -u2 + (A[0] ** 2 + B[0] ** 2))) + 8 * (
                    A[0] ** 2) * (np.sign(max(0, u1 - (A[0] ** 2 + B[0] ** 2))) + np.sign(
            max(0, -u2 + (A[0] ** 2 + B[0] ** 2))))
        pa2a2 = 4 * (-max(0, u1 - (A[1] ** 2 + B[1] ** 2)) + max(0, -u2 + (A[1] ** 2 + B[1] ** 2))) + 8 * (
                    A[1] ** 2) * (np.sign(max(0, u1 - (A[1] ** 2 + B[1] ** 2))) + np.sign(
            max(0, -u2 + (A[1] ** 2 + B[1] ** 2))))
        pa3a3 = 4 * (-max(0, u1 - (A[2] ** 2 + B[2] ** 2)) + max(0, -u2 + (A[2] ** 2 + B[2] ** 2))) + 8 * (
                    A[2] ** 2) * (np.sign(max(0, u1 - (A[2] ** 2 + B[2] ** 2))) + np.sign(
            max(0, -u2 + (A[2] ** 2 + B[2] ** 2))))

        pb1b1 = 4 * (-max(0, u1 - (A[0] ** 2 + B[0] ** 2)) + max(0, -u2 + (A[0] ** 2 + B[0] ** 2))) + 8 * (
                    B[0] ** 2) * (np.sign(max(0, u1 - (A[0] ** 2 + B[0] ** 2))) + np.sign(
            max(0, -u2 + (A[0] ** 2 + B[0] ** 2))))
        pb2b2 = 4 * (-max(0, u1 - (A[1] ** 2 + B[1] ** 2)) + max(0, -u2 + (A[1] ** 2 + B[1] ** 2))) + 8 * (
                    B[1] ** 2) * (np.sign(max(0, u1 - (A[1] ** 2 + B[1] ** 2))) + np.sign(
            max(0, -u2 + (A[1] ** 2 + B[1] ** 2))))
        pb3b3 = 4 * (-max(0, u1 - (A[2] ** 2 + B[2] ** 2)) + max(0, -u2 + (A[2] ** 2 + B[2] ** 2))) + 8 * (
                    B[2] ** 2) * (np.sign(max(0, u1 - (A[2] ** 2 + B[2] ** 2))) + np.sign(
            max(0, -u2 + (A[2] ** 2 + B[2] ** 2))))

        pa1b1 = 8 * (A[0] * B[0]) * (
                    np.sign(max(0, u1 - (A[0] ** 2 + B[0] ** 2))) + np.sign(max(0, -u2 + (A[0] ** 2 + B[0] ** 2))))
        pa2b1 = 8 * (A[1] * B[0]) * (
                    np.sign(max(0, u1 - (A[1] ** 2 + B[0] ** 2))) + np.sign(max(0, -u2 + (A[1] ** 2 + B[0] ** 2))))
        pa3b1 = 8 * (A[2] * B[0]) * (
                    np.sign(max(0, u1 - (A[2] ** 2 + B[0] ** 2))) + np.sign(max(0, -u2 + (A[2] ** 2 + B[0] ** 2))))
        pa1b2 = 8 * (A[0] * B[1]) * (
                    np.sign(max(0, u1 - (A[0] ** 2 + B[1] ** 2))) + np.sign(max(0, -u2 + (A[0] ** 2 + B[1] ** 2))))
        pa2b2 = 8 * (A[1] * B[1]) * (
                    np.sign(max(0, u1 - (A[1] ** 2 + B[1] ** 2))) + np.sign(max(0, -u2 + (A[1] ** 2 + B[1] ** 2))))
        pa3b2 = 8 * (A[2] * B[1]) * (
                    np.sign(max(0, u1 - (A[2] ** 2 + B[1] ** 2))) + np.sign(max(0, -u2 + (A[2] ** 2 + B[1] ** 2))))
        pa1b3 = 8 * (A[0] * B[2]) * (
                    np.sign(max(0, u1 - (A[0] ** 2 + B[2] ** 2))) + np.sign(max(0, -u2 + (A[0] ** 2 + B[2] ** 2))))
        pa2b3 = 8 * (A[1] * B[2]) * (
                    np.sign(max(0, u1 - (A[1] ** 2 + B[2] ** 2))) + np.sign(max(0, -u2 + (A[1] ** 2 + B[2] ** 2))))
        pa3b3 = 8 * (A[2] * B[2]) * (
                    np.sign(max(0, u1 - (A[2] ** 2 + B[2] ** 2))) + np.sign(max(0, -u2 + (A[2] ** 2 + B[2] ** 2))))

        H2 = np.zeros((9, 9))
        H2[0, 0], H2[0, 1], H2[0, 2], H2[0, 3], H2[0, 4], H2[0, 5], H2[0, 6], H2[0, 7], H2[
            0, 8] = pa1a1, pa1b1, 0, 0, pa1b2, 0, 0, pa1b3, 0
        H2[1, 0], H2[1, 1], H2[1, 2], H2[1, 3], H2[1, 4], H2[1, 5], H2[1, 6], H2[1, 7], H2[
            1, 8] = pa1b1, pb1b1, 0, pa2b1, 0, 0, pa3b1, 0, 0
        H2[2, 0], H2[2, 1], H2[2, 2], H2[2, 3], H2[2, 4], H2[2, 5], H2[2, 6], H2[2, 7], H2[
            2, 8] = 0, 0, 0, 0, 0, 0, 0, 0, 0
        H2[3, 0], H2[3, 1], H2[3, 2], H2[3, 3], H2[3, 4], H2[3, 5], H2[3, 6], H2[3, 7], H2[
            3, 8] = 0, pa2b1, 0, pa2a2, pa2b2, 0, 0, pa2b3, 0
        H2[4, 0], H2[4, 1], H2[4, 2], H2[4, 3], H2[4, 4], H2[4, 5], H2[4, 6], H2[4, 7], H2[
            4, 8] = pa1b2, 0, 0, pa2b2, pb2b2, 0, pa3b2, 0, 0
        H2[5, 0], H2[5, 1], H2[5, 2], H2[5, 3], H2[5, 4], H2[5, 5], H2[5, 6], H2[5, 7], H2[
            5, 8] = 0, 0, 0, 0, 0, 0, 0, 0, 0
        H2[6, 0], H2[6, 1], H2[6, 2], H2[6, 3], H2[6, 4], H2[6, 5], H2[6, 6], H2[6, 7], H2[
            6, 8] = 0, pa3b1, 0, 0, pa3b2, 0, pa3a3, pa3b3, 0
        H2[7, 0], H2[7, 1], H2[7, 2], H2[7, 3], H2[7, 4], H2[7, 5], H2[7, 6], H2[7, 7], H2[
            7, 8] = pa1b3, 0, 0, pa2b3, 0, 0, pa3b3, pb3b3, 0
        H2[8, 0], H2[8, 1], H2[8, 2], H2[8, 3], H2[8, 4], H2[8, 5], H2[8, 6], H2[8, 7], H2[
            8, 8] = 0, 0, 0, 0, 0, 0, 0, 0, 0

        return ss + H2

    def LM(self):
        lam = self.initialise_coefHQ()
        eps, n = 10 ** (-6), 0
        dlam, b = 1, 0.01
        nmax = 10
        N = len(self.x)
        A = lam[:, 0]
        B = lam[:, 1]
        C = lam[:, 2]
        while (dlam > eps) and (n < nmax):
            dj = self.gradJ(A, B, C, 0)
            H = self.Hessienne(A, B, C)
            # while(self.J(A1,B1,C1)>=self.J(A,B,C)):
            g = np.zeros((9, 9))
            for i in range(9):
                g[i, i] = 1
            G = H + b * g
            Gi = -np.linalg.inv(G)
            deltalam = np.matmul(Gi, dj)
            deltalam = np.reshape(deltalam, (3, 3))
            lam1 = lam[0:3, 0:3] + deltalam
            A1 = lam1[:, 0]
            B1 = lam1[:, 1]
            C1 = lam1[:, 2]
            while (self.J(A1, B1, C1) >= self.J(A, B, C)):
                b *= 10
                G = H + b * g
                Gi = -np.linalg.inv(G)
                deltalam = np.matmul(Gi, dj)
                deltalam = np.reshape(deltalam, (3, 3))
                lam1 = lam[0:3, 0:3] + deltalam
                A1 = lam1[:, 0]
                B1 = lam1[:, 1]
                C1 = lam1[:, 2]

            b *= 0.1
            # lam=lam[0:3,0:3]
            lam = lam[0:3, 0:3] + deltalam
            dlam = np.linalg.norm(deltalam)
            n += 1
        return lam

# phase 3
def phase3():
    phase3 = Phase3()
    lam = phase3.initialise_coefHQ()
    # print(lam)
    plt.figure()
    phase3.plotData()  # Data Scatter
    smin, smax = phase3.couronneNuage()
    phase3.circle(smin, smax)
    lam = phase3.LM()  # Renvoie lambda sous forme de matrice 3x3 contenant les élements Ak,Bk et Ck (Nh=3)
    print(lam)
    phase3.A = lam[:, 0]
    phase3.B = lam[:, 1]
    phase3.C = lam[:, 2]

    def f(x, y):
        s = 0
        for i in range(3):
            s = s + abs(phase3.A[i] * phase3.x[i] + phase3.B[i] * phase3.y[i] + phase3.C[i]) ** 4
        return s

    # phase3.graph("levenberg marquardt")
    xmin, xmax, nx = 25, 200, 4410
    tab_x = np.linspace(xmin, xmax, nx)
    ymin, ymax, ny = 40, 160, 4410
    tab_y = np.linspace(ymin, ymax, nx)

    x1d = np.linspace(xmin, xmax, nx)
    y1d = np.linspace(ymin, ymax, ny)
    x2d, y2d = np.meshgrid(x1d, y1d)

    # Tracé des isovaleurs de f1
    nIso = 21

    print("test")
    plt.figure()
    phase3.plotData()
    plt.contour(x2d, y2d, f(x2d, y2d), nIso)

phase3()










