from cmath import sin
from turtle import right

import sys
import math 
import numpy as np 
from numpy import sin, cos, tan, cosh, tanh, sinh
from scipy.integrate import odeint

class NonlinearEquation:
    x0 = 0
    xn = 3.14
    t0 = 0
    t1 = 1
    h = 0.01
    tau = 0.00005
    n = 100
    m = 20000
    v0formula = "0"
    vnformula = "0"
    u0formula = "np.sin(x)" 
    kformula = "1"
    fformula = "np.sin(x*t)"

    def __init__(self, _x0, _xn, _h, _tau, isKDiff, v0f, vnf, u0f, kf, ff):
        self.x0 = _x0
        self.xn = _xn
        self.h = _h
        self.tau = _tau
        self.n = int((_xn - _x0) / _h)
        self.m = int(1 / _tau)
        self.v0formula = v0f
        self.vnformula = vnf
        self.u0formula = u0f
        self.kformula = kf
        self.fformula = ff
        self.kIsDiff = isKDiff 
        #self.kvalues = np.zeros((self.m, self.n)) 

    def v0(self, x, t):
        return eval(self.v0formula, {"np": np, "math": math, "x": x, "t": t})

    def vn(self, x, t):
        return eval(self.vnformula, {"np": np, "math": math, "x": x, "t": t})

    def u0(self, x, t):
        #return np.sin(x)
        return eval(self.u0formula, {"np": np, "math": math, "x": x, "t": t})


    def u(self, xi, tj):
        if (tj == 0):
            return self.u0(xi, tj)
        if (xi == self.x0):
            return self.v0(xi, tj)
        if (xi == self.xn):
            return self.vn(xi, tj)
        return 1

    def f(self, x, t):
        return eval(self.fformula, {"np": np, "math": math, "x": x, "t": t})
    
    def k(self, u):
        return eval(self.kformula, {"np": np, "math": math, "u": u})

    def initGrid(self):

        N = self.n
        M = self.m
        tau = self.tau
        X0 = self.x0
        XN = self.xn

        gridX = np.zeros(N + 1)
        gridT = np.zeros(M + 1)
        gridY = np.zeros((M + 1, N + 1))

        for i in range(N + 1):
            gridX[i] = ((XN - X0) / N) * i

        for j in range(M + 1):
            gridT[j] = tau * j
            for i in range(N + 1):
                gridY[j][i] = self.f(gridX[i], gridT[j])

        return (gridX, gridT, gridY)


    def sweepMethod(self, leftMatrix, rightArr, N):
        a = np.zeros(N)
        b = np.zeros(N)
        result = np.zeros(N)
        y = leftMatrix[0][0]
        a[0] = -leftMatrix[0][1] / y
        b[0] = rightArr[0] / y

        for i in range(1, N - 1): 
            y = leftMatrix[i][i] + leftMatrix[i][i - 1] * a[i - 1]
            a[i] = -leftMatrix[i][i + 1] / y
            b[i] = (rightArr[i] - leftMatrix[i][i - 1] * b[i - 1]) / y

        result[N - 1] = (rightArr[N - 1] - leftMatrix[N - 1][N - 2] * b[N - 2]) / (leftMatrix[N - 1][N - 1] + leftMatrix[N - 1][N - 2] * a[N - 2])
        
        for i in range(N - 2, -1, -1):
            result[i] = a[i] * result[i + 1] + b[i]

        return result

    def calcMatrix(self, prevLayer):
        N = self.n
        M = self.m
        tau = self.tau
        h = self.h
        X0 = self.x0
        XN = self.xn
        alpha = np.zeros(N)
        beta = np.zeros(N)
        gamma = np.zeros(N)
        for i in range(1, N):
            kright = (self.k(prevLayer[i]) + self.k(prevLayer[i + 1])) / 2
            kleft = (self.k(prevLayer[i]) + self.k(prevLayer[i - 1])) / 2
            alpha[i] = -1 * tau * kleft / (h * h)
            beta[i] = (tau * (kleft + kright) + h * h) / (h * h)
            gamma[i] = -1 * tau * kright / (h * h)
        result = np.zeros((N, N))
        for i in range(N - 1):
                result[i][i] = beta[i + 1]
                result[i + 1][i] = gamma[i + 1]
                result[i][i + 1] = alpha[i + 1]
                result[N - 1][N - 1] = beta[N - 1]

        #print(result)
        return result

    def calcRightVector(self, y, gridF, j):
        N = self.n
        tau = self.tau
        result = np.zeros(N)
        for i in range(N):
            result[i] = y[i] + tau * gridF[j][i]
        return result

    def implicitSolution(self, gridX, gridT, gridF):
        N = self.n
        M = self.m
        tau = self.tau
        ans = np.zeros((M + 1, N + 1))
        rightArr = np.zeros(N)
        for i in range(N):
            ans[0][i] = self.u(gridX[i], 0)
            rightArr[i] = ans[0][i]
        for i in range(M):
            ans[i][0] = self.u(0, gridT[i])
            ans[i][N] = self.u(gridX[N], gridT[i])

        coeffMatrix = self.calcMatrix(ans[0])
        print(coeffMatrix)

        for j in range(1, M):
            coeffMatrix = self.calcMatrix(ans[j - 1])
            newLayer = self.sweepMethod(coeffMatrix, rightArr, N)
            print(newLayer)
            for i in range(len(newLayer)):
                ans[j][i] = newLayer[i]
                rightArr[i] = newLayer[i]

        return ans












