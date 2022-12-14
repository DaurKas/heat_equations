from cmath import sin
from turtle import right

import sys
import math 
import numpy as np 
from numpy import sin, cos, tan, cosh, tanh, sinh
from scipy.integrate import odeint

class linearEquation:
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
    
    def k(self, x, t):
        return eval(self.kformula, {"np": np, "math": math, "x": x, "t": t})

    def initGrid(self):

        N = self.n
        M = self.m
        tau = self.tau
        X0 = self.x0
        XN = self.xn

        gridX = np.zeros(N + 1)
        gridT = np.zeros(M + 1)
        gridY = np.zeros((M + 1, N + 1))
        kVal = np.zeros((M + 1, N + 1))

        for i in range(N + 1):
            gridX[i] = ((XN - X0) / N) * i

        for j in range(M + 1):
            gridT[j] = tau * j
            for i in range(N + 1):
                gridY[j][i] = self.f(gridX[i], gridT[j])
                kVal[j][i] = (self.k(gridX[i], gridT[j]))
        #print(N)
        #print(gridT, gridX, gridY)
        self.kValues = kVal
        print(self.kValues)
        return (gridX, gridT, gridY)

    def fourPointDiffScheme(self, gridX, gridT, gridF):
        N = self.n
        M = self.m
        tau = self.tau
        h = self.h
        X0 = self.x0
        XN = self.xn

        ans = np.zeros((M + 1, N + 1))
        print(M, N)
        if (not self.kIsDiff):
            coeff = ((self.k(0, 0) * tau) / (h * h))
            diffC = ((h * h) - 2 * (self.k(0, 0) * tau)) / (h * h)
            print("KF: ", diffC, coeff)
        for i in range(N + 1):
            ans[0][i] = self.u(gridX[i], 0)
        for i in range(M + 1):
            ans[i][0] = self.u(X0, gridT[i])
            ans[i][N] = self.u(XN, gridT[i])
        
        for j in range(M):
            for i in range(1, N):
                if (self.kIsDiff):
                    kright = (self.kValues[j][i] + self.kValues[j][i + 1]) / 2
                    kleft = (self.kValues[j][i] + self.kValues[j][i - 1]) / 2
                    diffRight = ans[j][i + 1] - ans[j][i]
                    diffLeft = ans[j][i] - ans[j][i - 1]
                    #print(diffRight, diffLeft)
                    ans[j + 1][i] = tau * (((kright * diffRight) - (kleft * diffLeft)) / (h * h)  + gridF[j][i]) + ans[j][i]
                else:
                    ans[j + 1][i] = (coeff * ans[j][i + 1]) + diffC * ans[j][i] + coeff * ans[j][i - 1] + tau * gridF[j][i]

        #print (ans)
                
        return ans

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

    def calcMatrix(self, t):
        N = self.n
        M = self.m
        tau = self.tau
        h = self.h
        X0 = self.x0
        XN = self.xn
        alpha = np.zeros(N)
        beta = np.zeros(N)
        if (not self.kIsDiff):
            K = self.k(0, 0)
            alpha[0] = -1 * (K * tau) / (h * h)
            beta[0] = ((h * h) + 2 * (K * tau)) / (h * h)
        else:
            for i in range(1, N):
                kright = (self.kValues[t][i] + self.kValues[t][i + 1]) / 2
                kleft = (self.kValues[t][i] + self.kValues[t][i - 1]) / 2
                alpha[i] = -1 * tau * kleft / (h * h)
                beta[i] = (tau * (kleft + kright) + h * h) / (h * h)
        result = np.zeros((N, N))
        for i in range(N - 1):
            if (self.kIsDiff):
                result[i][i] = beta[i + 1]
                result[i + 1][i] = alpha[i + 1]
                result[i][i + 1] = alpha[i + 1]
                result[N - 1][N - 1] = beta[N - 1]

            else:
                result[i][i] = beta[0]
                result[i + 1][i] = alpha[0]
                result[i][i + 1] = alpha[0]
                result[N - 1][N - 1] = beta[0]
        return result

    def calcRightVector(self, y, gridF, j):
        N = self.n
        tau = self.tau
        result = np.zeros(N)
        for i in range(N):
            result[i] = y[i] + tau * gridF[j][i]
        return result

    def linSolution2(self, gridX, gridT, gridF):
        N = self.n
        M = self.m
        tau = self.tau
        ans = np.zeros((M + 1, N + 1))
        rightArr = np.zeros(N)
        coeffMatrix = self.calcMatrix(0)
        #print(coeffMatrix)
        for i in range(N):
            ans[0][i] = self.u(gridX[i], 0)
            rightArr[i] = ans[0][i]

        for j in range(1, M):
            #for i in range(N):
                #rightArr[i] = rightArr[i] - tau * gridF[j][i] * 1
            if (self.kIsDiff):
                coeffMatrix = self.calcMatrix(j)
            newLayer = self.sweepMethod(coeffMatrix, rightArr, N)

            for i in range(len(newLayer)):
                ans[j][i] = newLayer[i]
                rightArr[i] = newLayer[i]

        return ans












