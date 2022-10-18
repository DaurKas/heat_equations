from cmath import sin
from turtle import right
from unittest import result
from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QLabel, QLineEdit, QRadioButton, QVBoxLayout, QFrame, QPushButton, QSlider

from PyQt5 import QtCore, QtWidgets, QtGui, QtMultimedia
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import sys
import math
import numpy as np 
from scipy.integrate import odeint

K = 1
a = 1
tau = 0.00005
h = 0.01
X0 = 0
XN = 3.14
T0 = 0
TN = 1
M = int((TN - T0) / tau)
N = int((XN - X0) / h)
def v0(x, t):
    return 0

def vn(x, t):
    return 0

def u0(x, t):
    return np.sin (x)

def u(xi, tj):
    if (xi == X0):
        return v0(xi, tj)
    if (xi == XN):
        return vn(xi, tj)
    if (tj == 0):
        return u0(xi, tj)
    return 1

def f(x, t):
    return np.cos(t * x)

def initGrid():
    gridX = np.zeros(N + 1)
    gridT = np.zeros(M + 1)
    gridY = np.zeros((M + 1, N + 1))
    print(M, N)
    for j in range(M + 1):
        gridT[j] = tau * j
        for i in range(N + 1):
            gridX[i] = ((XN - X0) / N) * i
            #print(i)
            gridY[j][i] = u(gridX[i], gridT[j])
    return (gridX, gridT, gridY)

def initMatrix(gridX, gridF):
    matrix = np.zeros((N + 1, N + 1))
    return matrix

def fourPointDiffScheme(gridX, gridT, gridF):
    ans = np.zeros((M + 1, N + 1))
    print(M, N)
    coeff = ((K * tau) / (h * h))
    diffC = ((h * h) - 2 * (K * tau)) / (h * h)
    print("KF: ", diffC, coeff)
    for i in range(N):
        ans[0][i] = u(gridX[i], 0)
    
    for j in range(M - 1):
        for i in range(N - 1):
            ans[j + 1][i] = (coeff * ans[j][i + 1]) + diffC * ans[j][i] + coeff * ans[j][i - 1] + tau * h * h * gridF[j][i]

            
    return ans

def sweepMethod(leftMatrix, rightArr, N):
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
        #print("CHECK: ", a[i], b[i])

    result[N - 1] = (rightArr[N - 1] - leftMatrix[N - 1][N - 2] * b[N - 2]) / (leftMatrix[N - 1][N - 1] + leftMatrix[N - 1][N - 2] * a[N - 2])
    
    for i in range(N - 2, -1, -1):
        result[i] = a[i] * result[i + 1] + b[i]

    return result

def calcMatrix():
    alpha = -1 * (K * tau) / (h * h)
    beta = ((h * h) + 2 * (K * tau)) / (h * h)
    result = np.zeros((N, N))
    for i in range(N - 1):
            result[i][i] = beta
            result[i + 1][i] = alpha
            result[i][i + 1] = alpha
    result[N - 1][N - 1] = beta
    return result

def calcRightVector(y, gridF, j):
    result = np.zeros(N)
    for i in range(N):
        result[i] = y[i] + tau * gridF[j][i]
    return result

def linSolution2(gridX, gridT, gridF):
    ans = np.zeros((M + 1, N + 1))
    rightArr = np.zeros(N)
    #coeffMatrix = np.zeros((N, N))
    coeffMatrix = calcMatrix()
    print(coeffMatrix)
    for i in range(N):
        ans[0][i] = u(gridX[i], 0)
        rightArr[i] = ans[0][i]
    for j in range(1, M):
        #rightArr = calcRightVector(rightArr, gridF, j)
        for i in range(N):
            rightArr[i] = rightArr[i] - tau * gridF[j][i] * 0
        newLayer = sweepMethod(coeffMatrix, rightArr, N)
        for i in range(len(newLayer)):
            ans[j][i] = newLayer[i]
            rightArr[i] = newLayer[i]
        #print(rightArr)
    return ans

def solveEquation():
    y = np.zeros(N)
    return y

def writeArray(nparray):
    f = open('results.txt', 'w')
    for i in range(len(nparray)):
        for j in range(len(nparray[i])):
            f.write(str(nparray[i][j]))
        f.write('\n')

def count(ax1, t):
    
    (gridX, gridT, gridF) = initGrid()
    print(gridX)
    #u = fourPointDiffScheme(gridX, gridT, gridF)
    u = linSolution2(gridX, gridT, gridF)
    u2 = fourPointDiffScheme(gridX, gridT, gridF)
    #print(u[t])
    plt1 = ax1.plot(gridX, u[t], 'r', label='u(x, t)')
    plt2 = ax1.plot(gridX, u2[t], 'g', label='u(x, t)')
    #ax1.axhline(y=u[1], color='black', linestyle = ':')
    ax1.legend()
    ax1.grid()
    ax1.set_xlim(0, 3.14)
    ax1.set_ylim(0, 2)
    ax1.set(xlabel='x', ylabel='y', title="1")  
    
    return (gridX, u, u2, plt1)

def drawPlt(ax1, t, gridX, u, u2, plt1):

    ax1.clear()
    ax1.plot(gridX, u[t], 'r', label='u(x)')
    ax1.plot(gridX, u2[t], 'g', label='u(x, t)')
    #(u[t])
    #ax1.axhline(y=u[1], color='black', linestyle = ':')
    ax1.grid()
    ax1.legend()
    ax1.set_xlim(0, 3.14)
    ax1.set_ylim(0, 2)    
    ax1.set(xlabel='x', ylabel='y', title="1")
    return

class Window(QMainWindow):
    def changeValue(self, value):
            print(value)
            self.t = value
            t_string = str(self.t)
            self.label4.setText("t: " + t_string)
            try:
                drawPlt(self.ax1, value, self.gridX, self.u, self.u2, self.plt1)
            except Exception as e:
                self.label_error.setText(str(e))
            #self.ax1.clear()
            #self.plot()
            #self.show()
            self.canvas.draw()

    def __init__(self):
        super().__init__()
        self.setWindowTitle('Уравнение теплопроводности (диффузии?)')
        self.setGeometry(150, 150, 1200, 720)
        self.setFixedWidth(1200)
        self.setFixedHeight(720)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)

        self.ax1 = self.figure.add_subplot(121)
        self.ax1.set_xlim(0, 5)
        self.ax1.set_ylim(-5, 5)
        #self.ax = plt.axes(xlim = (0, 5), ylim = (-5, 5))
        self.plot_widget = QWidget(self)
        self.plot_widget.setGeometry(0, 0, 1200, 600)
        plot_box = QVBoxLayout()
        plot_box.addWidget(self.canvas)
        self.plot_widget.setLayout(plot_box)

        self.label0 = QLabel("X1 =",self)
        self.label0.setGeometry(50,650, 50, 20)
        self.lineEdita1 = QLineEdit(self)
        self.lineEdita1.setGeometry(70, 650, 50, 20)
        self.lineEdita1.setText("0")
        
        self.label1 = QLabel("XN =", self)
        self.label1.setGeometry(130, 650, 50, 20)
        self.lineEditb1 = QLineEdit(self)
        self.lineEditb1.setGeometry(150, 650, 50, 20)
        self.lineEditb1.setText("3.14")
        
        self.label2 = QLabel("h =", self)
        self.label2.setGeometry(210, 650, 60, 20)
        self.lineEdit3 = QLineEdit(self)
        self.lineEdit3.setGeometry(230, 650, 50, 20)
        self.lineEdit3.setText("0.01")
        
        self.label3 = QLabel("tau =", self)
        self.label3.setGeometry(290, 650, 50, 20)
        self.lineEdit4 = QLineEdit(self)
        self.lineEdit4.setGeometry(310, 650, 50, 20)
        self.lineEdit4.setText("0.00005")

        self.t = 0
        t_string = str(self.t)
        self.label4 = QLabel("t: " + t_string, self)
        self.label4.setGeometry(370, 650, 50, 20)
        self.qsl = QSlider(QtCore.Qt.Horizontal, self)
        self.qsl.setGeometry(420, 650, 200, 20)
        self.qsl.valueChanged[int].connect(self.changeValue)
        self.qsl.setMaximum(20000)

        self.condLabel = QLabel("U_t = K * U_xx + f(x, t)\n U(x, 0) = sin(x)\n U(x_0, t) = 0\n U(x_n, t) = 0")
        self.condLabel.setGeometry(700, 650, 100, 100)
               
        self.label_error = QLabel("", self)
        self.label_error.setGeometry(370, 200, 400, 30)
        
        #self.button = QPushButton('Начать расчет', self)
        #self.button.setGeometry(800, 650, 200, 40)
        #self.button.clicked.connect(self.plot)
        self.plot()
        self.show()
        
        
    
    def plot(self):
        self.ax1.clear()
        self.label_error.setText("")
        
        if len(self.lineEdita1.text())!=0:
            self.a=float(self.lineEdita1.text())
        else:
            self.a = -0.02
            self.lineEdita1.setText("-0.02")
            
        if len(self.lineEditb1.text()) != 0:
            self.b=float(self.lineEditb1.text())
        else:
            self.b = 0.25
            self.lineEditb1.setText("0.25")
            
        if len(self.lineEdit3.text()) != 0:
            self.c=float(self.lineEdit3.text())
        else:
            self.c = 1.25
            self.lineEdit3.setText("1.25")
            
        if len(self.lineEdit4.text()) != 0:  
            self.dt=float(self.lineEdit4.text())
        else:
            self.dt = 0.01
            self.lineEdit4.setText("0.01")
                
        try:
            (self.gridX, self.u, self.u2, self.plt1) = count(self.ax1, 1)
        except Exception as e:
            self.label_error.setText(str(e))
        
        plt.tight_layout()
        plt.gcf().subplots_adjust(bottom=0.12, top=0.9, wspace=0.25)
        self.canvas.draw()



def main():
    (gridX, gridT, gridF) = initGrid()
    print(gridX)
    print(gridT)
    print(gridF)
    #writeArray(fourPointDiffScheme(gridX, gridT, gridF))
    u = linSolution2(gridX, gridT, gridF)
    #atrix = [[1, 2, 0], [2, 1, 2], [0, 2, 0]]
    #vec = [1, 2, 3]
    #sol = sweepMethod(matrix, vec, 3)
    #print(sol)
    writeArray(u)
    #print(u)
    #print(fourPointDiffScheme(gridX, gridT, gridF))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = Window()
    app.exec_()
    #main()

