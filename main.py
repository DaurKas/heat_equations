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

import LinearTEQ as lineq

def count(ax1, t, x0, xn, h, tau):
    
    linHeatEq = lineq.linearEquation(x0, xn, h, tau)
    (gridX, gridT, gridF) = linHeatEq.initGrid()
    u = linHeatEq.linSolution2(gridX, gridT, gridF)
    u2 = linHeatEq.fourPointDiffScheme(gridX, gridT, gridF)
    plt1 = ax1.plot(gridX, u[t], 'r', label='u(x, t)')
    plt2 = ax1.plot(gridX, u2[t], 'g', label='u(x, t)')
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
        self.x0 = float(self.lineEdita1.text())
        
        self.label1 = QLabel("XN =", self)
        self.label1.setGeometry(130, 650, 50, 20)
        self.lineEditb1 = QLineEdit(self)
        self.lineEditb1.setGeometry(150, 650, 50, 20)
        self.lineEditb1.setText("3.14")
        self.xn = float(self.lineEditb1.text())

        
        self.label2 = QLabel("h =", self)
        self.label2.setGeometry(210, 650, 60, 20)
        self.lineEdit3 = QLineEdit(self)
        self.lineEdit3.setGeometry(230, 650, 50, 20)
        self.lineEdit3.setText("0.01")                
        self.h = float(self.lineEdit3.text())

        
        self.label3 = QLabel("tau =", self)
        self.label3.setGeometry(290, 650, 50, 20)
        self.lineEdit4 = QLineEdit(self)
        self.lineEdit4.setGeometry(310, 650, 50, 20)
        self.lineEdit4.setText("0.00005")
        self.tau = float(self.lineEdit4.text())


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
        
        self.button = QPushButton('Расчет', self)
        self.button.setGeometry(800, 650, 200, 40)
        self.button.clicked.connect(self.plot)
        self.plot()
        self.show()
        
        
    
    def plot(self):
        self.ax1.clear()
        self.label_error.setText("")
        self.x0 = float(self.lineEdita1.text())
        self.xn = float(self.lineEditb1.text())
        self.h = float(self.lineEdit3.text())
        self.tau = float(self.lineEdit4.text())

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
            (self.gridX, self.u, self.u2, self.plt1) = count(self.ax1, 1, self.x0, self.xn, self.h, self.tau)
        except Exception as e:
            self.label_error.setText(str(e))
        
        plt.tight_layout()
        plt.gcf().subplots_adjust(bottom=0.12, top=0.9, wspace=0.25)
        self.canvas.draw()





if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = Window()
    app.exec_()
    #main()
