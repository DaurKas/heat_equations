from cmath import sin
from doctest import FAIL_FAST
from turtle import right
from unicodedata import ucd_3_2_0
from unittest import result
from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QLabel, QLineEdit, QRadioButton, QVBoxLayout, QFrame, QPushButton, QSlider, QCheckBox, QComboBox

from PyQt5 import QtCore, QtWidgets, QtGui, QtMultimedia
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import sys
import math
import numpy as np 
from scipy.integrate import odeint

import LinearTEQ as lineq
import runge_method
import NonLinearTEQ as nolineq





class Window(QMainWindow):
    def count(self, ax1, t, x0, xn, h, tau, qsl, isRungeRule):
        
        if (self.type == "Линейная задача"):
            linHeatEq = lineq.linearEquation(x0, xn, h, tau, self.kIsDiff,
                                            self.u0LineEdit.text(),
                                            self.unLineEdit.text(),
                                            self.v0LineEdit.text(),
                                            self.kLineEdit.text(),
                                            self.fLineEdit.text() )
            (gridX, gridT, gridF) = linHeatEq.initGrid()

            
            u2 = linHeatEq.fourPointDiffScheme(gridX, gridT, gridF)
            qsl.setMaximum(int(1 / tau))
            if (isRungeRule):
                new_linHeatEq = lineq.linearEquation(x0, xn, h, 0.1, self.kIsDiff,
                                            self.u0LineEdit.text(),
                                            self.unLineEdit.text(),
                                            self.v0LineEdit.text(),
                                            self.kLineEdit.text(),
                                            self.fLineEdit.text() )
                (new_gridX, new_gridT, new_gridF) = new_linHeatEq.initGrid()
                (u, self.numOfIter, self.resultRungeTau) = runge_method.rungeMethod(new_linHeatEq, new_gridX, new_gridT, new_gridF, self.eps, 1)
                self.precisionControlSteps.setText(f"Количество шагов: {self.numOfIter}")
                self.precisionControlResultTau.setText(f"Итоговое tau: {self.resultRungeTau}")
                self.gridX = gridX
                self.u = u
                self.explicitCheck.setChecked(False)
                qsl.setMaximum(int(1 / self.resultRungeTau))
                print("HELLO THERE")
            else:
                u = linHeatEq.linSolution2(gridX, gridT, gridF)
                self.plt2 = ax1.plot(gridX, u2[t], 'g', label='u(x, t)')
            self.plt1 = ax1.plot(gridX, u[t], 'r', label='u(x, t)')
        elif (self.type == "Нелинейная задача"):
            print("HE HE")
            nonLinHeatEq = nolineq.NonlinearEquation(x0, xn, h, tau, self.kIsDiff,
                                            self.u0LineEdit.text(),
                                            self.unLineEdit.text(),
                                            self.v0LineEdit.text(),
                                            self.nonLinKLineEdit.text(),
                                            self.fLineEdit.text() )
            
            print("HI HI1")

            (gridX, gridT, gridF) = nonLinHeatEq.initGrid()
            print("HI HI")
            self.gridX = gridX
            self.u3 = nonLinHeatEq.implicitSolution(gridX, gridT, gridF)
            qsl.setMaximum(int(1 / tau))
            self.plt3 = ax1.plot(gridX, self.u3[t], 'r', label='u(x, t)')
            u = self.u3
            u2 = self.u3


   
        ax1.legend()
        ax1.grid()
        ax1.minorticks_on()
        ax1.grid(which='minor', 
            color = 'k', 
            linestyle = ':')
        ax1.set_xlim(x0, xn)
        ax1.set_ylim(0, 5)
        ax1.set(xlabel='x', ylabel='y', title="1")  
        
        return (gridX, u, u2)

    def drawPlt(self, ax1, t, gridX, x0, xn, isExpl, isImpl):

        ax1.clear()
        if (self.type == "Линейная задача"):
            if (isImpl):
                self.plt1 = ax1.plot(gridX, self.u[t], 'r', label='u(x)')
            if (isExpl):
                self.plt2 = ax1.plot(gridX, self.u2[t], 'g', label='u(x, t)')
        elif (self.type == "Нелинейная задача"):
            self.plt3 = ax1.plot(gridX, self.u3[t], 'r', label='u(x)')

        ax1.grid()
        ax1.minorticks_on()
        ax1.legend()
        ax1.grid(which='minor', 
            color = 'k', 
            linestyle = ':')
        ax1.set_xlim(x0, xn)
        ax1.set_ylim(0, 5)   
        ax1.set(xlabel='x', ylabel='y', title="1")
        return

    def changeValue(self, value):
            print(value)
            self.t = value
            t_string = str(self.t)
            self.label4.setText("t: " + t_string)

            try:
                self.drawPlt(self.ax1, value, self.gridX, self.x0, self.xn, self.isShowExplicit, self.isShowImplicit)
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

        self.ax1 = self.figure.add_subplot(111)

        self.ax1.set_xlim(0, 5)
        self.ax1.set_ylim(-5, 5)
        self.plot_widget = QWidget(self)
        self.plot_widget.setGeometry(0, 0, 800, 600)
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
        self.lineEdit3.setText("0.1")                
        self.h = float(self.lineEdit3.text())

        
        self.label3 = QLabel("tau =", self)
        self.label3.setGeometry(290, 650, 50, 20)
        self.lineEdit4 = QLineEdit(self)
        self.lineEdit4.setGeometry(310, 650, 50, 20)
        self.lineEdit4.setText("0.005")
        self.tau = float(self.lineEdit4.text())


        self.t = 0
        t_string = str(self.t)
        self.label4 = QLabel("t: " + t_string, self)
        self.label4.setGeometry(370, 650, 50, 20)
        self.qsl = QSlider(QtCore.Qt.Horizontal, self)
        self.qsl.setGeometry(420, 650, 200, 20)
        self.qsl.valueChanged[int].connect(self.changeValue)
        self.qsl.setMaximum(int(1 / self.tau))

        self.equationTypeSelect = QComboBox(self)
        self.equationTypeSelect.setGeometry(820, 100, 140, 20)
        self.equationTypeSelect.addItem("Линейная задача")
        self.equationTypeSelect.addItem("Нелинейная задача")
        self.type = self.equationTypeSelect.currentText()
        self.equationTypeSelect.currentIndexChanged[int].connect(self.changeTaskType)


        self.explicitLabel = QLabel("Явная схема", self)
        self.explicitLabel.setGeometry(820, 200, 100, 20)
        self.explicitCheck = QCheckBox(self)
        self.explicitCheck.setGeometry(900, 200, 20, 20)
        self.explicitCheck.setChecked(1)
        self.isShowExplicit = 1
        self.explicitCheck.stateChanged.connect(self.showExplicit)

        self.implicitLabel = QLabel("Неявная схема", self)
        self.implicitLabel.setGeometry(940, 200, 100, 20)
        self.implicitCheck = QCheckBox(self)
        self.implicitCheck.setGeometry(1020, 200, 20, 20)
        self.implicitCheck.setChecked(1)
        self.isShowImplicit = 1
        self.implicitCheck.stateChanged.connect(self.showImplicit)

        self.precisionControlLabel = QLabel("Авто-поиск tau", self)
        self.precisionControlLabel.setGeometry(820, 300, 100, 20)
        self.precisionControlCheck = QCheckBox(self)
        self.precisionControlCheck.setGeometry(920, 300, 20, 20)
        self.isRungeRule = False
        self.precisionControlCheck.stateChanged.connect(self.useRungeRule)

        self.epsLabel = QLabel("eps =", self)
        self.epsLabel.setGeometry(960, 300, 50, 20)

        self.epsLineEdit = QLineEdit(self)
        self.epsLineEdit.setGeometry(1000, 300, 100, 20)
        self.epsLineEdit.setText("0.01")
        self.eps = float(self.epsLineEdit.text())

        self.precisionControlSteps = QLabel("Количество шагов:", self)
        self.precisionControlSteps.setGeometry(820, 320, 200, 20)
        self.precisionControlResultTau = QLabel("Итоговое tau:", self)
        self.precisionControlResultTau.setGeometry(820, 340, 200, 20)

        self.klabel = QLabel("K(x, t) = ", self)
        self.kLineEdit = QLineEdit(self)
        self.klabel.setGeometry(820, 380, 60, 20)
        self.kLineEdit.setGeometry(890, 380, 100, 20)
        self.kLineEdit.setText("1")
        self.kTypeLabel = QLabel("К переменный", self)
        self.kTypeCheck = QCheckBox(self)
        self.kTypeLabel.setGeometry(1020, 380, 100, 20)
        self.kTypeCheck.setGeometry(1100, 380, 20, 20)
        self.kIsDiff = False
        self.kTypeCheck.stateChanged.connect(self.changeKType)

        self.nonLinKLabel = QLabel("K(u) = ", self)
        self.nonLinKLineEdit = QLineEdit(self)
        self.nonLinKLabel.setGeometry(820, 380, 60, 20)
        self.nonLinKLineEdit.setGeometry(890, 380, 100, 20)
        self.nonLinKLabel.setHidden(True)
        self.nonLinKLineEdit.setHidden(True)

        self.u0Label = QLabel("U(x0, t) = ", self)
        self.u0LineEdit = QLineEdit(self)
        self.u0Label.setGeometry(820, 400, 60, 20)
        self.u0LineEdit.setGeometry(890, 400, 100, 20)
        self.u0LineEdit.setText("0")


        self.unLabel = QLabel("U(xn, t) = ", self)
        self.unLineEdit = QLineEdit(self)
        self.unLabel.setGeometry(820, 420, 60, 20)
        self.unLineEdit.setGeometry(890, 420, 100, 20)
        self.unLineEdit.setText("0")


        self.v0Label = QLabel("U(x, 0) = ", self)
        self.v0LineEdit = QLineEdit(self)
        self.v0Label.setGeometry(820, 440, 60, 20)
        self.v0LineEdit.setGeometry(890, 440, 100, 20)
        self.v0LineEdit.setText("np.sin(x)")

        self.fLabel = QLabel("f(x, t) = ", self)
        self.fLineEdit = QLineEdit(self)
        self.fLabel.setGeometry(820, 460, 60, 20)
        self.fLineEdit.setGeometry(890, 460, 100, 20)
        self.fLineEdit.setText("np.sin(x * t)")






                       
        self.label_error = QLabel("", self)
        self.label_error.setGeometry(370, 200, 400, 30)
        
        self.button = QPushButton('Расчет', self)
        self.button.setGeometry(800, 650, 200, 40)
        self.button.clicked.connect(self.plot)
        #self.plot()
        self.show()
        
        
    
    def plot(self):
        self.ax1.clear()
        self.label_error.setText("")
        self.x0 = float(self.lineEdita1.text())
        self.xn = float(self.lineEditb1.text())
        self.h = float(self.lineEdit3.text())
        self.tau = float(self.lineEdit4.text())
        self.eps = float(self.epsLineEdit.text())

                
        try:
            (self.gridX, self.u, self.u2) = self.count(self.ax1, 1, self.x0, self.xn, self.h, self.tau, self.qsl, self.isRungeRule)
        except Exception as e:
            self.label_error.setText(str(e))
        
        #plt.tight_layout()
        #plt.gcf().subplots_adjust(right=0.9, left=0.1, bottom=0.12, top=0.9)
        plt.gcf().subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.13)

        self.canvas.draw()
    
    def showExplicit(self, int):
        if (self.isShowExplicit):
            self.isShowExplicit = 0
        else:
            self.isShowExplicit = 1
        self.drawPlt(self.ax1, self.t, self.gridX, self.x0, self.xn, self.isShowExplicit, self.isShowImplicit)
        self.canvas.draw()

    def showImplicit(self, int):
        if (self.isShowImplicit):
            self.isShowImplicit = 0
        else:
            self.isShowImplicit = 1
        self.drawPlt(self.ax1, self.t, self.gridX, self.x0, self.xn, self.isShowExplicit, self.isShowImplicit)
        self.canvas.draw()

    def useRungeRule(self, int):
        self.isRungeRule = not(self.isRungeRule)

    def changeKType(self, int):
        self.kIsDiff = not(self.kIsDiff)

    def changeTaskType(self, int):
        self.type = self.equationTypeSelect.currentText()
        if (self.type == "Нелинейная задача"):
            self.nonLinKLabel.setHidden(False)
            self.nonLinKLineEdit.setHidden(False)
            self.explicitLabel.setHidden(True)
            self.explicitCheck.setHidden(True)
            self.implicitLabel.setHidden(True)
            self.implicitCheck.setHidden(True)
            self.precisionControlLabel.setHidden(True)
            self.precisionControlCheck.setHidden(True)
            self.epsLabel.setHidden(True)
            self.epsLineEdit.setHidden(True)
            self.precisionControlSteps.setHidden(True)
            self.precisionControlResultTau.setHidden(True)
            self.klabel.setHidden(True)
            self.kLineEdit.setHidden(True)
            self.kTypeLabel.setHidden(True)
            self.kTypeCheck.setHidden(True)
            #self.u0Label.setHidden(True)
            #self.u0LineEdit.setHidden(True)
            #self.unLabel.setHidden(True)
            #self.unLineEdit.setHidden(True)
            #self.v0Label.setHidden(True)
            #self.v0LineEdit.setHidden(True)
            #self.fLabel.setHidden(True)
            #self.fLineEdit.setHidden(True)
            #self.fLineEdit.setHidden(True)
        elif (self.type == "Линейная задача"):
            self.nonLinKLabel.setHidden(True)
            self.nonLinKLineEdit.setHidden(True)
            self.explicitLabel.setHidden(False)
            self.explicitCheck.setHidden(False)
            self.implicitLabel.setHidden(False)
            self.implicitCheck.setHidden(False)
            self.precisionControlLabel.setHidden(False)
            self.precisionControlCheck.setHidden(False)
            self.epsLabel.setHidden(False)
            self.epsLineEdit.setHidden(False)
            self.precisionControlSteps.setHidden(False)
            self.precisionControlResultTau.setHidden(False)
            self.klabel.setHidden(False)
            self.kLineEdit.setHidden(False)
            self.kTypeLabel.setHidden(False)
            self.kTypeCheck.setHidden(False)
            #self.u0Label.setHidden(False)
            #self.u0LineEdit.setHidden(False)
            #self.unLabel.setHidden(False)
            #self.unLineEdit.setHidden(False)
            #self.v0Label.setHidden(False)
            #self.v0LineEdit.setHidden(False)
            #self.fLabel.setHidden(False)
            #self.fLineEdit.setHidden(False)
            #self.fLineEdit.setHidden(False)







if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = Window()
    app.exec_()
    #main()
