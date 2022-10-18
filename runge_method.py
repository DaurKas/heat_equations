import numpy as np 
import LinearTEQ

EPS = 0.01

def quadDiff(u, v):
    diff = 0
    for i in range(len(u)):
        diff += (u[i] - v[i])**2
    return diff

def rungeMethod(linEq, gridX, gridT, gridF):
    u = linEq.linSolution2(gridX, gridT, gridF)
    linEq2 = LinearTEQ.linearEquation(linEq.x0, linEq.xn, linEq.h, (linEq.tau) / 2)
    (gridX2, gridT2, gridF2) = linEq2.initGrid()
    u2 = linEq2.linSolution2(gridX2, gridT2, gridF2)
    if (quadDiff(u[1], u2[1]) < EPS * EPS):
        return linEq
    return rungeMethod(linEq2, gridX2, gridT2, gridF2)




