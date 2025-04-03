#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright (c) 2022-2025 MRTech SK, s.r.o.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

###############################################################
# Application name: ldc
# Algebraic Lens Distortion Model Estimation.
###############################################################
import tkinter as tk
from tkinter import *  
from tkinter import filedialog
from tkinter import messagebox
import math
from PIL import Image, ImageTk
import os
import sys, getopt
import numpy as np
import cv2
import time
import random
import json

# Structures and Classes
#---------------------------------------------------
class CPoint:
    def __init__(self, x = 0, y = 0):
        self.x = x
        self.y = y
    def Norm(self, next):
        paso = (self.x - next.x) * (self.x - next.x) + (self.y - next.y) * (self.y - next.y)
        if paso > 0:
            return math.sqrt(paso)     
        else:
            return 0.0  

class GraphicalObjects:
    def __init__(self):
        self.grid_corner_circle = None
        self.grid_corner_v_line = None
        self.grid_corner_h_line = None
        self.grid_corner_correct_line = None

class CConfiguration():
    def __init__(self):
        self.output_file = "lens_profile.json"
        self.output_file_set = False 
        self.target_squares_rows = 19
        self.target_squares_columns = 19

        self.lable_font = "Arial"
        self.label_font_size = 14

        self.config_file = None

        self.target_searching_timeout = 10 #s

        self.optimize_center = False

        self.max_image_proc_resolution = 1024.

    def ParseConfigFile(self, config_file=None):
        self.config_file = config_file
        if self.config_file != None:
            try:
                with open(self.config_file, 'r+') as f:
                    file_content = f.read()
                    content = json.loads(file_content)
                    if "output file" in content:
                        if self.output_file_set == False:
                            self.output_file = content["output file"]    
                    if "target squares rows" in content:
                        self.target_squares_rows = int(content["target squares rows"]) - 1 
                    if "target squares columns" in content:
                        self.target_squares_columns = int(content["target squares columns"]) - 1 
                    if "label font" in content:
                        self.label_font = content["label font"]  
                    if "label font size" in content:
                        self.label_font_size = int(content["label font size"])      
                    if "target searching timeout" in content:
                        self.target_searching_timeout = int(content["target searching timeout"]) 
                    if "optimize center distortion" in content:
                        self.optimize_center = int(content["optimize center distortion"])         
                    if "max image processing resolution" in content:
                        self.max_image_proc_resolution = float(content["max image processing resolution"])    
                    return 0    
            except OSError as error: 
                messagebox.showwarning(title='Error', message=error)
                return 1
        return 1

#---------------------------------------------------
class CApp():
    def __init__(self, root):
        self.root = root
        self.root.geometry('%dx%d+0+0' % (int(self.root.winfo_screenwidth()), int(self.root.winfo_screenheight())))
        self.root.configure(bg="Black")
        self.root.title("Lens Distortion Correction")
   
        self.label = Label(self.root, text='Application initialization', anchor=CENTER, font=("Arial", 14), height=4, bg="Black", fg="White")
        self.label.grid(row=0, sticky="nsew")   
        self.canvas = tk.Canvas(self.root, bg="Black", width=int(self.root.winfo_screenwidth()), height=(int(self.root.winfo_screenheight())-self.label.winfo_height()))
        self.canvas.grid(row=1, sticky="nsew")
        self.lable_height = 0

        self.Input_image_file = None
        self.image_id = None
        self.proportion = 1.0

        self.InputImageTk = None
        self.resizedImage = None
        self.resizedImageTk = None
        self.background_id = None
        self.UndistoredImage = None
        self.DistoredImage = None
        self.UndistoredImageChanged = False

        self.InputImageOpenCV = None
        self.ImageOpenCVgray = None
        self.OpenCVfoundCorners = None

        self.config = None
        self.max_grid_rows = 20
        self.max_grid_columns = 20
        self.points_in_line = None
        self.lines = 0
        self.grid_corners = None
        self.grid_corners_ = None
        self.graphical_grid = None
        self.grid_corner_correct_line_color = None
        self.grid_corners_defined = False
        self.grid_point_pointed = False
        self.grid_point_pointed_coordinates = CPoint()
        self.max_grade_polinom = 4 # GRADE OF THE LENS DISTORTION POLYNOM
        self.center_x = 0.0
        self.center_y = 0.0
        self.grid_was_reset = False
        self.calculation_done = False
        self.resize_coefficient = 1.0

 
# Init
#---------------------------------------------------
root = tk.Tk()
app = CApp(root)
root.update_idletasks()
app.label.update_idletasks()
app.lable_height = app.label.winfo_height() - 4
app.canvas.update_idletasks()
geometry = root.geometry()
root_width = int(geometry[0:geometry.index("x")])
root_height = int(geometry[geometry.index("x")+1:geometry.index("+")])
app.canvas.config(width=root_width-4, height=root_height-app.label.winfo_height()-4)
app.canvas.update_idletasks()
canvas_width = app.canvas.winfo_width()
canvas_height = app.canvas.winfo_height() - app.lable_height

#---------------------------------------------------
# lens distortion correction functions
#---------------------------------------------------
#function to compute the real roots of a cubic polynomial. It returns the number of roots found sorted by magnitud
def ami_RootCubicPolynomial(a, N, x):
    if N!=3 or a[3]==0:
        return(-100000)

    a1 = a[2] / a[3]
    a2 = a[1] / a[3]
    a3 = a[0] / a[3]
    Q = (3 * a2 - a1 * a1) / 9.
    R = (9 * a1 * a2 - 27 * a3 - 2 * a1 * a1 * a1) / 54.
    D = Q * Q * Q + R * R

    if D > 0:
        S = R + math.sqrt(D)
        T = R - math.sqrt(D)
        if S > 0:
            S = math.pow(S, 1./3.)
        else:
            S = -math.pow(-S, 1./3.)
        if T > 0:
            T = math.pow(T, 1./3.)
        else:
            T = -math.pow(-T, 1./3.)
        x[0] = S + T - a1 / 3.
        return 1
    else:
        PI2 = math.acos(-1.)
        if Q != 0:
            A = math.acos(R / math.sqrt(-Q * Q * Q))
        else:
            A = 0

        Q = 2. * math.sqrt(-Q)
        x[0] = Q * math.cos(A / 3.) - a1 / 3.
        x[1] = Q * math.cos(A / 3. + 2. * PI2 / 3.) - a1 / 3.
        x[2] = Q * math.cos(A / 3 + 4. * PI2 / 3.) - a1 / 3.

        if math.fabs(x[0]) > math.fabs(x[1]):
            Q = x[1]
            x[1] = x[0]
            x[0] = Q
        if math.fabs(x[0]) > math.fabs(x[2]):
            Q = x[2]
            x[2] = x[0]
            x[0] = Q
        if math.fabs(x[1]) > math.fabs(x[2]):
            Q = x[2]
            x[2] = x[1]
            x[1] = Q

        return 3

#---------------------------------------------------
def test_compatibility_lens_distortion_model(a, Na, max_radius):
    b = [0.0 for i in range(Na + 1)]        
    # BUILD THE DERIVATIVE OF THE POLYNOMIAL
    for k in range(1, Na + 1):
        b[k] = (k + 1) * a[k]
    for r in range(1, int(max_radius + 0.5)):
        if ami_polynomial_evaluation(b, Na, r) < 0:
            return -1.
    return 0

#---------------------------------------------------
#function to add the information of a line point sequence to the 4 degree polynomial to compute the lens distortion model
def ami_lens_distortion_polynomial_update(coefficients, grid_corners, row, Np, x0, y0, k, pol):
    A = [0.0 for i in range(Np)]
    x2 = [0.0 for i in range(Np)]
    y2 = [0.0 for i in range(Np)]
    d = [0.0 for i in range(Np)]

    pol1 = [0.0 for i in range(5)]
    pol2 = [0.0 for i in range(5)]
    pol3 = [0.0 for i in range(5)]

    #COMPUTE THE DISTANCE TO THE IMAGE CENTER
    for i in range(Np):
        d[i] = math.sqrt((grid_corners[row][i][0] - x0) * (grid_corners[row][i][0] - x0) + (grid_corners[row][i][1] - y0) * (grid_corners[row][i][1] - y0))
    #COMPUTE THE POINT TRANSFORMATION WITH THE CURRENT LENS DISTORTION MODEL
    for i in range(Np):
        A[i] = ami_polynomial_evaluation(coefficients, app.max_grade_polinom, d[i])
        x2[i] = x0 + (grid_corners[row][i][0] - x0) * A[i]
        y2[i] = y0 + (grid_corners[row][i][1] - y0) * A[i]   

    #COMPUTE THE DISTANCE POWER k (THE COEFFICIENT OF THE LENS DISTORTION MODEL TO BE UPDATED     
    for i in range(Np):
        d[i] = math.pow(d[i], k)
    #COMPUTE THE VARIANCE OF THE TRANSFORMED POINTS
    x2_m = 0.0
    for i in range(Np):
        x2_m += x2[i]
    x2_m /= Np
    s_xx = 0.0
    for i in range(Np):
        s_xx += (x2[i] - x2_m) * (x2[i] - x2_m)
    s_xx /= Np
    y2_m = 0.0
    for i in range(Np):
        y2_m += y2[i]
    y2_m /= Np
    s_yy = 0.0
    for i in range(Np):
        s_yy += (y2[i] - y2_m) * (y2[i] - y2_m)
    s_yy /= Np
    #COMPUTE SOME AVERAGES WE NEED
    xA_m = 0.0
    for i in range(Np):
        xA_m += (grid_corners[row][i][0] - x0) * A[i]
    xA_m /= Np
    xd_m = 0
    for i in range(Np):
        xd_m += (grid_corners[row][i][0] - x0) * d[i]
    xd_m /= Np
    yA_m = 0
    for i in range(Np):
        yA_m += (grid_corners[row][i][1] - y0) * A[i]
    yA_m /= Np
    yd_m = 0
    for i in range(Np):
        yd_m += (grid_corners[row][i][1] - y0) * d[i]
    yd_m /= Np  
    #COMPUTE THE POLYNOMIAL TO MINIMIZE
    for i in range(5):
        pol1[i] = 0.0
        pol2[i] = 0.0
        pol3[i] = 0.0
    for i in range(Np):
        pol1[0] += ((grid_corners[row][i][0] - x0) * A[i] - xA_m) * ((grid_corners[row][i][0] - x0) * A[i] - xA_m)
        pol1[1] += 2. * ((grid_corners[row][i][0]-x0) * A[i] - xA_m)  *((grid_corners[row][i][0] - x0)*  d[i] - xd_m)
        pol1[2] += ((grid_corners[row][i][0] - x0) * d[i] - xd_m) * ((grid_corners[row][i][0] - x0) * d[i] - xd_m)
        pol2[0] += ((grid_corners[row][i][1] - y0) * A[i] - yA_m) * ((grid_corners[row][i][1]  -y0) * A[i] - yA_m)
        pol2[1] += 2. * ((grid_corners[row][i][1] - y0) * A[i] - yA_m) * ((grid_corners[row][i][1] - y0) * d[i] - yd_m)
        pol2[2] += ((grid_corners[row][i][1] - y0) * d[i] - yd_m) * ((grid_corners[row][i][1] - y0) * d[i] - yd_m)
        pol3[0] += ((grid_corners[row][i][1] - y0) * A[i] - yA_m) * ((grid_corners[row][i][0] - x0) * A[i] - xA_m)
        pol3[1] += ((grid_corners[row][i][1] - y0) * A[i] - yA_m) * ((grid_corners[row][i][0] - x0) * d[i] - xd_m) + ((grid_corners[row][i][1] - y0) * d[i] - yd_m) * ((grid_corners[row][i][0] - x0) * A[i] - xA_m)
        pol3[2] += ((grid_corners[row][i][1] - y0) * d[i] - yd_m) * ((grid_corners[row][i][0] - x0) * d[i] - xd_m)

    for i in range(3):
        for j in range(3):
            pol[i + j] += (pol1[i] * pol2[j] - pol3[i] * pol3[j]) / 1.

#---------------------------------------------------
def ami_lens_distortion_model_update(coefficients, grid_corners, k, pol, max_radius):
    M = app.max_grade_polinom
    p = [0.0 for i in range(3)]             
    x = [0.0 for i in range(3)]         
    b = [0.0 for i in range(4)]         
    b2 = [0.0 for i in range(app.max_grade_polinom + 8)]      
    #FILL THE AUXILIARY LENS DISTORTION MODEL   
    for i in range(M):
        b2[i] = coefficients[i]

    b[0] = pol[1]
    b[1] = 2 * pol[2]
    b[2] = 3. * pol[3]
    b[3] = 4. * pol[4]   

    M = ami_RootCubicPolynomial(b, 3, x)  

    for i in range(M):
        p[i] = ami_polynomial_evaluation(pol, 4, x[i]) 
    Error = 1e30 
    j = M  
    for i in range(M):
        b2[k] = coefficients[k] + x[i]
        if test_compatibility_lens_distortion_model(b2, app.max_grade_polinom, max_radius) == 0:
            if p[i] < Error:
                j = i
                Error = p[i]
    if j < M:
        coefficients[k] += x[j]   

    return 0 

#---------------------------------------------------
#function to multiply polynoms of 2 variables
def ami_2v_polynom_multiplication(p1, N1, p2, N2, p3):
    for i in range(N1 + 1): 
        for j in range(N1 + 1): 
            if p1[i][j] != 0:
                for k in range(N2 + 1):
                    for l in range(N2 + 1):
                        if p2[k][l] != 0:
                            p3[i + k][j + l] += p1[i][j] * p2[k][l]


#---------------------------------------------------
#function To add the information of a line point sequence to the 4 degree polynomial to compute the lens distortion model
def ami_lens_distortion_polynomial_update_distance_2v(coefficients, grid_corners, row, Np, x0, y0, k1, k2, pol, alpha): 
    A = [0.0 for i in range(Np)]        
    x2 = [0.0 for i in range(Np)]        
    y2 = [0.0 for i in range(Np)]        
    d1 = [0.0 for i in range(Np)]        
    d2 = [0.0 for i in range(Np)]   
    pol1 = [[0.0 for i in range(2)] for j in range(2)]
    pol2 = [[0.0 for i in range(2)] for j in range(2)]
    p_xx = [[0.0 for i in range(3)] for j in range(3)]
    p_xy = [[0.0 for i in range(3)] for j in range(3)]
    p_yy = [[0.0 for i in range(3)] for j in range(3)]  

    #COMPUTE THE DISTANCE TO THE IMAGE CENTER    
    for i in range(Np): 
        d1[i] = math.sqrt((grid_corners[row][i][0] - x0) * (grid_corners[row][i][0] - x0) + (grid_corners[row][i][1] - y0) * (grid_corners[row][i][1] - y0))
    #COMPUTE THE POINT TRANSFORMATION WITH THE CURRENT LENS DISTORTION MODEL                         
    for i in range(Np): 
        A[i] = ami_polynomial_evaluation(coefficients, app.max_grade_polinom, d1[i])
        x2[i] = x0 + (grid_corners[row][i][0] - x0) * A[i]
        y2[i] = y0 + (grid_corners[row][i][1] - y0) * A[i]
    #COMPUTE THE POLYNOMS CORRESPONDING TO THE DISTANCE ERROR
    for i in range(3): 
        for j in range(3):
            p_xx[i][j] = 0. 
    for i in range(Np): 
        aux = 0.0
        for k in range(app.max_grade_polinom + 1): 
            if (k != k1 and k != k2):
                aux += coefficients[k] * math.pow(d1[i], k)
        pol1[0][0] = aux * d1[i]
        pol1[1][0] = math.pow(d1[i], k1 + 1)
        pol1[0][1] = math.pow(d1[i], k2 + 1)
        ami_2v_polynom_multiplication(pol1, 1, pol1, 1, p_xx)
    for i in range(3): 
        for j in range(3):
            p_xx[i][j] = alpha * p_xx[i][j] / Np    
    #UPDATE THE ERROR POLYNOM
    ami_2v_polynom_multiplication(p_xx, 2, p_xx, 2, pol)            

#---------------------------------------------------
#function To add the information of a line point sequence to the 4 degree polynomial to compute the lens distortion model
def ami_lens_distortion_polynomial_update_2v(coefficients, grid_corners, row, Np, x0, y0, k1, k2, pol):  
    A = [0.0 for i in range(Np)]        
    x2 = [0.0 for i in range(Np)]        
    y2 = [0.0 for i in range(Np)]        
    d1 = [0.0 for i in range(Np)]        
    d2 = [0.0 for i in range(Np)] 
    pol1 = [[0.0 for i in range(2)] for j in range(2)]
    pol2 = [[0.0 for i in range(2)] for j in range(2)]
    p_xx = [[0.0 for i in range(3)] for j in range(3)]
    p_xy = [[0.0 for i in range(3)] for j in range(3)]
    p_yy = [[0.0 for i in range(3)] for j in range(3)]      

    #COMPUTE THE DISTANCE TO THE IMAGE CENTER
    for i in range(Np):
        d1[i] = math.sqrt((grid_corners[row][i][0] - x0) * (grid_corners[row][i][0] - x0) + (grid_corners[row][i][1] - y0) * (grid_corners[row][i][1] - y0))
    #COMPUTE THE POINT TRANSFORMATION WITH THE CURRENT LENS DISTORTION MODEL
    for i in range(Np):
        A[i] = ami_polynomial_evaluation(coefficients, app.max_grade_polinom, d1[i])
        x2[i] = x0 + (grid_corners[row][i][0] - x0) * A[i]
        y2[i] = y0 + (grid_corners[row][i][1] - y0) * A[i]       
    #COMPUTE THE DISTANCE POWER k1 AND k2 (THE COEFFICIENT OF THE LENS DISTORTION MODEL TO BE UPDATED    
    for i in range(Np): 
        paso = d1[i]
        d1[i] = math.pow(paso, k1)
        d2[i] = math.pow(paso, k2)
    #COMPUTE THE VARIANCE OF THE TRANSFORMED POINTS
    x2_m = 0.0
    for i in range(Np): 
        x2_m += x2[i]
    x2_m /= Np
    s_xx = 0.0
    for i in range(Np): 
        s_xx += (x2[i] - x2_m) * (x2[i] - x2_m)
    s_xx /= Np
    y2_m = 0.0
    for i in range(Np): 
        y2_m += y2[i]
    y2_m /= Np
    s_yy = 0.0
    for i in range(Np): 
        s_yy += (y2[i] - y2_m) * (y2[i] - y2_m)
    s_yy /= Np
    #COMPUTE SOME AVERAGES WE NEED    
    xA_m = 0.0
    for i in range(Np):
        xA_m += (grid_corners[row][i][0] - x0) * A[i]
    xA_m /= Np
    xd1_m = 0.0
    for i in range(Np):
        xd1_m += (grid_corners[row][i][0] - x0) * d1[i]
    xd1_m /= Np
    xd2_m = 0.0
    for i in range(Np):
        xd2_m += (grid_corners[row][i][0] - x0) * d2[i]
    xd2_m /= Np
    yA_m = 0.0
    for i in range(Np):
        yA_m += (grid_corners[row][i][1] - y0) * A[i]
    yA_m /= Np
    yd1_m = 0.0
    for i in range(Np):
        yd1_m += (grid_corners[row][i][1] - y0) * d1[i]
    yd1_m /= Np
    yd2_m = 0
    for i in range(Np):
        yd2_m += (grid_corners[row][i][1] - y0) * d2[i]
    yd2_m /= Np

    #COMPUTE THE POLYNOMS OF THE SECOND ORDER MOMENT OF THE POINT p_xx p_xy AND p_yy DISTRIBUTION
    for i in range(Np):
        pol1[0][0]= (grid_corners[row][i][0] - x0) * A[i] - xA_m
        pol1[1][0]= (grid_corners[row][i][0] - x0) * d1[i] - xd1_m
        pol1[0][1]= (grid_corners[row][i][0] - x0) * d2[i] - xd2_m
        pol2[0][0]= (grid_corners[row][i][1] - y0) * A[i] - yA_m
        pol2[1][0]= (grid_corners[row][i][1] - y0) * d1[i] - yd1_m
        pol2[0][1]= (grid_corners[row][i][1] - y0) * d2[i] - yd2_m
        ami_2v_polynom_multiplication(pol1, 1, pol1, 1, p_xx)
        ami_2v_polynom_multiplication(pol1, 1, pol2, 1, p_xy)
        ami_2v_polynom_multiplication(pol2, 1, pol2, 1, p_yy)
    for i in range(3):    
        for j in range(3):
            p_xx[i][j] /= 1.
    ami_2v_polynom_multiplication(p_xx, 2, p_yy, 2, pol)
    for i in range(3):    
        for j in range(3):
             p_xx[i][j] = -p_xy[i][j]/1.
    ami_2v_polynom_multiplication(p_xy, 2, p_xx, 2, pol)

#---------------------------------------------------
#function to compute the partial derivatives of a 2 variable polynom. The degree of the derivative polynoms is assumed to be the same that the original one
def ami_2v_polynom_derivatives(p, N, p_x, p_y):   
    for i in range(N + 1):    
        for j in range(N + 1):
            p_x[i][j] = 0.0
            p_y[i][j] = 0.0
    for i in range(N + 1):    
        for j in range(N + 1):
            p_x[i - 1][j] = i * p[i][j]   
    for i in range(N + 1):    
        for j in range(N + 1):
            p_y[i][j - 1] = j * p[i][j] 

#function to multiply polinoms of 1 variable. the result is added to the output polynom COEFFICIENTs
def ami_1v_polynom_multiplication(p1, N1, p2, N2, p3):
    for i in range(N1 + 1):
        if (p1[i] != 0):
            for j in range(N2 + 1):
                if (p2[j] != 0):
                    p3[i + j] += p1[i] * p2[j]
    return p3

#---------------------------------------------------
#function to compute the determinant of a polynom matrix
def ami_polynom_determinant(p, Np, Nd, q): 
    p2 = [[[0.0 for i in range(19)] for j in range(6)] for m in range(6)]   
    q2 = [0.0 for i in range(Np + 1)]  
    if (Nd == 1):
        for i in range(19):
            q[i] = p[0][0][i]
        return
    cont = -1
    for i in range(Nd):
        for k in range(Np + 1):
            q2[k] = 0.0
        cont *= -1
        for k in range(Nd - 1):    
            for l in range(Nd - 1):
                for m in range(Np + 1):
                    if l >= i:
                        p2[k][l][m] = p[k + 1][l + 1][m]
                    else:
                        p2[k][l][m] = p[k + 1][l][m] 
        ami_polynom_determinant(p2, Np, Nd - 1, q2)    
        if (cont < 0):
            for m in range(Np + 1):
                q2[m] = -q2[m]               
        q = ami_1v_polynom_multiplication(p[0][i], Np, q2, Np, q) 

#---------------------------------------------------
#function to evaluate a polynomial and its derivate using Horner Algorithm
def ami_horner(pol, degree, x, fp):               
    PPX = pol[degree]
    PX = pol[degree]
    for k in range(degree - 1, 0, -1):
        PX = PX * x + pol[k]
        PPX = PPX * x + PX
    PX = PX * x + pol[0]
    fp = PPX
    return PX

#---------------------------------------------------
#function to compute the polynomial root in an interval. We use a combination of bisection and Newton-Raphson tecniques
def ami_root_bisection(pol, degree, a, b, TOL):
    maxiter = 100000    
    iterr = 0
    fp = 0
    a2 = a
    b2 = b

    #evaluate the polynomial at interval extrema
    fa = ami_horner(pol, degree, a, fp)
    fb = ami_horner(pol, degree, b, fp)
    #check if there is a root in the interval
    if (fa * fb > 0):
        return (0.5 * (a + b))
    #evaluate the interval middle point
    c = (a2 + b2) * 0.5
    c0 = c + 1e30
    #perform Newton-Raphson or Bisection iterations
    while((b2 - a2) > (TOL * math.fabs(b2) + 1e-10) and math.fabs(c - c0) > (TOL * math.fabs(c0) + 1e-10) and (iterr) < maxiter):
        iterr += 1
        fc = ami_horner(pol, degree, c, fp)
        #try Newton-Raphson
        c0 = c
        if (fp != 0.):
            c = c - fc / fp
            if(c >= a2 and c <= b2):
                continue
        #if Newton-Raphson fails, we try bisection
        c = (a2 + b2) * 0.5
        fc = ami_horner(pol, degree, c, fp)
        c0 = c + 1e30
        #check if we have arrive to the root
        if (math.fabs(fc) < TOL):
            break
        #perform bisection iteration
        if ((fa * fc) < 0):
            b2 = c
            fb = fc
        else:
            a2 = c
            fa = fc
        c = (a2 + b2) * 0.5
    if iterr >= maxiter:
        return ((a2 + b2) * 0.5)
    return c    

#---------------------------------------------------
#function to estimate real polynomial roots. We use a recursive procedure computing the derivative polynomial roots. This function only works for relative low polynomial degree
def ami_polynomial_root(pol, degree, root_r, root_i):
    TOL=1e-14
    N = degree
    fp = 0.0
    a = 0.0
    b = 0.0
    fa = 0.0
    fb = 0.0
    maxx = 0.0
    n_factor = 0.0
    Nr = 0

    p = [0.0 for i in range(N + 2)] 
    p_aux = [0.0 for i in range(N + 2)]
    ap = [0.0 for i in range(N + 1)] 
    f = [0.0 for i in range(N + 1)] 
    pol2 = [0.0 for i in range(N + 2)]  

    #POLYNOMIAL NORMALIZATION
    if (pol[0] != 1.):
        for k in range(N + 1):
            pol[k] /= pol[0]
    #REORDER THE POLYNOM
    for k in range(N + 1):
        pol2[k] = pol[N - k]
    pol2[N + 1] = 0.0   
    #NORMALIZE POLINOMIAL COEFFICIENTS TO AVOID POLYNOMIAL EVALUATION IN LARGE NUMBERS
    n_factor = 0.0     
    for k in range(N):
        if (math.fabs(pol2[k]) > 10.):
            maxx = math.log(math.fabs(pol2[k]) / 10.) / (N - k)
            if (maxx > n_factor):
                n_factor = maxx
    n_factor = math.exp(n_factor)   
    maxx = n_factor
    for k in range(N - 1, -1, - 1):
        pol2[k] /= maxx
        maxx *= n_factor
    #COMPUTE FACTORIAL NUMBERS
    f[0] = 1.0
    for k in range(1, N + 1):
        f[k] = f[k - 1] * k
    #COMPUTE THE INITIAL INTERVAL WHERE ALL ROOTS ARE INCLUDED
    maxx = math.fabs(pol2[0])
    for k in range(1, N + 1):
        if (math.fabs(pol2[k]) > maxx):
            maxx = math.fabs(pol2[k])  
    #INITIALIZE THE INTERVALS   
    p[0] = -1. - maxx / math.fabs(pol2[N])      
    #COMPUTE THE POLYNOMIAL DEGREE-1 DERIVATE ROOT
    p[1] = -(pol2[N - 1] * f[N - 1]) / (pol2[N] * f[N])
    for k in range(2, N + 2):
        p[k] =- p[0]
    for k in range(N + 1):
        p_aux[k] = p[k]   
    p_aux[N + 1] = 0.0    
    #COMPUTE THE DERIVATIVE POLINOMIAL ROOTS IN A RECURSIVE WAY 
    for k in range(N - 2, -1, - 1):
        #compute polynomial derivative coefficients
        for l in range(N - k + 1):
            ap[l] = pol2[l + k] * (f[k + l] / f[l])
        #check the valid intervals to compute derivative polynomial roots
        m = 1
        fa = ami_horner(ap, N - k, p_aux[0], fp) 
        for l in range(1, N - k + 1): 
            fb = ami_horner(ap, N-k, p_aux[l], fp)
            if ((fa * fb) <= 0 or math.fabs(fb) <= TOL):
                p[m] = p_aux[l]
                m += 1
                fa = fb  
        for l in range(m, N):
            p[l] = p[N]    
        #compute the derivative polynomial roots in each interval        
        for l in range(N - k + 1): 
            p_aux[l] = ami_root_bisection(ap, N-k, p[l-1], p[l], TOL)
    #STORE THE ROOTS
    root_i[0] = 0. #fit the complex component of the root to 0
    Nr = 0
    for k in range(1, N + 1):
        #check if the polynomial has a root in the point          
        a = (p_aux[k] + p_aux[k - 1]) * 0.5
        b = (p_aux[k] + p_aux[k + 1]) * 0.5
        fa = ami_horner(pol2, degree, a, fp)
        fb = ami_horner(pol2, degree, b, fp)
        if (fa * fb < 0 or math.fabs(ami_horner(pol2, degree, p_aux[k], fp)) < TOL):
            root_i[Nr] = 0 #we fit the complex component of the root to 0 
            root_r[Nr] = p_aux[k] * n_factor #we denormalize the root
            Nr += 1
    return Nr    

#---------------------------------------------------
#function to evaluate a 2 variable polynom in one of the variable value. The output is a 1 degree polynom
def ami_2v_polynom_to_1v_polynom(p1, N1, p3, z, flat):
    p = [0.0 for i in range(N1 + 1)] 
    if flat == 1:
        for i in range(N1 + 1):
            for j in range(N1 + 1):
                p[j] = p1[i][j]
            p3[i] = ami_polynomial_evaluation(p, N1, z)
    else:
        for i in range(N1 + 1):
            for j in range(N1 + 1):
                    p[j] = p1[j][i]
            p3[i]=ami_polynomial_evaluation(p, N1, z)

#---------------------------------------------------
#function to evaluate a 2 variable polynom in one point
def ami_2v_polynom_evaluation(p1, N1, x, y):
    p = [0.0 for i in range(N1 + 1)] 
    q = [0.0 for i in range(N1 + 1)] 
    for i in range(N1 + 1):
        for j in range(N1 + 1):
            p[j] = p1[i][j]
        q[i] = ami_polynomial_evaluation(p, N1, y)
    evall = ami_polynomial_evaluation(q, N1, x)
    return(evall)

#---------------------------------------------------
#function to update the lens distortion model by minimizing a 4 degree 2 variable polynom
def ami_lens_distortion_model_update_2v(coefficients, k1, k2, pol, max_radius):
    p_r3 = [[[0.0 for i in range(19)] for j in range(6)] for m in range(6)]
    x = [0.0 for i in range(3)] 
    pol_x = [[0.0 for i in range(5)] for j in range(5)]
    pol_y = [[0.0 for i in range(5)] for j in range(5)]
    p3 = [0.0 for i in range(5)] 
    b = [0.0 for i in range(app.max_grade_polinom + 1)] 
    # FILL THE AUXILIARY LENS DISTORTION MODEL
    for i in range(app.max_grade_polinom):
        b[i] = coefficients[i]
    #NORMALIZE POLYNOM COEFFICIENT
    sx = math.pow(pol[4][0], 0.25)
    sy = math.pow(pol[0][4], 0.25)
    if sx == 0 or sy == 0:
        FatalError('Processing Fatal Error')
    for i in range(5):    
        for j in range(5):
            if i > 0:
                pol[i][j] /= math.pow(sx, i)
            if j > 0:
                pol[i][j] /= math.pow(sy, j)
    #COMPUTE THE DERIVATIVES OF THE POLYNOM
    ami_2v_polynom_derivatives(pol, 4, pol_x, pol_y)  
    #FILL THE MATRIX TO COMPUTE THE DETERMINANT
    for i in range(4):    
        for m in range(5):  
            p_r3[2][i+2][m] = pol_x[3 - i][m]
            p_r3[1][i+1][m] = pol_x[3 - i][m]
            p_r3[0][i][m] = pol_x[3 - i][m]
            p_r3[5][i+2][m] = pol_y[3 - i][m]
            p_r3[4][i+1][m] = pol_y[3 - i][m]
            p_r3[3][i][m] = pol_y[3 - i][m]  
    #COMPUTE THE RESOLVENT POLYNOM
    pol_r = [0.0 for i in range(19)]               
    ami_polynom_determinant(p_r3, 18, 6, pol_r)
    #COMPUTE THE RESOLVENT POLYNOM DEGREE
    for i in range(19):   
        if (pol_r[i] != 0):
            Nr = i 
    #COMPUTE THE ROOT OF THE RESOLVENT POLYNOM
    rx = [0.0 for i in range(Nr)]                  
    ry = [0.0 for i in range(Nr)]                  
    b2 = [0.0 for i in range(Nr + 1)]    
    for i in range(Nr + 1): 
        b2[i] = pol_r[Nr - i]                
    Nr = ami_polynomial_root(b2, Nr, rx, ry) 
    #COMPUTE THE X COMPONENT BY REPLACING THE ROOTS IN THE DERIVATIVES OF THE POLYNOM   
    xr = 0.0
    yr = 0.0
    Emin = 10e90
    for i in range(Nr): 
        if (math.fabs(ry[i])> 0.000000000000001):
            continue
        ami_2v_polynom_to_1v_polynom(pol_x, 4, p3, rx[i], 1)
        M = ami_RootCubicPolynomial(p3, 3, x)
        for m in range(M):
            Energy = ami_2v_polynom_evaluation(pol, 4, x[m], rx[i])
            if (Energy < Emin):
                b[k1] = coefficients[k1] + (x[m] / sx)
                b[k2] = coefficients[k2] + (rx[i] / sy)
                if (test_compatibility_lens_distortion_model(b, app.max_grade_polinom, max_radius) == 0):
                    Emin = Energy
                    xr = rx[i]
                    yr = x[m]
        ami_2v_polynom_to_1v_polynom(pol_y, 4, p3, rx[i], 1)
        M = ami_RootCubicPolynomial(p3, 3, x)
        for m in range(M):
            Energy = ami_2v_polynom_evaluation(pol, 4, x[m], rx[i])
            if (Energy < Emin):
                b[k1] = coefficients[k1] + (x[m] / sx)
                b[k2] = coefficients[k2] + (rx[i] / sy)
                if (test_compatibility_lens_distortion_model(b, app.max_grade_polinom, max_radius) == 0):
                    Emin = Energy
                    xr = rx[i]
                    yr = x[m]
    #UPDATE THE DISTORSION POLYNOMIAL MODEL
    coefficients[k1] += (yr / sx)
    coefficients[k2] += (xr / sy)
    return 0

#---------------------------------------------------
#function to update the lens distortion polynomial model for 2 variables. If alpha>0, we adapt a[0] to minimize the square distence between distorted and undistorted points and we add a term to the polynomial also minimizing such distance with weight alpha
def ami_lens_distortion_estimation_2v(coefficients, grid_corners, x0, y0, k1, k2, alpha, max_radius):
    Error = 0.0
    pol_v2 = [[0.0 for i in range(5)] for j in range(5)]

    # UPDATE a[0] BY MINIMIZING THE DISTANCE OF THE DISTORTED POINTS TO THE UNDISTORTED POINTS
    if alpha > 0:
        sum_dd = 0.0
        sum_Ad = 0.0
        for m in range(app.lines):
            for i in range(app.points_in_line[m]):
                d = math.sqrt((grid_corners[m][i][0] - x0) * (grid_corners[m][i][0] - x0) + (grid_corners[m][i][1] - y0) * (grid_corners[m][i][1] - y0))
                A = 0.0
                for k in range(1, app.max_grade_polinom + 1):
                    A += coefficients[k] * math.pow(d, k + 1.0)
                sum_dd += d * d
                sum_Ad += A * d
        coefficients[0] = 1 - sum_Ad / sum_dd
    for m in range(app.lines):
        #WE UPDATE DE POLYNOM TO MINIMIZE
        ami_lens_distortion_polynomial_update_2v(coefficients, grid_corners, m, app.points_in_line[m], x0, y0, k1, k2, pol_v2)
        ami_lens_distortion_polynomial_update_distance_2v(coefficients, grid_corners, m, app.points_in_line[m], x0, y0, k1, k2, pol_v2, alpha)
    #UPDATE THE POLYNOMIAL LENS DISTORTION MODEL
    ami_lens_distortion_model_update_2v(coefficients, k1, k2, pol_v2, max_radius)  
    for i in range(app.lines):
        Error += ami_LensDistortionEnergyError(grid_corners[i], app.points_in_line[i], x0, y0, coefficients)

    return Error / app.lines 

#---------------------------------------------------
#function to compute the lens distortion model
def ami_lens_distortion_estimation(coefficients, grid_corners, x0, y0, k, alpha, max_radius):
    Error = 0.0
    pol = [0.0 for i in range(5)]

    #ADAPT a[0] TO MINIMIZE THE SQUARE OF THE DISTANCE BEWTEEN DISTORTED AND UNDISTORDED POINTS
    if alpha > 0:
        sum_dd = 0.0
        sum_Ad = 0.0 
        for m in range(app.lines):
            for i in range(app.points_in_line[m]):
                d = math.sqrt((grid_corners[m][i][0] - x0) * (grid_corners[m][i][0] - x0) + (grid_corners[m][i][1] - y0) * (grid_corners[m][i][1] - y0))
                A = 0
                for j in range(1, app.max_grade_polinom + 1):
                    A += coefficients[j] * math.pow(d, j + 1.0)
                sum_dd += d * d
                sum_Ad += A * d
        coefficients[0] = 1 - sum_Ad / sum_dd
    #COMPUTE THE LENS DISTORTION MODEL
    for i in range(app.lines):
        ami_lens_distortion_polynomial_update(coefficients, grid_corners, i, app.points_in_line[i], x0, y0, k, pol) 
    ami_lens_distortion_model_update(coefficients, grid_corners, k, pol, max_radius)     
    for i in range(app.lines):
        Error += ami_LensDistortionEnergyError(grid_corners[i], app.points_in_line[i], x0, y0, coefficients)
    return Error / app.lines    
 
#---------------------------------------------------
#function to calculate the best center of distortion
def search_for_best_center(coefficients, grid_corners, max_radius, width, height):
    x_c = math.floor(coefficients[app.max_grade_polinom + 1])   
    y_c = math.floor(coefficients[app.max_grade_polinom + 2])

    aux_image = [[-1.0 for i in range(height)] for j in range(width)]

    #RUN A SAFE PREVIOUS ITERATION TO AVOID CONVERGENCE PROBLEMS
    initial_Emin = ami_lens_distortion_estimation(coefficients, grid_corners, x_c, y_c, 2, 0., max_radius)
    initial_Emin = ami_lens_distortion_estimation(coefficients, grid_corners, x_c, y_c, 4, 0., max_radius)
    #RUN THE ALGEBRAIC METHOD FOR BOTH PARAMETERS IN ONE ITERATION
    initial_Emin = ami_lens_distortion_estimation_2v(coefficients, grid_corners, x_c, y_c, 2, 4, 0., max_radius)
  
    #initial_Emin has the algebraic solution for the given center of distortion.
    #The corresponding solution is in vector a
    #This initial_Emin solution is the initial one to compare with, for the patch search
    #scan iteratively the image looking for the best center of distortion
    last_Emin = 1e100
    best_x_c = x_c
    best_y_c = y_c
    patch_half = int(math.floor(20 / 2.0))
    n_iterations = 0
    while (math.fabs(initial_Emin - last_Emin) > 1.0e-6):
        last_Emin = initial_Emin
        n_iterations += 1
        lim_inf_x = x_c - patch_half
        lim_sup_x = x_c + patch_half
        lim_inf_y = y_c - patch_half
        lim_sup_y = y_c + patch_half
        #scan the rectangular patch: pixel precission, subpixel is reached through gradient
        for i in range(lim_inf_x, lim_sup_x + 1):
            if (i >= 0 and i < width): #check that we are inside the image (x-axis)
                for j in range(lim_inf_y, lim_sup_y + 1):
                    if (j >= 0 and j < height): #check that we are inside the image (y-axis)
                        if (aux_image[i][j] == -1.0): #check if this value is already calculated
                            coefficients[2] = 0.0
                            coefficients[4] = 0.0      #reset the a vector
                            #REMARK: ACTIVATE THE NEXT TWO LINES IF THERE ARE CONVERGENCE PROBLEMS
                            Emin = ami_lens_distortion_estimation(coefficients, grid_corners, i, j, 2, 0.,max_radius)
                            Emin = ami_lens_distortion_estimation(coefficients, grid_corners, i, j, 4, 0.,max_radius)
                            Emin = ami_lens_distortion_estimation_2v(coefficients, grid_corners, i, j, 2, 4, 0.,max_radius)
                            aux_image[i][j] = Emin #save the value for possible use in new iterations
                        else:
                            Emin = aux_image[i][j]
                        if (Emin < initial_Emin):
                            initial_Emin = Emin #We save the best solution 
                            best_x_c = i
                            best_y_c = j #save the best center of distortion found 
        x_c = best_x_c
        y_c = best_y_c
        if (n_iterations == 20):
            break
        patch_half *= 2    
    coefficients[app.max_grade_polinom + 1] = best_x_c   
    coefficients[app.max_grade_polinom + 2] = best_y_c

    return 

#---------------------------------------------------
#function to estimate the position of 2D points (pixels) for the actual
def calculate_points(amin, Points2D, x0, y0, points_in_line):
    for c in range(points_in_line):
        d1 = math.sqrt(math.pow(Points2D[c][0] - x0, 2.0) + math.pow(Points2D[c][1] - y0, 2.0))
        sol = amin[app.max_grade_polinom]
        for j in range(app.max_grade_polinom - 1, -1, -1):
             sol = sol * d1 + amin[j]
        Points2D[c][0] = (Points2D[c][0] - x0) * sol + x0     
        Points2D[c][1] = (Points2D[c][1] - y0) * sol + y0
    #print(Points2D)    

#---------------------------------------------------
#function to compute a line equation by minimizing the distance to a point collection
def ami_line2d_calculation(line, Points2D, points_in_line):
    zero = 0.00000000000001
    r = [[0.0 for i in range(3)] for j in range(4)]
    suu = 0.0
    suv = 0.0
    svv = 0.0
    um = 0.0
    vm = 0.0 

    for i in range(points_in_line):
        um = um + Points2D[i][0]   
        vm = vm + Points2D[i][1]

    um = um / points_in_line
    vm = vm / points_in_line

    for i in range(points_in_line):
        suu = suu + (Points2D[i][0] - um) * (Points2D[i][0] - um)
        svv = svv + (Points2D[i][1] - vm) * (Points2D[i][1] - vm)
        suv = suv + (Points2D[i][0] - um) * (Points2D[i][1] - vm)
    suu = suu / points_in_line    
    svv = svv / points_in_line
    suv = suv / points_in_line

    if math.fabs(suv) <= zero:
        if suu < svv and svv > zero:
            line[0] = 1.0
            line[1] = 0.0
            line[2] = -um
            return
        if svv < suu and suu > zero:
            line[0] = 0.0
            line[1] = 1.0
            line[2] = -vm
            return    
    r[2][1] = 1.0        
    r[3][1] = 1.0
    r[0][0] = 1.0
    r[1][0] = 1.0
    h = 0.5 * (suu - svv) / suv

    if h > 0:
        r[0][1] = -h - math.sqrt(1. + h * h)
        r[0][2] = -(um + r[0][1] * vm)
        r[1][1] = -1. / r[0][1]
        r[1][2] = -(um + r[1][1] * vm)
        r[2][0] = h + math.sqrt(1. + h * h)
        r[2][2] = -(r[2][0] * um + vm)
        r[3][0] = -1. / r[2][0]
        r[3][2] = -(r[3][0] * um + vm)
    else:
        r[0][1] = -h + math.sqrt(1 + h * h)
        r[0][2] = -(um + r[0][1] * vm)
        r[1][1] = -1. / r[0][1]
        r[1][2] = -(um + r[1][1] * vm)
        r[2][0] = h - math.sqrt(1 + h * h)
        r[2][2] = -(r[2][0] * um + vm)
        r[3][0] = -1. / r[2][0]
        r[3][2] = -(r[3][0] * um + vm)

    for j in range(4):
        norm = math.sqrt(r[j][0] * r[j][0] + r[j][1] * r[j][1])
        for i in range(3):
            r[j][i] = r[j][i] / norm

    minn = 0.0
    k = 0        

    for c in range(points_in_line):
        aux = r[0][0] * Points2D[c][0] + r[0][1] * Points2D[c][1] + r[0][2]
        minn = minn + math.fabs(aux)

    for j in range(1, 4):
        h = 0.0    
        for c in range(points_in_line):
            aux = r[j][0] * Points2D[c][0] + r[j][1] * Points2D[c][1] + r[j][2]
            h = h + math.fabs(aux)
        if h < minn:
            k = j
            minn = h
    line[0] = r[k][0]
    line[1] = r[k][1]
    line[2] = r[k][2] 

#---------------------------------------------------
#function to minimize in one dimension (searching lambda)
def find_lambda(lambda1, lambda2, lambda3, f_1, f_2, f_3, amin_copy, amin, grid_corners, grad_f, change_k):
    Naa = 7
    f_tmp1 = f_3
    f_tmp2 = f_2
    f = f_1
    lambda_ = 0.0

    lambda_tmp1 = (lambda2 + lambda1) / 2.0
    for i in range(Naa):
        if change_k[i] == 1:
            amin_copy[i] = amin[i] - lambda_tmp1 * grad_f[i]
    f_tmp1 = distance_function(amin_copy, grid_corners)
    if f_tmp1 < f_1:
        return lambda_tmp1

    lambda_tmp2 = (lambda3 + lambda2) / 2.0           
    for i in range(Naa):
        if change_k[i] == 1:
            amin_copy[i] = amin[i] - lambda_tmp2 * grad_f[i]
    f_tmp2=distance_function(amin_copy, grid_corners)
    if f_tmp2 < f_1:
        return lambda_tmp2
    f = f_1 
    while 1:
        f_min = f
        lambda_ = lambda1 + 1e-8 
        for i in range(Naa):
            if change_k[i] == 1:
                amin_copy[i] = amin[i] - lambda_ * grad_f[i]          
        f = distance_function(amin_copy, grid_corners)
        if f <= f_min:
            break
    return lambda_        

#---------------------------------------------------
#function to build and minimize the cuadratic TPP polynom
def minimize_cuadratic_polynom(lambda1, lambda2, lambda3, f_1, f_2, f_3, amin_copy, amin, grid_corners, grad_f, change_k):
    a12 = lambda1 - lambda2
    a23 = lambda2 - lambda3
    a31 = lambda3 - lambda1
    b12 = math.pow(lambda1, 2.0) - math.pow(lambda2, 2.0)
    b23 = math.pow(lambda2, 2.0) - math.pow(lambda3, 2.0)
    b31 = math.pow(lambda3, 2.0) - math.pow(lambda1,2.0)
 
    min_lambda = 0.5 * (b23 * f_1 + b31 * f_2 + b12 * f_3)
    if (a23 * f_1 + a31 * f_2 + a12 * f_3) == 0:
        FatalError('Processing Fatal Error')
    min_lambda = min_lambda / (a23 * f_1 + a31 * f_2 + a12 * f_3)    

    if ((lambda1 < min_lambda) and (min_lambda < lambda3)):
        return min_lambda
    else:
        min_lambda = find_lambda(lambda1, lambda2, lambda3, f_1, f_2, f_3, amin_copy, amin, grid_corners, grad_f, change_k)    
        return min_lambda

#---------------------------------------------------
#function to find the minimum of the interpolating polynom
def cuadratic_fitting(amin_copy, coefficients, grid_corners, lambda1, lambda2, lambda3, f_1, f_2, f_3, grad_f, change_k):
    minimo_lambda = 0.0
    f_minimo = 0.0
    error = 1e100           
    iterations_lambda = 0 
    tol_lambda = 1e-8     
    Naa = 7
    #We loop till getting the minimum
    while error > tol_lambda:
        minimo_lambda = minimize_cuadratic_polynom(lambda1, lambda2, lambda3, f_1, f_2, f_3, amin_copy, coefficients, grid_corners, grad_f, change_k)                    
        for i in range(Naa):
            if change_k[i] == 1:
                amin_copy[i] = coefficients[i] - minimo_lambda * grad_f[i]
        f_minimo = distance_function(amin_copy, grid_corners)
        if minimo_lambda > lambda2:     
            if f_minimo > f_2:
                lambda3 = minimo_lambda
                f_3 = f_minimo
            else:
                lambda1 = lambda2
                f_1 = f_2
                lambda2 = minimo_lambda
                f_2 = f_minimo 
        else:
            if f_minimo >= f_2:
                lambda1 = minimo_lambda
                f_1 = f_minimo
            else:
                lambda3 = lambda2
                f_3 = f_2
                lambda2 = minimo_lambda
                f_2 = f_minimo 
        error = math.fabs(lambda3 - lambda1)
        if f_minimo == f_2:
            lambda2 = lambda2 + tol_lambda        
        iterations_lambda = iterations_lambda + 1
        if iterations_lambda == 10:
            return lambda2    
    return lambda2

#---------------------------------------------------
#function Unidimensional lambda miminization
def minimize_lambda(coefficients, grid_corners, grad_f, f_objective, change_k):
    Naa = 7   
    amin_copy = [0.0 for i in range(7)]             
    f_1 = f_objective
    tol_ff = 1.0e-10
    f_3 = 0.0
    #search the TTP points
    lambda1 = 0.0
    #search lambda2
    lambda2 = math.fabs(grad_f[2])
    for i in range(Naa):
        amin_copy[i] = coefficients[i] - lambda2 * grad_f[i]
    f_2 = distance_function(amin_copy, grid_corners)

    if f_2 > f_1:
        lambda3 = lambda2
        f_3 = f_2
        #search lambda2 by dividing the (lambda1,lambda3) interval
        lambda2 = lambda3 / 2.0
        while 1:
            last_f2 = f_2
            for i in range(Naa):
                if change_k[i] == 1:
                    amin_copy[i] = coefficients[i] - lambda2 * grad_f[i]
                f_2 = distance_function(amin_copy, grid_corners)   
                if f_2 < f_1:
                    break
                if last_f2 != 0:
                    if math.fabs((f_2 - last_f2) / last_f2) <= tol_ff: 
                        return lambda2
                lambda2 = lambda2 / 2.0                    
    else:
        #search lambda3 by increasing the (lambda1,lambda2) interval
        lambda3 = lambda2 * 2.0 
        while 1:
            last_f3 = f_3
            for i in range(Naa):               
                if change_k[i] == 1:
                    amin_copy[i] = coefficients[i] - lambda3 * grad_f[i]
            f_3 = distance_function(amin_copy, grid_corners) 
            if f_3 > f_2:
                break
            if last_f3 != 0:    
                if math.fabs((f_3 - last_f3) / last_f3) <= tol_ff: 
                    return lambda3
            lambda3 = 2 * lambda3
    #We have the points satisfying the TTP condition lambda1,f_1lambda_2,f_2lambda3,f_3 minimize the cuadratic polynom
    lambda_ = cuadratic_fitting(amin_copy, coefficients, grid_corners, lambda1, lambda2, lambda3, f_1, f_2, f_3, grad_f, change_k)
    return lambda_

#---------------------------------------------------
#function to apply the zoom strategy
def ami_lens_distortion_zoom_normalization(coefficients, grid_corners):
    N = 0
    x0 = coefficients[app.max_grade_polinom + 1]   
    y0 = coefficients[app.max_grade_polinom + 2]  
    #UPDATE a BY ESTIMATING A ZOOM FACTOR Z
    sum_Ad = 0.0
    for r in range(app.lines):
        for c in range(app.points_in_line[r]):     
            N = N + 1
            d = math.sqrt((grid_corners[r][c][0] - x0) * (grid_corners[r][c][0] - x0) + (grid_corners[r][c][1] - y0) * (grid_corners[r][c][1] - y0))
            A = coefficients[0]
            for k in range(1, app.max_grade_polinom + 1):
                A = A + coefficients[k] * math.pow(d, k)
            sum_Ad = sum_Ad + A * A
    Z = math.sqrt(N / sum_Ad)        
    for k in range(0, app.max_grade_polinom + 1):
        coefficients[k] = coefficients[k] * Z

#---------------------------------------------------
#function to be optimized by the gradient (objective distance function)
def distance_function(coefficients, grid_corners):
    line = [0.0 for i in range(3)] 
    amin = [0.0 for i in range(app.max_grade_polinom + 1)] 
    for i in range(app.max_grade_polinom + 1):
        amin[i] = coefficients[i]
    x0 = coefficients[app.max_grade_polinom + 1]   
    y0 = coefficients[app.max_grade_polinom + 2]  

    f_objective = 0.0 

    points_2D_modified = [[0.0, 0.0] for i in range(app.max_grid_columns)]
    for r in range(app.lines):
        for c in range(app.points_in_line[r]):
            points_2D_modified[c][0] = grid_corners[r][c][0]
            points_2D_modified[c][1] = grid_corners[r][c][1]
        #print(points_2D_modified)    
        calculate_points(amin, points_2D_modified, x0, y0, app.points_in_line[r])
        ami_line2d_calculation(line, points_2D_modified, app.points_in_line[r])  
        a = line[1]
        b = line[0]
        c = line[2]
        #print("line", line)
        tmp = math.pow(a, 2.0) + math.pow(b, 2.0)
        sum_ = 0.0
        for j in range(app.points_in_line[r]):
            sum_ = sum_ + math.pow(b * points_2D_modified[j][0] + a * points_2D_modified[j][1] + c, 2.0)
        sum_ = sum_ / tmp
        sum_ = sum_ / app.points_in_line[r]
        f_objective = f_objective + sum_ 
    f_objective = f_objective / app.lines
    return f_objective

#---------------------------------------------------
#function to minimize the distance function by gradient
def gradient_method(coefficients, grid_corners, change_k, zoom):
    grad_f = [0.0 for i in range(7)]  
    n_iterations = 0
    last_f = 1.0e100
    kk = 0.0
    f_objective = distance_function(coefficients, grid_corners) 
    while math.fabs((f_objective - last_f) / last_f) > 1.0e-6:
        last_f = f_objective
        for i in range(7):
            #Move along each axis and incremental step to calculate the derivative
            if change_k[i] == 1:
                kk = coefficients[i]
                coefficients[i] = kk + 1.0e-10
                grad_f[i] = (distance_function(coefficients, grid_corners) - f_objective) / 1.0e-10
                #print("grad_f[%d] = %f" % (i, grad_f[i]))
                coefficients[i] = kk 
        #Activate to stop the gradient when the gradient_norm<tol_norma gradient_norm=0.0
        lambda_ = minimize_lambda(coefficients, grid_corners, grad_f, f_objective, change_k)
        #print("lambda_", lambda_)
        for i in range(7):
            if change_k[i] == 1:
                coefficients[i] = coefficients[i] - lambda_ * grad_f[i]
                #print("coefficients[%d] = %f" % (i, coefficients[i]))
        n_iterations = n_iterations + 1
        f_objective=distance_function(coefficients, grid_corners)   
        #Activate to have the trace of the execution
        if n_iterations == 100:
            break 
    if zoom == 1:
        ami_lens_distortion_zoom_normalization(coefficients, grid_corners)           
    
#---------------------------------------------------
#function to calculate the factor_n, needed for transforming k
def calculate_factor_n(grid_corners, x0, y0):
    factor_n = 0.0
    cont = 0
    for r in range(app.lines):
        for c in range(app.points_in_line[r]):
            factor_n = factor_n + (grid_corners[r][c][0] - x0) * (grid_corners[r][c][0] - x0) + (grid_corners[r][c][1] - y0) * (grid_corners[r][c][1] - y0)
            cont += 1
    factor_n = math.sqrt(factor_n / cont)   
    return factor_n

#---------------------------------------------------
#function to evaluate a polynom using the Horner algorithm
def ami_polynomial_evaluation(a, Na, x):
    sol = a[Na]
    for i in range(Na - 1, -1, -1):
        sol = sol * x + a[i]
    return sol    

#---------------------------------------------------
#function to compute the lens distortion energy error (the residual variance of the point distribution)
def ami_LensDistortionEnergyError(grid_corners, Np, x0, y0, amin):
    x2 = [0.0 for i in range(Np)] 
    y2 = [0.0 for i in range(Np)] 

    #COMPUTE THE POINT TRANSFORMATION USING THE LENS DISTORTION MODEL
    for i in range(Np):
        d = math.sqrt((grid_corners[i][0] - x0) * (grid_corners[i][0] - x0) + (grid_corners[i][1] - y0) * (grid_corners[i][1] - y0))
        A = ami_polynomial_evaluation(amin, app.max_grade_polinom, d)
        x2[i] = x0 + (grid_corners[i][0] - x0) * A
        y2[i] = y0 + (grid_corners[i][1] - y0) * A
    #COMPUTE THE VARIANCE OF THE TRANSFORMED POINTS
    x2_m = 0.0
    for i in range(Np):
        x2_m = x2_m + x2[i]
    x2_m = x2_m / Np
    s_xx = 0.0
    for i in range(Np):
        s_xx = s_xx + (x2[i] - x2_m) * (x2[i] - x2_m)
    s_xx= s_xx / Np
    y2_m = 0.0
    for i in range(Np):
        y2_m = y2_m + y2[i]
    y2_m = y2_m / Np
    s_yy = 0.0
    for i in range(Np):
        s_yy = s_yy + (y2[i] - y2_m) * (y2[i] - y2_m)
    s_yy = s_yy / Np
    s_xy = 0.0
    for i in range(Np):
        s_xy = s_xy + (y2[i] - y2_m) * (x2[i] - x2_m)
    s_xy = s_xy / Np

    return s_xx * s_yy + s_xy * s_xy  

#---------------------------------------------------
#function to compute the lens distortion vmin energy error of the point distribution
def ami_LensDistortionEnergyError_Vmin(grid_corners, Np, x0, y0, amin): 
    x2 = [0.0 for i in range(Np)] 
    y2 = [0.0 for i in range(Np)]   

    #COMPUTE THE POINT TRANSFORMATION USING THE LENS DISTORTION MODEL
    for i in range(Np):
        d = math.sqrt((grid_corners[i][0] - x0) * (grid_corners[i][0] - x0) + (grid_corners[i][1] - y0) * (grid_corners[i][1] - y0))
        A = ami_polynomial_evaluation(amin, app.max_grade_polinom, d)
        x2[i] = x0 + (grid_corners[i][0] - x0) * A
        y2[i] = y0 + (grid_corners[i][1] - y0) * A  
    #COMPUTE THE VARIANCE OF THE TRANSFORMED POINTS
    x2_m = 0.0
    for i in range(Np):
        x2_m = x2_m + x2[i]
    x2_m = x2_m / Np
    s_xx = 0.0
    for i in range(Np):
        s_xx = s_xx + (x2[i] - x2_m) * (x2[i] - x2_m)
    s_xx= s_xx / Np
    y2_m = 0.0
    for i in range(Np):
        y2_m = y2_m + y2[i]
    y2_m = y2_m / Np
    s_yy = 0.0
    for i in range(Np):
        s_yy = s_yy + (y2[i] - y2_m) * (y2[i] - y2_m)
    s_yy = s_yy / Np
    s_xy = 0.0
    for i in range(Np):
        s_xy = s_xy + (y2[i] - y2_m) * (x2[i] - x2_m)
    s_xy = s_xy / Np
    if s_xx > s_yy:
        s_max = s_xx
    else:        
        s_max = s_yy
    if s_max == 0:
        FatalError('Processing Fatal Error')        
    return ((s_xx * s_yy + s_xy * s_xy) / s_max) 

#---------------------------------------------------
#function to execute the gradient method
def optimize(coefficients, grid_corners, grid_corners_, factor_n, zoom):
    x00 = coefficients[app.max_grade_polinom + 1]
    y00 = coefficients[app.max_grade_polinom + 2]
    n_iterations=0

    #TO CALCULATE IN NORMALIZED COORDINATES
    coefficients[app.max_grade_polinom + 1] = 0
    coefficients[app.max_grade_polinom + 2] = 0

    change_k = [0 for i in range(7)]
    for i in range(app.max_grade_polinom + 1):
        change_k[i] = 1
   
    if app.config.optimize_center == True:
        #OPTIMIZE THE CENTER OF DISTORTION
        change_k[5] = 1
        change_k[6] = 1
    gradient_method(coefficients, grid_corners, change_k, zoom)
    
    #RECOVER THE CENTER OF DISTORTION, OPTIMIZED OR NOT
    x0 = coefficients[app.max_grade_polinom + 1]
    y0 = coefficients[app.max_grade_polinom + 2]
    #print(coefficients)

    #CALCULATE THE NEW CENTER OF DISTORTION IN PIXELS
    if app.config.optimize_center == True:
        x0 = x0 + x00
        y0 = y0 + y00
    else:
        x0 = x00    
        y0 = y00
    #TRANSFORM THE SOLUTION IN NORMALIZED COORDINATES TO PIXEL UNITS
    if app.config.optimize_center == False:
        paso = 1.0
        for i in range(app.max_grade_polinom + 1):
            coefficients[i] = coefficients[i] * paso
            paso = paso / factor_n
    else:
        factor_n_new = calculate_factor_n(grid_corners_, x0, y0)    
        paso = 1.0
        for i in range(app.max_grade_polinom + 1):
            #print("paso", paso)
            coefficients[i] = coefficients[i] * paso
            paso = paso / factor_n_new
    #print("f", factor_n, factor_n_new)        
    #UPDATE THE SOLUTION (CENTER OF DISTORTION)        
    coefficients[app.max_grade_polinom + 1] = x0        
    coefficients[app.max_grade_polinom + 2] = y0

    amin = [0.0 for i in range(app.max_grade_polinom + 1)]
    for i in range(app.max_grade_polinom + 1):
        amin[i] = coefficients[i]

    Emin=0.0
    #print(amin)
    for r in range(app.lines):   
        Emin = Emin + ami_LensDistortionEnergyError(grid_corners_[r], app.points_in_line[r], x0, y0, amin) 
        #print(Emin)
    Emin = Emin / app.lines   
    Vmin = 0.0
    for i in range(app.lines): 
        Vmin = Vmin + ami_LensDistortionEnergyError_Vmin(grid_corners_[r], app.points_in_line[r], x0, y0, amin)
        #print(Vmin)   
    Vmin = Vmin / app.lines   
    D = distance_function(coefficients, grid_corners_) 

    #for i in range(app.max_grade_polinom + 1):
        #print(i, amin[i]) 

#---------------------------------------------------
#function to calculate the solution for the gradient method applied from the algebraic method solution
def algebraic_method_pre_gradient(coefficients, grid_corners, grid_corners_, factor_n, zoom, max_radius):
    x00 = coefficients[app.max_grade_polinom + 1]
    y00 = coefficients[app.max_grade_polinom + 2]
    if app.config.optimize_center == False:
        coefficients[0] = 1.0
        for i in range(1, app.max_grade_polinom + 1):
            coefficients[i] = 0.0
        #RUN A SAFE PREVIOUS ITERATION TO AVOID CONVERGENCE PROBLEMS
        ami_lens_distortion_estimation(coefficients, grid_corners, 0., 0., 2, 0.,max_radius)
        ami_lens_distortion_estimation(coefficients, grid_corners, 0., 0., 4, 0.,max_radius)
        #RUN THE ALGEBRAIC METHOD FOR BOTH PARAMETERS IN ONE ITERATION
        ami_lens_distortion_estimation_2v(coefficients, grid_corners, 0., 0., 2, 4, 0.,max_radius)
        #SET THE ORIGINAL CENTER OF DISTORTION
        coefficients[app.max_grade_polinom + 1] = x00
        coefficients[app.max_grade_polinom + 2] = y00
        # APPLY THE GRADIENT METHOD FROM THE SOLUTION 
        optimize(coefficients, grid_corners, grid_corners_, factor_n, zoom)    
    else:
        paso=1.0
        for i in range(app.max_grade_polinom + 1):
            coefficients[i] = coefficients[i] * paso
            paso = paso / factor_n
        coefficients[app.max_grade_polinom + 1] = x00
        coefficients[app.max_grade_polinom + 2] = y00   
        # APPLY THE GRADIENT METHOD FROM THE SOLUTION 
        optimize(coefficients, grid_corners, grid_corners_, factor_n, zoom)    
    return

#---------------------------------------------------
#Build an intermediate vector with values of L-1(r) for d = (0 to max_distance_corner)
def build_l1r_vector(l1r, max_distance_corner, max_grade_polinom, coefficients):
    root = 1.
    b = [0.0 for i in range(max_grade_polinom + 2)]
    b2 = [0.0 for i in range(max_grade_polinom + 2)]

    #DEFINE THE POLYNOM WE NEED TO COMPUTE ROOTS AND THE DERIVATIVE
    for i in range (1, max_grade_polinom + 2):
        b[i] = coefficients[i - 1]
    for i in range (0, max_grade_polinom + 1):
        b2[i] = coefficients[i] * (i + 1)
    #BUILD THE VECTOR OF INVERSE LENS DISTORTION FUNCTION        
    for dist in range (1, int(max_distance_corner + 1.5)):
        #DEFINE THE POLYNOM WE NEED TO COMPUTE ROOTS AND THE DERIVATIVE
        b[0] = -dist

        #NEWTON-RAPHSON TO COMPUTE THE POLYNOMIAL ROOT
        for k in range(10000):
            pol_eval = ami_polynomial_evaluation(b, max_grade_polinom + 1, root)
            pol_der = ami_polynomial_evaluation(b2, max_grade_polinom, root)
            if pol_der == 0.0:
                break
            root2 = root - pol_eval / pol_der  
            if (math.fabs(root - root2) < (math.fabs(root) * 1e-8)):
                root = root2
                break
            root = root2   
        l1r[dist] = root / dist 
    l1r[0] = l1r[1]

#---------------------------------------------------
#ESTIMATE AN UNDISTORTED IMAGE USING a DISTORTION MODEL (inverse method, accelerated)
def undistort_image_inverse_fast(input_image, output_image, width0, height0, coefficients, cx, cy, image_amplification_factor):
    size = width0 * height0
    width = width0 * image_amplification_factor
    height = height0 * image_amplification_factor
    sizenew = int(width * height)

    dc = CPoint(cx, cy)
    corner = CPoint()
    max_distance_corner = CPoint()
    d1 = max_distance_corner.Norm(corner)
    corner.y = height0
    distance_corner = CPoint()
    d2 = distance_corner.Norm(corner)
    if d2 > d1:
        d1 = d2
    corner.x = width0    
    d2 = distance_corner.Norm(corner)
    if d2 > d1:
        d1 = d2
    corner.y = 0.0
    d2 = distance_corner.Norm(corner)
    if d2 > d1:
        d1 = d2

    #BUILD INTERMEDIATE VECTOR    
    l1r = [0 for i in range(int(d1 + 1.5))]
    if build_l1r_vector(l1r, d1, app.max_grade_polinom, coefficients) == -1:
        return
    for nc in range(3):
        n2 = nc * size
        n2new = nc * sizenew
        for i in range(int(height)):
            for j in range(int(width)):
                temp = CPoint(j / image_amplification_factor, i / image_amplification_factor)
                d = CPoint(dc.x, dc.y)
                distance_center = d.Norm(temp)
                #INTERPOLATION
                ind = int(distance_center)
                dl1r = l1r[ind] + (distance_center - ind) * (l1r[ind + 1] - l1r[ind])
                p = CPoint()
                p.x = dc.x + (temp.x - dc.x) * dl1r
                p.y = dc.y + (temp.y - dc.y) * dl1r
                m = int(p.y)
                n = int(p.x)

                if (0 <= m and m < height0 and 0 <= n and n < width0):
                    #COLOUR INTERPOLATION
                    di = p.y - m
                    dj = p.x - n
                    k = int(i * width + j)
                    k0 = m * width0 + n
                    accum = 0.0
                    w_accum = 0.0
                    w = ((1. - di) * (1. - dj))

                    accum += w * input_image[k0 + n2]
                    w_accum += w

                    if ((di * (1. - dj)) > 0. and (m + 1) < height0):
                        k0 = (m + 1) * width0 + n
                        w = (di * (1. - dj))
                        accum += w * input_image[k0 + n2]
                        w_accum += w
                    if (((1 - di) * (dj)) > 0. and (n + 1) < width0):
                        k0 = (m) * width0 + n + 1
                        w = (1. - di) * (dj)
                        accum += w * input_image[k0 + n2]
                        w_accum += w
                    if (((di) * (dj)) > 0. and (n + 1) < width0 and (m + 1) < height0):
                        k0 = (m + 1) * width0 + n + 1
                        w = (di) * (dj)
                        accum += w * input_image[k0 + n2]
                        w_accum += w
                  

                    if w_accum > 0.:
                        output_image[k + n2new] = int(accum / w_accum)

#---------------------------------------------------
def ComputeSolution():
    if app.grid_corners_defined == False:
        return

    width = app.InputImageTk.width
    height = app.InputImageTk.height

    #TO AVOID A BAD SOLUTION WHEN OPTIMIZING THE CENTER OF DISTORTION USING THE A MINIMUM NUMBER OF PRIMITIVES    
    if app.lines == 1:
        app.config.optimize_center = False
    
    coefficients = [0.0 for i in range(app.max_grade_polinom + 8)]
    coefficients[0] = 1.0
    coefficients[app.max_grade_polinom + 1] = app.center_x
    coefficients[app.max_grade_polinom + 2] = app.center_y

    for r in range(app.lines):
        for c in range(app.points_in_line[r]):
            app.grid_corners[r][c][0] = app.grid_corners_[r][c][0]
            app.grid_corners[r][c][1] = app.grid_corners_[r][c][1]
            #print("x[%d][%d] = %f; y[%d][%d] = %f;" % (r, c, app.grid_corners[r][c][0], r, c, app.grid_corners[r][c][1]))                 

    #SELECT THE BIGGEST VALUE (FOR THE INPUT IMAGE): THE POINT IN A CORNER, TO BE USED IN THE LOCATION OF THE ALGEBRAIC ROOTS
    xtmp = width - width / 2.0   
    ytmp = height - height / 2.0
    max_radius = math.sqrt(math.pow(xtmp, 2.0) + math.pow(ytmp, 2.0))

    #OPTIMIZED CENTER: WE SCAN A PATCH CENTERED AT THE GIVEN CENTER OF DISTORTION:
    #SELECT THE BEST SOLUTION AND THEN APPLIED THE GRADIENT METHOD,
    #TO GET THE SOLUTIONS USING ZOOM AND NOT USING ZOOM
    if app.config.optimize_center:
        search_for_best_center(coefficients, app.grid_corners, max_radius, width, height)
        app.center_x = coefficients[app.max_grade_polinom + 1]
        app.center_y = coefficients[app.max_grade_polinom + 2]
  
    #ALGEBRAIC METHOD & GRADIENT METHOD WORK WITH NORMALIZED UNITS
    factor_n = 0.0
    cont = 0

    #print(app.center_x, app.center_y)
    
    for r in range(app.lines):
        for c in range(app.points_in_line[r]):
            #print(r, c, app.grid_corners[r][c][0], app.grid_corners[r][c][1])
            app.grid_corners[r][c][0] -= app.center_x
            app.grid_corners[r][c][1] -= app.center_y
            cont += 1
            factor_n += (app.grid_corners[r][c][0] * app.grid_corners[r][c][0]) + (app.grid_corners[r][c][1] * app.grid_corners[r][c][1])
            #print(grid_corners[r][c][0], grid_corners[r][c][1], factor_n)
    if factor_n == 0:
        FatalError('Processing Fatal Error')        
    factor_n = math.sqrt(factor_n / cont)       
    for r in range(app.lines):
        for c in range(app.points_in_line[r]):
             app.grid_corners[r][c][0] = app.grid_corners[r][c][0] / factor_n
             app.grid_corners[r][c][1] = app.grid_corners[r][c][1] / factor_n

    xtmp = xtmp / factor_n   
    ytmp = ytmp / factor_n   
    max_radius = math.sqrt(math.pow(xtmp, 2.0) + math.pow(ytmp, 2.0))    
    #print(factor_n, max_radius)

    #ALGEBRAIC METHOD FROM TRIVIAL SOLUTION + GRADIENT TO IMPROVE THE ALGEBRAIC SOLUTION; WITHOUT ZOOM
    #ALGEBRAIC METHOD & GRADIENT METHOD WORK WITH NORMALIZED UNITS
    #NO ZOOM APPLIED
    algebraic_method_pre_gradient(coefficients, app.grid_corners, app.grid_corners_, factor_n, 0, max_radius)
    #ALGEBRAIC METHOD FROM TRIVIAL SOLUTION + GRADIENT TO IMPROVE THE ALGEBRAIC SOLUTION; WITH ZOOM
    algebraic_method_pre_gradient(coefficients, app.grid_corners, app.grid_corners_, factor_n, 1, max_radius)

    undistort_image_inverse_fast(app.DistoredImage, app.UndistoredImage, app.InputImageTk.width, app.InputImageTk.height, coefficients, app.center_x, app.center_y , 1.0)
    app.UndistoredImageChanged = True
    line = app.lines
    app.lines = 0
    ResizeImage()
    DrawObjects()
    app.label["text"] = 'Calculated successfully \n Undistorted image \n Left button double click to exit' 

    for i in range(1, app.max_grade_polinom + 1):
        coefficients[i] *= math.pow(app.resize_coefficient, i)
    
    #write results
    dictionary = {
        "distortion": {
            "type": "radial_simple",
            "centerX": app.center_x / app.resize_coefficient,
            "centerY": app.center_y / app.resize_coefficient,
            "coefficients": [ coefficients[i] for i in range(app.max_grade_polinom + 1) ]
        }
    }
    with open(app.config.output_file, "w") as outfile:
        json.dump(dictionary, outfile)

    print(dictionary)    

    #with open(app.config.output_file, 'r') as openfile:
        # Reading from json file
    #    json_object = json.load(openfile)
     
    #print(json_object)

    app.calculation_done = True

    #Uncomment next row if OpenCV coefficients not needed 
    app.grid_was_reset = True
    if app.grid_was_reset == True:   
        return

    #OpenCV test of undistorion    
    nY = app.points_in_line[0]    
    nX = line
    square_size = 0.023
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)
    object_points = []
    image_points = []
    object_points_3D = np.zeros((nX * nY, 3), np.float32)       
    # These are the x and y coordinates                                              
    object_points_3D[:,:2] = np.mgrid[0:nY, 0:nX].T.reshape(-1, 2) 
    object_points_3D = object_points_3D * square_size

    object_points.append(object_points_3D)

    corners_2 = cv2.cornerSubPix(app.ImageOpenCVgray, app.OpenCVfoundCorners, (11,11), (-1,-1), criteria) 

    image_points.append(corners_2)

    #cv2.drawChessboardCorners(gray, GRID, app.OpenCVfoundCorners, app.grid_corners_defined)
    #cv2.imshow('img_captured_corners', gray) 

    ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(object_points, 
                                                image_points, 
                                                app.ImageOpenCVgray.shape[::-1], 
                                                None, 
                                                None)

    # Get the dimensions of the image 
    height, width = app.InputImageOpenCV.shape[:2]
     
    # Refine camera matrix
    # Returns optimal camera matrix and a rectangular region of interest
    optimal_camera_matrix, roi = cv2.getOptimalNewCameraMatrix(mtx, dist, 
                                                            (width,height), 
                                                            1, 
                                                            (width,height))

    # Undistort the image
    undistorted_image = cv2.undistort(app.InputImageOpenCV, mtx, dist, None, 
                                    optimal_camera_matrix)

    #undistorted_image1 = cv2.resize(undistorted_image, (int(app.proportion * width * 1.2), int(app.proportion * height * 1.2)))
    cv2.imshow('Image undistored by OpenCV', undistorted_image)      
    return True


###################################################
# Opens input image file
# file_name - Path and name of image to open
###################################################
def OpenImageFile(file_name):
    if file_name[-3:] == 'tif' or file_name[-4:] == 'tiff':
        app.label["text"] = 'Opening image...'
        ResizeImage()
        try:
            app.InputImageTk = Image.open(file_name)
            app.resize_coefficient = app.InputImageTk.width / app.config.max_image_proc_resolution
            if app.resize_coefficient < app.InputImageTk.height / app.config.max_image_proc_resolution:
                app.resize_coefficient = app.InputImageTk.height / app.config.max_image_proc_resolution
            app.resize_coefficient = 1.0 / app.resize_coefficient    
            if app.resize_coefficient > 1.0:
                app.resize_coefficient = 1.0
            #app.resize_coefficient = 1.0    
        except ValueError:
            messagebox.showwarning(title='Error', message='Couldn"t open file')    
            return 1

        app.InputImageOpenCV = cv2.imread(file_name)

        app.ImageOpenCVgray = cv2.cvtColor(app.InputImageOpenCV, cv2.COLOR_BGR2GRAY)  

        if app.resize_coefficient < 1.0:
             #resize big image
            app.InputImageTk = app.InputImageTk.resize((int(app.InputImageTk.width * app.resize_coefficient), int(app.InputImageTk.height * app.resize_coefficient)))
            app.ImageOpenCVgray = cv2.resize(app.ImageOpenCVgray, (0, 0), fx = app.resize_coefficient, fy = app.resize_coefficient)

        app.center_x = app.InputImageTk.width / 2.0
        app.center_y = app.InputImageTk.height / 2.0    

        #app.ImageOpenCVgray = cv2.convertScaleAbs(app.ImageOpenCVgray, alpha=1.5, beta=10)
                        
        #ALLOCATE MEMORY FOR AUXILARY POINTS
        if app.grid_corners == None:
            app.grid_corners = [[[0.0, 0.0] for i in range(app.max_grid_columns)] for j in range(app.max_grid_rows)]
        if app.grid_corners_ == None:
            app.grid_corners_ = [[[0.0, 0.0] for i in range(app.max_grid_columns)] for j in range(app.max_grid_rows)]
        if app.graphical_grid == None:
           app.graphical_grid = [[GraphicalObjects() for i in range(app.max_grid_columns)] for j in range(app.max_grid_rows)]
        app.grid_corner_correct_line_color = [None for i in range(app.max_grid_rows)]
        for r in range(app.max_grid_rows):   
            app.grid_corner_correct_line_color[r] = "#{}{}{}".format(hex(random.randrange(0, 255)).lstrip("0x").zfill(2).upper(), hex(random.randrange(0, 255)).lstrip("0x").zfill(2).upper(), hex(random.randrange(0, 255)).lstrip("0x").zfill(2).upper())
        #print(app.grid_corner_correct_line_color) 

        if app.points_in_line != None:
            for r in range(app.config.target_squares_rows):
                app.points_in_line[r] = 0
        else:
            app.points_in_line = [0 for i in range(app.config.target_squares_rows)]        

        app.lines = 0       
    else:   
        messagebox.showwarning(title='Error', message='Not supported image format')
        return 1

    size = app.InputImageTk.width * app.InputImageTk.height    
    app.UndistoredImage = [0 for i in range(3 * size)] 
    app.DistoredImage = [0 for i in range(3 * size)] 
    pixelMap = app.InputImageTk.load()
    for y in range(app.InputImageTk.height):
        for x in range(app.InputImageTk.width):
            #R
            app.DistoredImage[y * app.InputImageTk.width + x] = pixelMap[x, y][0]
            app.UndistoredImage[y * app.InputImageTk.width + x] = pixelMap[x, y][0]
            #G
            app.DistoredImage[y * app.InputImageTk.width + x + size] = pixelMap[x, y][1]
            app.UndistoredImage[y * app.InputImageTk.width + x + size] = pixelMap[x, y][1]
            #B
            app.DistoredImage[y * app.InputImageTk.width + x + 2 * size] = pixelMap[x, y][2]
            app.UndistoredImage[y * app.InputImageTk.width + x + 2 * size] = pixelMap[x, y][2]
    app.label["text"] = 'Searching for a target...'
    ResizeImage()
    DrawObjects()
    app.canvas.update_idletasks()

    timeout_start_point = time.time()
    timeout = False
    for c in range(app.config.target_squares_columns, 2, -1):
        for r in range(app.config.target_squares_rows, 2, -1):
            GRID = (c, r)

            app.grid_corners_defined, app.OpenCVfoundCorners = cv2.findChessboardCorners(
                app.ImageOpenCVgray, GRID, cv2.CALIB_CB_ADAPTIVE_THRESH + cv2.CALIB_CB_NORMALIZE_IMAGE + cv2.CALIB_CB_FAST_CHECK)
            #print(r, c, app.grid_corners_defined, app.OpenCVfoundCorners)
            if app.grid_corners_defined:
                app.lines = r
                break
            ts = time.time() 
            if ts - timeout_start_point > app.config.target_searching_timeout:
                #target searching timeout
                timeout = True
                break    
        if app.grid_corners_defined:
            for r in range(app.lines):
                app.points_in_line[r] = c
            break
        if timeout:
            break    

    
    if app.grid_corners_defined:
        FineTuneTargetPoints(True)
        for r in range(app.lines):
            for c in range(app.points_in_line[r]):
                app.grid_corners[r][c][0] = np.float64(app.OpenCVfoundCorners[r * app.points_in_line[r] + c][0][0])
                app.grid_corners[r][c][1] = np.float64(app.OpenCVfoundCorners[r * app.points_in_line[r] + c][0][1])
                app.grid_corners_[r][c][0] = np.float64(app.OpenCVfoundCorners[r * app.points_in_line[r] + c][0][0])
                app.grid_corners_[r][c][1] = np.float64(app.OpenCVfoundCorners[r * app.points_in_line[r] + c][0][1])
                #print("x[%d][%d] = %f; y[%d][%d] = %f;" % (r, c, app.grid_corners[r][c][0], r, c, app.grid_corners[r][c][1])) 
    else:
        ManualDefineTargetPoints(False)            
    DrawObjects()                
        
    return 0 

###################################################
# Image resizer
###################################################
def ResizeImage():
    app.canvas.update_idletasks()
    if canvas_width <= 1 or canvas_height <= 1:
        return

    if app.InputImageTk == None:
        return
  
    useHeight = canvas_height < canvas_width
    if useHeight:
        measurement = canvas_height
    else:
        measurement = canvas_width

    if useHeight:
        app.proportion = measurement / app.InputImageTk.height
    else:
        app.proportion = measurement / app.InputImageTk.width

    size = app.InputImageTk.width * app.InputImageTk.height   
    if app.UndistoredImageChanged == True: 
        pixelMap = app.InputImageTk.load()
        for y in range(app.InputImageTk.height):
            for x in range(app.InputImageTk.width):
                pixelMap[x, y] = (app.UndistoredImage[y * app.InputImageTk.width + x], app.UndistoredImage[y * app.InputImageTk.width + x + size], app.UndistoredImage[y * app.InputImageTk.width + x + 2 * size])
            
    app.resizedImage = app.InputImageTk.resize((int(app.InputImageTk.width*app.proportion), int(app.InputImageTk.height*app.proportion)))

    app.resizedImageTk = ImageTk.PhotoImage(app.resizedImage)    

###################################################
# Delete all graphical object to Canvas
###################################################
def DeleteGraphicalObject():
     for r in range(app.config.target_squares_rows):
        for c in range(app.config.target_squares_columns):
            if app.graphical_grid[r][c].grid_corner_circle != None:
                app.canvas.delete(app.graphical_grid[r][c].grid_corner_circle)
                app.graphical_grid[r][c].grid_corner_circle = None
            if app.graphical_grid[r][c].grid_corner_v_line != None:
                app.canvas.delete(app.graphical_grid[r][c].grid_corner_v_line)
                app.graphical_grid[r][c].grid_corner_v_line = None
            if app.graphical_grid[r][c].grid_corner_h_line != None:
                app.canvas.delete(app.graphical_grid[r][c].grid_corner_h_line)
                app.graphical_grid[r][c].grid_corner_h_line = None      
            if app.graphical_grid[r][c].grid_corner_correct_line != None:
                app.canvas.delete(app.graphical_grid[r][c].grid_corner_correct_line)
                app.graphical_grid[r][c].grid_corner_correct_line = None          
   
###################################################
# Draws all graphical object to Canvas
################################################### 
def DrawObjects():
    if app.resizedImageTk == None:
        return

    #update background
    width = app.canvas.winfo_width()
    height = app.canvas.winfo_height()
    if app.background_id != None:
        app.canvas.delete(app.background_id)
    app.background_id = app.canvas.create_rectangle(0, 0, width, height, fill="black")

    #display image
    if app.image_id != None:
        app.canvas.delete(app.image_id)
    
    offset_x = int((width - app.resizedImageTk.width() + 0.5) / 2)
    app.image_id = app.canvas.create_image((offset_x, 0), image=app.resizedImageTk, anchor="nw")
    
    DeleteGraphicalObject()
    #draw grid corners
    for r in range(app.lines):
        for c in range(app.points_in_line[r]):
            x = app.grid_corners_[r][c][0] * app.proportion + offset_x
            y = app.grid_corners_[r][c][1] * app.proportion
            app.graphical_grid[r][c].grid_corner_circle = app.canvas.create_oval(x - width / 400, y - width / 400, x + width / 400, y + width / 400, outline="red", width = 2)
            app.graphical_grid[r][c].grid_corner_v_line = app.canvas.create_line(x, y - width / 400, x, y + width / 400, fill="red", width = 1)
            app.graphical_grid[r][c].grid_corner_h_line = app.canvas.create_line(x - width / 400, y, x + width / 400, y, fill="red", width = 1)
            if c > 0:
                xx = app.grid_corners_[r][c - 1][0] * app.proportion + offset_x
                yy = app.grid_corners_[r][c - 1][1] * app.proportion
                app.graphical_grid[r][c].grid_corner_correct_line = app.canvas.create_line(xx, yy, x, y, fill=app.grid_corner_correct_line_color[r], width = 1, dash=(3,5))


def ManualDefineTargetPoints(reset):
    app.grid_corners_defined = False
    if reset:
        app.label["text"] = 'Target was reset: \n Hold Shift and press left mouse button to set target points, \n Hold Control to start new line, \n Left button double click to apply, right click to reset'      
    else:
        app.label["text"] = 'Target was not found: \n Hold Shift and press left mouse button to set target points, \n Hold Control to start new line, \n Left button double click to apply, right click to reset'              

def FineTuneTargetPoints(found):        
    if found:
        app.label["text"] = 'Target was found: \n Hold left mouse button and move cursor to adjust grid points if necessary, \n Left button double click to apply, right click to reset'      
    else:
        app.label["text"] = 'Hold left mouse button and move cursor to adjust grid points if necessary, \n Left button double click to apply, right click to reset'              


def FatalError(msg):
    messagebox.showwarning(title='Error', message=msg)
    root.destroy()
    os._exit(1)      

###################################################
# Handle to mouse button events
###################################################
def handle_mouse(event):
    if app.image_id == None:
        return
    image_pos = app.canvas.bbox(app.image_id)
    if event.x > image_pos[0] and event.y > image_pos[1] and event.x < image_pos[2] and event.y < image_pos[3]:
        x = event.x - image_pos[0]
        y = event.y
        if "Motion" in str(event):
            button1_pressed = False
            if "Button1" in str(event):
                button1_pressed = True
            if  button1_pressed and app.grid_point_pointed == True: #left button pressed  
                app.grid_corners_[app.grid_point_pointed_coordinates.x][app.grid_point_pointed_coordinates.y][0] = x / app.proportion 
                app.grid_corners_[app.grid_point_pointed_coordinates.x][app.grid_point_pointed_coordinates.y][1] = y / app.proportion
        elif "ButtonPress" in str(event):
            app.grid_point_pointed = False
            if app.grid_corners_defined:
                for r in range(app.lines):
                    for c in range(app.points_in_line[r]):
                        xx = app.grid_corners_[r][c][0] * app.proportion
                        yy = app.grid_corners_[r][c][1] * app.proportion
                        if math.sqrt((xx - x) * (xx - x) + (yy  - y) * (yy  - y)) < 10:
                            app.grid_point_pointed = True
                            app.grid_point_pointed_coordinates.x = r    
                            app.grid_point_pointed_coordinates.y = c
            else:
                new_line = False
                if app.lines == 0 or (app.lines < app.max_grid_rows and "Control" in str(event)):
                    if app.lines:
                        if app.points_in_line[app.lines - 1] <= 2:
                            #in line can't be less than 2 points
                            return
                    app.lines = app.lines + 1
                    new_line = True
                if new_line or "Shift" in str(event):    
                    if app.points_in_line[app.lines - 1] < app.max_grid_columns:
                        app.points_in_line[app.lines - 1] = app.points_in_line[app.lines - 1] + 1    
                        app.grid_corners_[app.lines - 1][app.points_in_line[app.lines - 1] - 1][0] = x / app.proportion     
                        app.grid_corners_[app.lines - 1][app.points_in_line[app.lines - 1] - 1][1] = y / app.proportion    
                        #print("a", app.lines - 1, app.points_in_line[app.lines - 1] - 1, x, y) 
        DrawObjects()   

###################################################
# Handle to mouse double cleck
################################################### 
def handle_mouse_double(event):
    if app.calculation_done:
        root.destroy()
        os._exit(0)  
    if app.grid_corners_defined:
        app.label["text"] = 'Processing....'
        app.canvas.update_idletasks()
        ComputeSolution()
    else:
        if app.lines:
            app.grid_corners_defined = True
            FineTuneTargetPoints(False)
    

###################################################
# Handle to mouse right button releasing
################################################### 
def handle_mouse_rb_released(event):
    #reset grid
    for r in range(app.lines):
        for c in range(app.points_in_line[r]):
            app.grid_corners[r][c][0] = 0
            app.grid_corners[r][c][1] = 0
            app.grid_corners_[r][c][0] = 0
            app.grid_corners_[r][c][1] = 0
        app.points_in_line[r] = 0    
    app.lines = 0
    app.grid_was_reset = True
    ManualDefineTargetPoints(True)            
    DrawObjects()               

###################################################
# Handle application events
###################################################
def handle_configure(event):
    DrawObjects() 
    
# Event handlers 
#---------------------------------------------------    
root.bind("<Configure>", handle_configure)  
app.canvas.bind("<Motion>", handle_mouse)
app.canvas.bind('<ButtonPress-1>', handle_mouse)  
app.canvas.bind('<Double-Button-1>', handle_mouse_double) 
app.canvas.bind('<ButtonRelease-3>', handle_mouse_rb_released) 

#main
#---------------------------------------------------
def Usage():
    print('ldc [options]')
    print('options:')
    print('-i Input image (supported formats .tif .tiff)')
    print('-c Configuration file (supported formats .json)')
    print('-o Output file (supported formats .json)')
    print('-h Help')
#--------------------------------------------------- 
def main(argv):
    app.config = CConfiguration() 

    try:
        opts, args = getopt.getopt(argv,"i:c:o:")
    except getopt.GetoptError:
        Usage()
        root.destroy()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            root.destroy()
            sys.exit(0)
        elif opt in ("-i"):
            app.Input_image_file = arg 
        elif opt in ("-o"):
            app.config.output_file_set = True
            app.config.output_file = arg    
        elif opt in ("-c"):
            if app.config.ParseConfigFile(arg) == 1:
                root.destroy()
                os._exit(1)
            app.label.config(font=(app.config.label_font, app.config.label_font_size))        

    if app.Input_image_file:
        if OpenImageFile(app.Input_image_file) == 1:
            root.destroy()
            os._exit(1)  

def close_window():
    root.destroy()
    os._exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])

 
root.protocol("WM_DELETE_WINDOW", close_window)
root.mainloop()
