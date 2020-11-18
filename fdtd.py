# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 16:39:40 2017

@author: Owner
"""

# -*- coding: utf-8 -*-
"""
Created on Sun May  7 20:08:42 2017

@author: Owner

"""
import math
import sys
import numpy as np
import matplotlib.pyplot as plt

#### original module ####
import file_selecter as fs

class FDTD:
    def __init__(self):
        #initial condition
        self.t = 0.

        self.set_physical_const()
        self.set_simulation_param()
        self.set_mesh_num()
        #self.set_pml()

        self.set_start_idx(1)
        self.set_end_idx()

        self.make_diff_table(self.n_w, self.n_w2, self.n_L, self.n_H)
        self.make_eps_table(self.eps_0)

        self.set_differential()
        self.set_E_H_for_plot()
        self.set_cfl_condition()

        print("Initializing Done.")

    def run(self):
        self.struct_eh_matrix()
        self.struct_yee_grid()
        self.solve()
        return

    def set_physical_const(self):
        #const
        self.c = 3.000 * 10**8
        self.eps_0 = 8.854 * 10**(-12)
        self.mu_0 = 1/(self.c**2 * self.eps_0) # VacuumPermeability

    def set_simulation_param(self):
        #parameter
        self.n_co = 3.6
        self.n_cl = 3.24
        self.w = 0.3 * 10**(-6)
        self.n_H = 1.2
        self.n_L = 1.
        self.lambda_0 = 6.0 *10**(-6) # wavelength
        self.omega_0 = 2.*math.pi*self.c/self.lambda_0
        self.area_x = 30.0 * 10**(-6) # and -50.0 um
        self.area_y = 30.0 * 10**(-6) # and -50.0 um
        self.area_z = 30.0 * 10**(-6) # and 100.0 um
        self.t_max = 1.0 * 10**(-9)
        self.m = 0.8
        self.R_0 = 0.01 *10**(-2)
        self.mu = self.mu_0 #Not magnetic body

    def set_mesh_num(self):
        #mesh
        self.n_meshx = 30
        self.n_meshy = 30
        self.n_meshz = 90
        self.n_w = 5
        self.n_w2= 30

    def set_pml(self):
        self.d_pml_x = self.n_meshx/10.
        self.d_pml_y = self.n_meshy/10.
        self.d_pml_z = self.n_meshz/10.
    
    def set_start_idx(self, start):
        #Start & End
        self.xs = start
        self.ys = start
        self.zs = start

    def set_end_idx(self):
        self.xe = self.n_meshx - self.xs
        self.ye = self.n_meshy - self.ys
        self.ze = self.n_meshz - self.zs

    def make_diff_table(self, n_w, n_w2, n_L, n_H):
        n_meshx = self.n_meshx
        n_meshy = self.n_meshy
        n_meshz = self.n_meshz
        n_w2 = self.n_w2
        n_L = self.n_L
        n_H = self.n_H
        self.diff_table = np.zeros([n_meshx, n_meshy, n_meshz])
        self.diff_table[0:n_meshx, 0:n_meshy, 0:n_meshz-n_w2] = n_L 
        self.diff_table[0:n_meshx, 0:n_meshy, n_meshz-n_w2:n_meshz] = n_H

    def make_eps_table(self, eps_0):
        n_meshx = self.n_meshx
        n_meshy = self.n_meshy
        n_meshz = self.n_meshz
        self.eps = np.zeros([n_meshx, n_meshy, n_meshz])
        self.eps[0:n_meshx, 0:n_meshy, 0:n_meshz] = np.power(self.diff_table[0:n_meshx, 0:n_meshy, 0:n_meshz], 2) * self.eps_0       

    def set_differential(self):
        n_meshx = self.n_meshx
        n_meshy = self.n_meshy
        n_meshz = self.n_meshz
        area_x = self.area_x
        area_y = self.area_y
        area_z = self.area_z
        #dt = 1.0 * 10**(-18) #CFLでdt決めるからここでは不要
        self.dx = 2*area_x/n_meshx
        self.dy = 2*area_y/n_meshy
        self.dz = 2*area_z/n_meshz 

    def set_E_H_for_plot(self):
        n_meshx = self.n_meshx
        n_meshy = self.n_meshy
        n_meshz = self.n_meshz
        #for display
        self.H_ = np.zeros([n_meshx, n_meshy, n_meshz])
        self.E_max = []
        self.H_max = []
        self.Exy_for_plt = np.zeros([n_meshx, n_meshy])
        self.Exz_for_plt = np.zeros([n_meshx, n_meshz])
        self.Eyz_for_plt = np.zeros([n_meshy, n_meshz])

        #flag
        self.flag_plt = False
        self.flag_start = False
        
    def set_cfl_condition(self):
        #CFL condition
        #v*dt <= 1./sqrt((1/dx)^2+(1/dy)^2+(1/dz)^2)
        dx = self.dx
        dy = self.dy
        dz = self.dz
        n_min = min([self.n_co, self.n_cl, self.n_H, self.n_L])
        v_max = self.c/n_min
        self.dt = 1/5.*1/math.sqrt((1/dx)**2+(1/dy)**2+(1/dz)**2)/v_max

        if not(v_max * self.dt <= 1/math.sqrt((1/dx)**2+(1/dy)**2+(1/dz)**2)):
            raise "Runtime Error: Parametars don't meet the CFL condition."

    def struct_eh_matrix(self):
        #Structing E&H matrix
        n_meshx = self.n_meshx
        n_meshy = self.n_meshy
        n_meshz = self.n_meshz
        self.E_0 = np.zeros([n_meshx, n_meshy, n_meshz])
        self.H_0 = np.zeros([n_meshx, n_meshy, n_meshz])

    def struct_pml(self):
        #Structing PML
        d_pml_x = self.d_pml_x
        m = self.m
        sigma_x_max = -(m+1)*math.ln(self.R_0)/(2*self.eta_0*self.eps_0*d_pml_x) # It shoud be modified later.
        for xi_x in d_pml_x: # 謎コード。勉強しなおそう
            sigma = (xi_x/d_pml_x)**m * sigma_x_max

        print("Initializing Done.")    

    def struct_yee_grid(self):
        """
        Yee Grid (Staggered grid)

        E[n-1] H[n-1/2] E[n] H[n+1/2] E[n+1]
        E[0               1             2     ]
        H[       0             1              ]
        E_n[t, x, y, z]
        t= n*dt
        E_upd = f(E)
        """
        n_meshx = self.n_meshx
        n_meshy = self.n_meshy
        n_meshz = self.n_meshz
        self.E_nx = self.H_0
        self.E_ny = self.H_0
        self.E_nz = self.H_0
        self.H_nx = self.E_0
        self.H_ny = self.E_0
        self.H_nz = self.E_0
        self.E_nx1 = np.zeros([n_meshx, n_meshy, n_meshz])
        self.E_ny1 = np.zeros([n_meshx, n_meshy, n_meshz])
        self.E_nz1 = np.zeros([n_meshx, n_meshy, n_meshz])
        self.H_nx1 = np.zeros([n_meshx, n_meshy, n_meshz])
        self.H_ny1 = np.zeros([n_meshx, n_meshy, n_meshz])
        self.H_nz1 = np.zeros([n_meshx, n_meshy, n_meshz])

        #for i in range(20, 30):
        #    E_nx[i,25,25] = 1
        #    E_nx[i,24,25] = 1
        #    H_ny1[i+1,25,25] = 0.5
        #    H_ny1[i+1,24,25] = 0.5
        #E_nx[25,25,25] = 1.

        self.E_nx[self.xs:self.xe, self.ys:self.ye, self.zs:self.ze] = 1.
        #H_ny[25,25,25]=E_nx[25,25,25] /pow(mu/eps, 1/2.)
        #Gaussian
        for i in range(1, n_meshx-1): #Calc new E
            for j in range(1, n_meshy-1):
                for k in range(1, n_meshz-1):  
                    xc = n_meshx/2
                    yc = n_meshy/2
                    zc = n_meshz/2
                    x = i 
                    y = j
                    z = k
                    #r2 = (x - xc)**2 + (y - yc)**2 + (z - zc)**2
                    w2 = 10.
                    #E_nx[i, j, k] = math.exp(-r2 / w2)
                    self.E_nx[i, j, k] = math.exp(-((x - xc)**2 + (y - yc)**2 + (z  - zc)**2) / w2)
                    self.E_ny[i, j, k] = math.exp(-((x - xc)**2 + (y - yc)**2 + (z - zc)**2) / w2)
                    self.E_nz[i, j, k] = math.exp(-((x - xc)**2 + (y - yc)**2 + (z - zc)**2) / w2)
        #    

    def solve(self):
        #Difference method (FDTD)
        dt = self.dt
        for t in range(1,int(self.t_max/dt)-1):
            print("n = "+ str(t))
            print("t = "+ str(round(t*dt*10**15,3)) + " [fs]")
            self.set_forced_occilation(t)
            self.set_bc_mur1(t)
            self.update(t)
            self.plot_result(t)
        
    def set_forced_occilation(self, t):
        #Forced Oscillation
    #    for i in range(n_w, n_meshx-1-n_w): #Calc new E
    #          for j in range(n_w, n_meshy-1-n_w):
    #              E_nx[i, j, int(n_meshz/2.)] = math.sin(omega_0*t*dt)
        n_meshx = self.n_meshx
        n_meshy = self.n_meshy
        n_meshz = self.n_meshz
        xs = self.xs
        xe = self.xe
        ys = self.ys
        ye = self.ye
        zs = self.zs
        ze = self.ze
        dt = self.dt
        dx = self.dx
        dy = self.dy
        dz = self.dz
        eps = self.eps
        self.E_nx[self.n_w:xe-self.n_w, self.n_w:ye-self.n_w, int(n_meshz/2.)] = np.sin(self.omega_0*t*dt)
        
        self.E_nx1[xs:xe, ys:ye, zs:ze] = self.E_nx[xs:xe, ys:ye, zs:ze] \
                + dt / (eps[xs:xe, ys:ye, zs:ze]*dy) * (self.H_nz[xs:xe, ys:ye, zs:ze] - self.H_nz[xs:xe, ys-1:ye-1, zs:ze]) \
                - dt / (eps[xs:xe, ys:ye, zs:ze]*dz) * (self.H_ny[xs:xe, ys:ye, zs:ze] - self.H_ny[xs:xe, ys:ye, zs-1:ze-1]) 
                
        self.E_ny1[xs:xe, ys:ye, zs:ze] = self.E_ny[xs:xe, ys:ye, zs:ze] \
                + dt / (eps[xs:xe, ys:ye, zs:ze]*dz) * (self.H_nx[xs:xe, ys:ye, zs:ze] - self.H_nx[xs:xe, ys:ye, zs-1:ze-1]) \
                - dt / (eps[xs:xe, ys:ye, zs:ze]*dx) * (self.H_nz[xs:xe, ys:ye, zs:ze] - self.H_nz[xs-1:xe-1, ys:ye, zs:ze]) 
                
        self.E_nz1[xs:xe, ys:ye, zs:ze] = self.E_nz[xs:xe, ys:ye, zs:ze] \
                + dt / (eps[xs:xe, ys:ye, zs:ze]*dx) * (self.H_ny[xs:xe, ys:ye, zs:ze] - self.H_ny[xs-1:xe-1, ys:ye, zs:ze]) \
                - dt / (eps[xs:xe, ys:ye, zs:ze]*dy) * (self.H_nx[xs:xe, ys:ye, zs:ze] - self.H_nx[xs:xe, ys-1:ye-1, zs:ze])       
    
    def set_bc_mur1(self, t):
        #Boundary Condition (Mur 1) for E-field
        n_meshx = self.n_meshx
        n_meshy = self.n_meshy
        n_meshz = self.n_meshz
        xs = self.xs
        xe = self.xe
        ys = self.ys
        ye = self.ye
        zs = self.zs
        ze = self.ze
        c = self.c
        dt = self.dt
        dx = self.dx
        dy = self.dy
        dz = self.dz        

        self.E_nx1[0, ys:ye, zs:ze] = self.E_nx[1, ys:ye, zs:ze] +(c*dt-dx)/(c*dt+dx)*(self.E_nx1[1, ys:ye, zs:ze] - self.E_nx[0, ys:ye, zs:ze])
        self.E_ny1[0, ys:ye, zs:ze] = self.E_ny[1, ys:ye, zs:ze] +(c*dt-dx)/(c*dt+dx)*(self.E_ny1[1, ys:ye, zs:ze] - self.E_ny[0, ys:ye, zs:ze])
        self.E_nz1[0, ys:ye, zs:ze] = self.E_nz[1, ys:ye, zs:ze] +(c*dt-dx)/(c*dt+dx)*(self.E_nz1[1, ys:ye, zs:ze] - self.E_nz[0, ys:ye, zs:ze])
        self.E_nx1[n_meshx-1, ys:ye, zs:ze] = self.E_nx[n_meshx-2, ys:ye, zs:ze] +(c*dt-dx)/(c*dt+dx)*(self.E_nx1[n_meshx-2, ys:ye, zs:ze] - self.E_nx[n_meshx-1, ys:ye, zs:ze])
        self.E_ny1[n_meshx-1, ys:ye, zs:ze] = self.E_ny[n_meshx-2, ys:ye, zs:ze] +(c*dt-dx)/(c*dt+dx)*(self.E_ny1[n_meshx-2, ys:ye, zs:ze] - self.E_ny[n_meshx-1, ys:ye, zs:ze])
        self.E_nz1[n_meshx-1, ys:ye, zs:ze] = self.E_nz[n_meshx-2, ys:ye, zs:ze] +(c*dt-dx)/(c*dt+dx)*(self.E_nz1[n_meshx-2, ys:ye, zs:ze] - self.E_nz[n_meshx-1, ys:ye, zs:ze])
        
        self.E_nx1[xs:xe, 0, zs:ze] = self.E_nx[xs:xe, 1, zs:ze] +(c*dt-dy)/(c*dt+dy)*(self.E_nx1[xs:xe, 1, zs:ze] - self.E_nx[0, ys:ye, zs:ze])
        self.E_ny1[xs:xe, 0, zs:ze] = self.E_ny[xs:xe, 1, zs:ze] +(c*dt-dy)/(c*dt+dy)*(self.E_ny1[xs:xe, 1, zs:ze] - self.E_ny[0, ys:ye, zs:ze])
        self.E_nz1[xs:xe, 0, zs:ze] = self.E_nz[xs:xe, 1, zs:ze] +(c*dt-dy)/(c*dt+dy)*(self.E_nz1[xs:xe, 1, zs:ze] - self.E_nz[0, ys:ye, zs:ze])
        self.E_nx1[xs:xe, n_meshy-1, zs:ze] = self.E_nx[xs:xe, n_meshy-2, zs:ze] +(c*dt-dy)/(c*dt+dy)*(self.E_nx1[xs:xe, n_meshy-2, zs:ze] - self.E_nx[n_meshy-1, ys:ye, zs:ze])
        self.E_ny1[xs:xe, n_meshy-1, zs:ze] = self.E_ny[xs:xe, n_meshy-2, zs:ze] +(c*dt-dy)/(c*dt+dy)*(self.E_ny1[xs:xe, n_meshy-2, zs:ze] - self.E_ny[n_meshy-1, ys:ye, zs:ze])
        self.E_nz1[xs:xe, n_meshy-1, zs:ze] = self.E_nz[xs:xe, n_meshy-2, zs:ze] +(c*dt-dy)/(c*dt+dy)*(self.E_nz1[xs:xe, n_meshy-2, zs:ze] - self.E_nz[n_meshy-1, ys:ye, zs:ze])

        self.E_nx1[xs:xe, ys:ye, 0] = self.E_nx[xs:xe, ys:ye, 1] +(c*dt-dz)/(c*dt+dz)*(self.E_nx1[xs:xe, ys:ye, 1] - self.E_nx[xs:xe, ys:ye, 0])
        self.E_ny1[xs:xe, ys:ye, 0] = self.E_ny[xs:xe, ys:ye, 1] +(c*dt-dz)/(c*dt+dz)*(self.E_ny1[xs:xe, ys:ye, 1] - self.E_ny[xs:xe, ys:ye, 0])
        self.E_nz1[xs:xe, ys:ye, 0] = self.E_nz[xs:xe, ys:ye, 1] +(c*dt-dz)/(c*dt+dz)*(self.E_nz1[xs:xe, ys:ye, 1] - self.E_nz[xs:xe, ys:ye, 0])
        self.E_nx1[xs:xe, ys:ye, n_meshz-1] = self.E_nx[xs:xe, ys:ye, n_meshz-2] +(c*dt-dz)/(c*dt+dz)*(self.E_nx1[xs:xe, ys:ye, n_meshz-2] - self.E_nx[xs:xe, ys:ye, n_meshy-1])
        self.E_ny1[xs:xe, ys:ye, n_meshz-1] = self.E_ny[xs:xe, ys:ye, n_meshz-2] +(c*dt-dz)/(c*dt+dz)*(self.E_ny1[xs:xe, ys:ye, n_meshz-2] - self.E_ny[xs:xe, ys:ye, n_meshy-1])
        self.E_nz1[xs:xe, ys:ye, n_meshz-1] = self.E_nz[xs:xe, ys:ye, n_meshz-2] +(c*dt-dz)/(c*dt+dz)*(self.E_nz1[xs:xe, ys:ye, n_meshz-2] - self.E_nz[xs:xe, ys:ye, n_meshy-1])
    
    def update(self, t):
        n_meshx = self.n_meshx
        n_meshy = self.n_meshy
        n_meshz = self.n_meshz
        xs = self.xs
        xe = self.xe
        ys = self.ys
        ye = self.ye
        zs = self.zs
        ze = self.ze
        dt = self.dt

        ### update E
        self.E_nx = self.E_nx1
        self.E_ny = self.E_ny1
        self.E_nz = self.E_nz1
        
        self.calc_maxwell_equation()

        #Update
        self.H_nx = self.H_nx1
        self.H_ny = self.H_ny1
        self.H_nz = self.H_nz1
        
    def calc_maxwell_equation(self):
        xs = self.xs
        xe = self.xe
        ys = self.ys
        ye = self.ye
        zs = self.zs
        ze = self.ze
        dt = self.dt
        dx = self.dx
        dy = self.dy
        dz = self.dz        
        mu = self.mu

        self.H_nx1[xs:xe, ys:ye, zs:ze] = self.H_nx[xs:xe, ys:ye, zs:ze] \
                + dt / (mu*dz) * (self.E_ny[xs:xe, ys:ye, zs+1:ze+1] - self.E_ny[xs:xe, ys:ye, zs:ze]) \
                - dt / (mu*dy) * (self.E_nz[xs:xe, ys+1:ye+1, zs:ze] - self.E_nz[xs:xe, ys:ye, zs:ze]) 
                
        self.H_ny1[xs:xe, ys:ye, zs:ze] = self.H_ny[xs:xe, ys:ye, zs:ze] \
                + dt / (mu*dx) * (self.E_nz[xs+1:xe+1, ys:ye, zs:ze] - self.E_nz[xs:xe, ys:ye, zs:ze]) \
                - dt / (mu*dz) * (self.E_nx[xs:xe, ys:ye, zs+1:ze+1] - self.E_nx[xs:xe, ys:ye, zs:ze]) 
                
        self.H_nz1[xs:xe, ys:ye, zs:ze] = self.H_nz[xs:xe, ys:ye, zs:ze] \
                + dt / (mu*dy) * (self.E_nx[xs:xe, ys+1:ye+1, zs:ze] - self.E_nx[xs:xe, ys:ye, zs:ze]) \
                - dt / (mu*dx) * (self.E_ny[xs+1:xe+1, ys:ye, zs:ze] - self.E_ny[xs:xe, ys:ye, zs:ze]) 

    def plot_result(self, t):
        n_meshx = self.n_meshx
        n_meshy = self.n_meshy
        n_meshz = self.n_meshz
        xs = self.xs
        xe = self.xe
        ys = self.ys
        ye = self.ye
        zs = self.zs
        ze = self.ze
        dt = self.dt
        dx = self.dx
        dy = self.dy
        dz = self.dz     
        area_x = self.area_x
        area_y = self.area_y
        area_z = self.area_z

        #Plotting
        self.Exy_for_plt[0:n_meshx, 0:n_meshy] = np.power( \
                np.power(self.E_nx[0:n_meshx, 0:n_meshy, int(n_meshz/2.)], 2.) + \
                np.power(self.E_ny[0:n_meshx, 0:n_meshy, int(n_meshz/2.)], 2.) + \
                np.power(self.E_nz[0:n_meshx, 0:n_meshy, int(n_meshz/2.)], 2.) , 1./2.) 
        self.Exz_for_plt[0:n_meshx, 0:n_meshz] = np.power( \
                np.power(self.E_nx[0:n_meshx, int(n_meshy/2.), 0:n_meshz], 2.) + \
                np.power(self.E_ny[0:n_meshx, int(n_meshy/2.), 0:n_meshz], 2.) + \
                np.power(self.E_nz[0:n_meshx, int(n_meshy/2.), 0:n_meshz], 2.) , 1./2.)
        self.Eyz_for_plt[0:n_meshy, 0:n_meshz] = np.power( \
                np.power(self.E_nx[int(n_meshx/2.), 0:n_meshy, 0:n_meshz], 2.) + \
                np.power(self.E_ny[int(n_meshx/2.), 0:n_meshy, 0:n_meshz], 2.) + \
                np.power(self.E_nz[int(n_meshx/2.), 0:n_meshy, 0:n_meshz], 2.) , 1./2.)   
        self.E_max.append([t*dt*10**15, max(self.Exy_for_plt.max(), self.Exz_for_plt.max(), self.Eyz_for_plt.max())]) 
        self.H_[xs:xe, ys:ye, zs:ze] = np.power( np.power(self.H_nx[xs:xe, ys:ye, zs:ze], 2.) + \
                                            np.power(self.H_ny[xs:xe, ys:ye, zs:ze], 2.) + \
                                                np.power(self.H_nz[xs:xe, ys:ye, zs:ze], 2.), 1./2.)  

        self.H_max.append([t*dt*10**15, np.max(self.H_)*1000]) 
        print("Emax = "+str(self.Exy_for_plt.max()))
        #Plotting
        x_plt = np.arange(-(area_x)*10**6, area_x*10**6, dx*10**6)
        y_plt = np.arange(-(area_y)*10**6, area_y*10**6, dy*10**6)
        z_plt = np.arange(-(area_z)*10**6, area_z*10**6, dz*10**6)
        self.E_max_plt = np.array(self.E_max)
        self.H_max_plt = np.array(self.H_max)
        
        Y1, X1 = np.meshgrid(y_plt, x_plt)
        Z, X2 = np.meshgrid(z_plt, x_plt)
        Z, Y = np.meshgrid(z_plt, y_plt)
        
        if self.flag_plt == False:
            plt.figure(figsize =(12, 10))
        
        plt.subplot(2,2,1)
        plt.pcolor(X1, Y1, self.Exy_for_plt)
        plt.title('|E (X-Y plane)|')
        plt.xlabel("y [um]")
        plt.ylabel("x [um]")
        plt.colorbar()
        
        plt.subplot(2,2,2)
        plt.pcolor(X2, Z, self.Exz_for_plt)
        plt.title('|E (X-Z plane)|')
        plt.xlabel("x [um]")
        plt.ylabel("z [um]")
        plt.colorbar()
        
        plt.subplot(2,2,3)
        plt.pcolor(Y, Z, self.Eyz_for_plt)
        plt.title('|E (Y-Z plane)|')
        plt.xlabel("y [um]")
        plt.ylabel("z [um]")
        plt.colorbar()
        
        plt.subplot(2,2,4)
        plt.plot(self.E_max_plt[:, 0], self.E_max_plt[:, 1])
        plt.plot(self.E_max_plt[:, 0], self.H_max_plt[:, 1])
        plt.title('Emax (t)')
        plt.xlabel("t [fs]")
        plt.ylabel("E [V/m]")
        
        if self.flag_start == False:
            command = input("Press Enter, then start.")
            self.flag_start = True
            
        plt.pause(0.001)
        plt.clf()
        self.flag_plt=True


def main():
    fdtd = FDTD()
    fdtd.run()


if __name__ == "__main__":
    main()