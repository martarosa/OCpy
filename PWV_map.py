#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:31:46 2020

@author: davide
"""


from scipy.signal import hilbert
from scipy.integrate import trapz, simps
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import auxiliary_functions as af

class PseudoWVD():
    
    def __init__(self,xData,yData):
        self.x=xData
        self.y=yData
        self.already_evaluated=False
        self.time_broadening = None
        self.freq_broadening = None
    
    def pseudoWVD(self, make_analytic=True):
        '''Calculates the Pseudo-Wigner-Ville Distribution using the analytical signal associated to the real-valued input signal,
        in case the input signal is already complex valued choose 'make_analytic=False'
        '''
        if self.already_evaluated!=True:
            max_t=self.x[-1]
            min_t=self.x[0]
            length_t=len(self.x)
            time_axes=np.linspace(min_t,max_t,length_t)
            time_step=time_axes[1]-time_axes[0]    
            
            
            if make_analytic==True:
                y_analytic=hilbert(self.y)
            else:
                y_analytic=self.y
            
            y_analytic_expanded=np.append(np.zeros_like(y_analytic),y_analytic)
            
            
            N=len(y_analytic)
            produkt_matrix=np.ones((len(self.x),2*N),dtype=complex)
            for ti in range(0,len(self.x),1):
                produkt_ti=[]
                for k in range(-N+1,N+1,1):
                    produkt_ti.append(y_analytic_expanded[ti+k]*np.conjugate(y_analytic_expanded[ti-k]))  
                #produkt_matrix[ti]=np.asarray(2*np.asarray(produkt_ti))
                produkt_matrix[ti]=np.asarray(np.asarray(produkt_ti))
        
            
            FFT_matrix=[]
            freq_matrix=[]
            for row in produkt_matrix:
                n=int(len(row))        
                if n%2==0:
                    FFT_row=np.fft.fft(row)[0:int(n/2)]
                    FFT_matrix.append(FFT_row)
                    freq_matrix.append(np.fft.fftfreq(n,time_step)[0:int(n/2)])
                elif n%2==1:
                    FFT_row=np.fft.fft(row)[0:int((n-1)/2)]
                    FFT_matrix.append(FFT_row)
                    freq_matrix.append(np.fft.fftfreq(n,time_step)[0:int((n-1)/2)])
        
            pseudoWVD_matrix=[np.real(row) for row in FFT_matrix]
            #pseudoWVD_freq_axes=np.array(freq_matrix[0])*0.5
            pseudoWVD_freq_axes=np.array(freq_matrix[0])*np.pi
            self.pseudoWVD_freq_axes=pseudoWVD_freq_axes
            self.time_axes=time_axes
            self.pseudoWVD_matrix=np.transpose(pseudoWVD_matrix)
            self.already_evaluated=True
            
     
    def set_time_broadening(self, time_broadening):
        self.time_broadening = time_broadening
        
    def set_freq_broadening(self, freq_broadening):
        self.freq_broadening = freq_broadening
        
    def pseudoWVD_npFT(self):
        '''Calculates the Pseudo-Wigner-Ville Distribution using numpy native implementation of FT without performing Hilbert transform'''
        if self.already_evaluated!=True:
            max_t=self.x[-1]
            min_t=self.x[0]
            length_t=len(self.x)
            time_axes=np.linspace(min_t,max_t,length_t)
            time_step=time_axes[1]-time_axes[0]    
        
            y_expanded=np.append(np.zeros_like(self.y),self.y)
        
        
            N=len(self.y)
            produkt_matrix=np.ones((len(self.x),2*N),dtype=complex)
        
            for ti in range(0,len(self.x),1):
                produkt_ti=[]
                for k in range(-N+1,N+1,1):
                    produkt_ti.append(y_expanded[ti+k]*np.conjugate(y_expanded[ti-k]))  
            #produkt_matrix[ti]=np.asarray(2*np.asarray(produkt_ti))
                produkt_matrix[ti]=np.asarray(np.asarray(produkt_ti))
    
        
            FFT_matrix=[]
            freq_matrix=[]
            for row in produkt_matrix:
                n=int(len(row))        
                if n%2==0:
                    FFT_row=np.fft.fft(row)[0:int(n/2)]
                    FFT_matrix.append(FFT_row)
                    freq_matrix.append(np.fft.fftfreq(n,time_step)[0:int(n/2)])
                elif n%2==1:
                    FFT_row=np.fft.fft(row)[0:int((n-1)/2)]
                    FFT_matrix.append(FFT_row)
                    freq_matrix.append(np.fft.fftfreq(n,time_step)[0:int((n-1)/2)])
                    
                    
            pseudoWVD_matrix=[np.real(row) for row in FFT_matrix]
            #pseudoWVD_freq_axes=np.array(freq_matrix[0])*0.5
            pseudoWVD_freq_axes=np.array(freq_matrix[0])*np.pi
            self.pseudoWVD_freq_axes=pseudoWVD_freq_axes
            self.time_axes=time_axes
            self.pseudoWVD_matrix=np.transpose(pseudoWVD_matrix)
            self.already_evaluated=True
            
            
    
    def pseudoWVD_differentFT(self):
        '''Calculates the Pseudo-Wigner-Ville Distribution using the analytical signal associated to the real-valued input signal,
        in case the input signal is already complex valued choose 'make_analytic=False'
        '''
        if self.already_evaluated!=True:
            max_t=self.x[-1]
            min_t=self.x[0]
            length_t=len(self.x)
            time_axes=np.linspace(min_t,max_t,length_t)
            time_step=time_axes[1]-time_axes[0]    
            
            y_expanded=np.append(np.zeros_like(self.y),self.y)
            
            
            N=len(self.y)
            produkt_matrix=np.ones((len(self.x),2*N),dtype=complex)
            
            for ti in range(0,len(self.x),1):
                produkt_ti=[]
                for k in range(-N+1,N+1,1):
                    produkt_ti.append(y_expanded[ti+k]*np.conjugate(y_expanded[ti-k]))  
                #produkt_matrix[ti]=np.asarray(2*np.asarray(produkt_ti))
                produkt_matrix[ti]=np.asarray(np.asarray(produkt_ti))
        
            
            FFT_matrix=[]
            freq_matrix=[]
            for row in produkt_matrix:
                n=int(len(row))        
               # FFT_zeile=1/np.sqrt(n)*np.fft.fft(row)[0:int(n/2)]
                FFT_row = af.smooth_fft_field(time_step, n, 0, row)[1][0:int(n/2)]
                #FFT_zeile=np.fft.fft(row)[0:int(n/2)]
                FFT_matrix.append(FFT_row)
                freq_matrix.append(af.smooth_fft_field(time_step, n, 0, row)[0][0:int(n/2)])
            
            
            pseudoWVD_matrix=FFT_matrix
            #pseudoWVD_freq_axes=np.array(freq_matrix[0])*0.5
            pseudoWVD_freq_axes=np.array(freq_matrix[0])*0.5
            self.pseudoWVD_freq_axes=pseudoWVD_freq_axes
            self.time_axes=time_axes
            self.pseudoWVD_matrix=np.transpose(pseudoWVD_matrix)
            self.already_evaluated=True
    
    
    def gaussian_filter(self, variable, mean, std):
        return np.exp((-np.power(variable-mean,2))/(2*std))
            
    def apply_gaussian_filter(self): 
        self.smoothed_map = np.zeros_like(self.pseudoWVD_matrix)
        for t in range(len(self.time_axes)):
            time_filter = self.gaussian_filter(self.time_axes[t], self.time_axes, self.time_broadening)
            for omega in range(len(self.pseudoWVD_freq_axes)):
                frequency_filter = self.gaussian_filter(self.pseudoWVD_freq_axes[omega], self.pseudoWVD_freq_axes, self.freq_broadening)
                to_integrate = np.outer(frequency_filter, np.transpose(time_filter)) @ self.pseudoWVD_matrix
           #     first_integration = np.trapz(to_integrate, self.time_axes, a)
                self.smoothed_map[t, omega] = simps(simps(to_integrate, self.time_axes, axis=0), self.pseudoWVD_freq_axes, axis=0) 
                
        
    def contour_plot(self, fig_title="default fig_title",axes_titles=["default_x","default_y"],colormap=cm.viridis):
        '''Produces a standard filled contour-plot using matplotlibs contourf
        '''
        plt.figure()
        plt.title(fig_title)
        plt.xlabel(axes_titles[0])
        plt.ylabel(axes_titles[1])
        plt.contourf(self.time_axes,self.pseudoWVD_freq_axes,abs(self.pseudoWVD_matrix),1000,cmap=colormap)
        plt.xlim(self.time_axes[0],self.time_axes[-1])
        plt.ylim(self.pseudoWVD_freq_axes[0],self.pseudoWVD_freq_axes[-1])
      #  plt.clim(vmin = 0, vmax = 1)
        plt.colorbar()
        plt.show()
        
        


