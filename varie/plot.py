import matplotlib.pyplot as plt
import tkinter as tk
import GUI_system as sy
import numpy as np
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
from matplotlib import ticker
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class Process_box(tk.Toplevel):
    def __init__(self, figure):
        tk.Toplevel.__init__(self)
        self.geometry('400x400')
        self.title("Plot_field")
        self.canv_field = FigureCanvasTkAgg(figure, self)
        self.canv_field.show()
        self.canv_field.get_tk_widget().pack()
        self.canv_field.pack()
        #self.grab_set()
        #self.attributes('-topmost', 'true')






class Plot_field():
    def __init__(self):

        self.tau=None
        self.x_lim=None
        self.y_lim=None
        self.xtics=None
        self.xaxis=None





        self.fig = plt.Figure()
        self.fieldx = self.fig.add_subplot(311)
        self.fieldy = self.fig.add_subplot(312, sharex=self.fieldx)
        self.fieldz = self.fig.add_subplot(313, sharex=self.fieldx)


    def set_xaxis(self,nstep):
        self.xaxis = np.arange(0, self.tau, self.tau/nstep)


    def set_tau(self,nstep,dt):
        self.tau=nstep*dt

    def plot(self,field):
        for axis in ['top','bottom','left','right']:
            self.fieldx.spines[axis].set_linewidth(0.3)
        self.fieldx.set_xlabel('time [a.u.]', fontsize=5,color='black')
        self.fieldx.set_ylabel('field', fontsize=5, color='black',)
        self.fieldx.tick_params(axis='both', labelsize=4, width=0.3,length=2)

        #self.x_lim=[0,self.tau]
        #self.xtics=np.arange(0, self.x_lim[-1]+1,10)
        #self.fieldx.set_xlim(self.x_lim)
        #self.fieldx.xaxis.set_ticks(self.xtics)
        #self.fieldx.yaxis.set_ticks(self.ytics)
        self.fieldx.cla()
        self.fieldy.cla()
        self.fieldz.cla()
        self.fieldx.plot(self.xaxis,field[:,0])
        self.fieldy.plot(self.xaxis, field[:,1])
        self.fieldz.plot(self.xaxis, field[:,2])











class Plot_final_pop():
    def __init__(self):
        self.x_lim=None
        self.y_lim=None
        self.xtics=None
        self.ytics=None
        self.fig=None
        self.lines=None

        self.fig= plt.Figure(figsize=(1.4,1.2),dpi=300)
        self.ax = self.fig.add_subplot(111)
        self.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
        for axis in ['top','bottom','left','right']:
            self.ax.spines[axis].set_linewidth(0.3)
        self.ax.set_xlabel('iteration [a.u.]', fontsize=5,color='black')
        self.ax.set_ylabel('population', fontsize=5, color='black',)
        self.ax.tick_params(axis='both', labelsize=4, width=0.3,length=2)


    def plot(self,x,vect_y,index_states):
        self.initialize_plot()
        self.ax.plot(x,vect_y[:,index_states],'--o', markersize=2)



    def empty_plot(self):
        self.ax.set_xlim([0,1])
        self.ax.set_ylim([0,1])
        self.ax.xaxis.set_ticks([0])
        self.ax.yaxis.set_ticks([0])
        self.ax.plot(0, 0, 'o' ,markersize=0.001)



    def initialize_plot(self):
        self.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
        for axis in ['top','bottom','left','right']:
            self.ax.spines[axis].set_linewidth(0.3)
        self.ax.set_xlabel('iteration [a.u.]', fontsize=5,color='black')
        self.ax.set_ylabel('population', fontsize=5, color='black',)
        self.ax.tick_params(axis='both', labelsize=4, width=0.3,length=2)
        self.x_lim=[0,sy.system.oc.iterations]
        self.y_lim=[0,1]
        self.ytics=[0,0.25,0.50,0.75,1]
        self.xtics=np.arange(0, self.x_lim[-1]+1,10)
        self.ax.set_xlim(self.x_lim)
        self.ax.set_ylim(self.y_lim)
        self.ax.xaxis.set_ticks(self.xtics)
        self.ax.yaxis.set_ticks(self.ytics)