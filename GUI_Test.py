# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 12:29:05 2022

@author: User
"""

from tkinter import * 
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
import numpy as np
import time


###Set initial conditions

plotStart = False
msDelay = 500

### Plot data

def plot_data():
    global plotStart, after_id
    if plotStart:
        y = [np.random.random() for i in range(100)]
        ax1.cla()
        set_axes('ax1')
        ax1.plot(y)
        y = [np.random.random() for i in range(100)]
        ax2.cla()
        set_axes('ax2')
        ax2.plot(y)
        y = [np.random.random() for i in range(100)]
        ax3.cla()
        set_axes('ax3')
        ax3.plot(y)
        y = [np.random.random() for i in range(100)]
        ax4.cla()
        set_axes('ax4')
        ax4.plot(y)
        
        canvas.draw()
        
    after_id = root.after(msDelay,plot_data)
    
def plot_start():
    global plotStart
    plotStart = True
    
def plot_stop():
    global plotStart
    plotStart = False
    
def quit_all():
    """Cancel all scheduled callbacks and quit."""
    root.after_cancel(after_id)
    root.destroy()
    
def set_axes(ax):
    if ax == 'ax1':
        ax1.set_title("Test Data 1")
        ax1.set_xlabel("Time")
        ax1.set_ylabel("Variable 1")
        ax1.set_xlim(0,100)
        ax1.set_ylim(-2,2)
    if ax == 'ax2':
        ax2.set_title("Test Data 2")
        ax2.set_xlabel("Time")
        ax2.set_ylabel("Variable 2")
        ax2.set_xlim(0,100)
        ax2.set_ylim(-2,2)
    if ax == 'ax3':
        ax3.set_title("Test Data 3")
        ax3.set_xlabel("Time")
        ax3.set_ylabel("Variable 3")
        ax3.set_xlim(0,100)
        ax3.set_ylim(-2,2)
    if ax == 'ax4':
        ax4.set_title("Test Data 4")
        ax4.set_xlabel("Time")
        ax4.set_ylabel("Variable 4")
        ax4.set_xlim(0,100)
        ax4.set_ylim(-2,2)


### Main Window

root = Tk()
root.after_cancel(plot_data)
root.title("Real Time Plot")
root.geometry("1280x720")

### Set up main figure widget
fig = Figure()
axesAll = ['ax1','ax2','ax3','ax4']
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)
fig.tight_layout()

for ax in axesAll:
    set_axes(ax)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(side=LEFT,
                            ipadx=100,
                            ipady=100,
                            fill='y')
canvas.draw()


###Buttons
root.update();
start = Button(root,text = "Start",
                  font = ("calibri",12),
                  command = lambda: plot_start())
start.place(x=650,y=1)

root.update();
stop = Button(root,text = "Stop",
                  font = ("calibri",12),
                  command = lambda: plot_stop())
stop.place(x=start.winfo_x()+start.winfo_reqwidth()+20,y=1)

root.update();
quitWindow = Button(root,text = "Quit",
                  font = ("calibri",12),
                  command = lambda: quit_all())
quitWindow.place(x=stop.winfo_x()+stop.winfo_reqwidth()+20,y=1)

### Textboxes




root.after(msDelay,plot_data)
root.protocol('WM_DELETE_WINDOW', quit_all)
root.mainloop()


