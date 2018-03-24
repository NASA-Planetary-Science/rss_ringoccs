# plot_geo.py
#
# plot selected entries from a geometry file (eventually, a geo instance)
# with option to click on/off individual plots, using checkboxes
#
# Revisions:
#   2018 Mar 23 - rfrench - original version

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
from numpy import genfromtxt
import os
import datetime

program = 'plot_geo.py'

infile = '../output/Rev07E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_GEO.TAB'
spm,ret,scet,rkm,londeg,azdeg,Bdeg,D,ripdot,ripazdot,Fkm,impact,X,Y,Z,Xdot,Ydot,Zdot= np.loadtxt(infile,delimiter=',',unpack=True)

t = spm
s0 = X
s1 = Y
s2 = Z
s3 = D

s4 = Xdot
s5 = Ydot
s6 = Zdot

s7 = londeg
s8 = azdeg

s9 = Bdeg
s10 = Fkm

plt.close("all")

fig,ax=plt.subplots(4)
fig.set_size_inches(8.,10.)

l0, = ax[0].plot(t, s0, visible=True, color='r', lw=2)
l1, = ax[0].plot(t, s1, visible=True, color='g', lw=2)
l2, = ax[0].plot(t, s2, visible=True, color='b', lw=2)
l3, = ax[0].plot(t, s3, visible=True, color='k', lw=2)
ax[0].set_ylabel('km')
ax[0].set_title(os.path.basename(infile))
xloc_legend = 1.02
yloc_legend = 0.4
#ax[0].legend((l0,l1,l2,l3),(r'$X$',r'$Y$',r'$Z$',r'$D$'),fontsize='small',loc=[xloc_legend,yloc_legend])

l4, = ax[1].plot(t, s4, visible=True, color='r', lw=2)
l5, = ax[1].plot(t, s5, visible=True, color='g', lw=2)
l6, = ax[1].plot(t, s6, visible=True, color='b', lw=2)
ax[1].set_ylabel('km/s')
#ax[1].legend((l4,l5,l6),(r'$\dot X$',r'$\dot Y$',r'$\dot Z$'),fontsize='small',loc=[xloc_legend,yloc_legend])

l7, = ax[2].plot(t, s7, visible=True, color='r', lw=2)
l8, = ax[2].plot(t, s8, visible=True, color='g', lw=2)
ax[2].set_ylim(-50,360)
ax[2].set_ylabel('degrees')
#ax[2].legend((l7,l8),('Az(deg)','Phi(deg)'),fontsize='small',loc=[xloc_legend,yloc_legend])

l9, = ax[3].plot(t, s9, visible=True, color='r', lw=2)
l10,= ax[3].plot(t, s10,visible=True, color='g', lw=2)
#ax[3].legend((l7,l8),('B(deg)','F(km)'),fontsize='small',loc=[xloc_legend,yloc_legend])

plt.subplots_adjust(left=.27,right=.87)

plt.xlabel('Seconds past midnight')

rax0 = plt.axes([0.05, 0.75, 0.1, 0.1])
plt.title('Position')

#check0 = CheckButtons(rax0, (r' $X$', r' $Y$', r' $Z$', r' $D$'), (True, True, True, True))
check0 = CheckButtons(rax0, (' X', ' Y', ' Z', ' D'), (True, True, True, True))

rax1 = plt.axes([0.05, 0.55, 0.1, 0.1])
#check1 = CheckButtons(rax1, (r' $\dot X$', r' $\dot Y$', r' $\dot Z$'), (True, True, True))
check1 = CheckButtons(rax1, (' X-dot', ' Y-dot', 'Z-dot'), (True, True, True))
plt.title('Velocity')

rax2 = plt.axes([0.05, 0.35, 0.12, 0.1])
check2 = CheckButtons(rax2, ('Az(deg)', 'Phi(deg)'), (True, True))
plt.title('Angles')

rax3 = plt.axes([0.05, 0.15, 0.12, 0.1])
#check3 = CheckButtons(rax3, (r'$B^\circ$', r'$F$(km)'), (True, True))
check3 = CheckButtons(rax3, (' B(deg)', ' F(km)'), (True, True))

def func0(label):
    print(label)
    if label == ' X':
        l0.set_visible(not l0.get_visible())
    elif label == ' Y':
        l1.set_visible(not l1.get_visible())
    elif label == ' Z':
        l2.set_visible(not l2.get_visible())
    elif label == ' D':
        l3.set_visible(not l3.get_visible())
    plt.draw()
def func1(label):
    if label == ' X-dot':
        l4.set_visible(not l4.get_visible())
    elif label == ' Y-dot':
        l5.set_visible(not l5.get_visible())
    elif label == ' Z-dot':
        l6.set_visible(not l6.get_visible())
    plt.draw()
def func2(label):
    print(label)
    if label == 'Az(deg)':
        l7.set_visible(not l7.get_visible())
    elif label == 'Phi(deg)':
        l8.set_visible(not l8.get_visible())
    plt.draw()
def func3(label):
    if label == 'B(deg)':
        l9.set_visible(not l9.get_visible())
    elif label == 'F(km)':
        l10.set_visible(not l10.get_visible())
    plt.draw()
check0.on_clicked(func0)
check1.on_clicked(func1)
check2.on_clicked(func2)
check3.on_clicked(func3)

plotfile = '../figs/plot_geo.jpg'
now = datetime.datetime.now()
str_now = str(now)
plt.show()
ax = fig.add_axes([0,0,1,1])
plt.text(0.5,0.03,str_now[:-7]+'    '+program+' '+os.getcwd()+'/'+plotfile,verticalalignment='bottom',\
        horizontalalignment='center',transform=ax.transAxes, fontsize=8)
ax.set_axis_off()

fig.savefig(plotfile)
