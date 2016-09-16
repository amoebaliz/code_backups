import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt

def setup_backend(backend='TkAgg'):
    import sys
    del sys.modules['matplotlib.backends']
    del sys.modules['matplotlib.pyplot']
    import matplotlib as mpl
    mpl.use(backend)  # do this before importing pyplot
    import matplotlib.pyplot as plt
    return plt


my_data=genfromtxt('data.txt', usecols=np.arange(0,2))
my_data = my_data[110:685+1,:]

np.set_printoptions(suppress=True)

soi = my_data[:,1]
soi_la = np.copy(soi); soi_la[soi_la>0]= np.nan
soi_el = np.copy(soi); soi_el[soi_el<0]= np.nan

print np.max(soi)
print np.min(soi)

date = my_data[:,0]

ind = np.arange(len(soi))
width = 1

t_dif = 25
plot_ind = np.arange(t_dif)

y_max = np.ceil(np.max(np.abs(soi)))
y_min = -1*y_max

def animate():
    t_0=0; t_end = t_0 + t_dif
    rect = ax.bar(plot_ind,soi_la[t_0:t_end],width,color='b',edgecolor=[0,0,.5])
    ax.plot([t_dif/2+.5,(t_dif)/2+.5],[-4,4],'--k',linewidth=2.0)
    ax.set_ylim(y_min,y_max)
    ax.set_xlim(0,t_dif)
    ax.text(t_dif-2.5,y_max-.5,'El Nino', color='r',fontsize=14)
    ax.text(t_dif-2.5,y_min+.2,'La Nina', color='b',fontsize=14)
    for i in range(5):
        t_0=i;t_end = t_0+t_dif    
        h = soi_la[t_0:t_end]
        rect.set_height(h)
        fig.canvas.draw()

#fig, ax = plt.subplots(figsize = (10,5))
plt = setup_backend()
fig,ax = plt.subplots()

win = fig.canvas.manager.window
win.after(10, animate)

#ax.bar(plot_ind,soi_la[t_0:t_end],width,color='b',edgecolor=[0,0,.5])
#ax.bar(plot_ind,soi_el[t_0:t_end],width,color='r',edgecolor=[.5,0,0])

#ax.plot([t_dif/2+.5,(t_dif)/2+.5],[-4,4],'--k',linewidth=2.0)
#ax.set_ylim(y_min,y_max)
#ax.set_xlim(0,t_dif)
#ax.text(t_dif-2.5,y_max-.5,'El Nino', color='r',fontsize=14)
#ax.text(t_dif-2.5,y_min+.2,'La Nina', color='b',fontsize=14)

plt.show()



