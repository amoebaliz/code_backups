# FORCING THE COLORBAR TO A GIVEN RANGE

## PCOLOR EXAMPLE ##
import matplotlib.pylab as plt
import numpy as np
A = np.random.randint(2,20,size =(30,30)) 
B = np.random.randint(10,50,size =(30,30))

plt.subplot(121)
plt.pcolor(A,vmin=0,vmax=20)
# ALT: plt.clim(0,20)
plt.colorbar()

plt.subplot(122)
plt.pcolor(B,vmin=0,vmax=20)
# ALT:plt.clim(0,20)
plt.colorbar()


## CONTOURF EXAMPLE ##
v = np.linspace(0,20,11,endpoint=True)

plt.subplot(121)
plt.contourf(A,v)
plt.clim = (0,20)
cbar1 = plt.colorbar(ticks=v)


plt.subplot(122)
cp = plt.contourf(B,levels=np.arange(0,21),cmap=plt.cm.jet)
cp.cmap.set_over('k')
cp.set_clim(0,20)
#cbar2 = plt.colorbar(ticks=v)
cbar2 = plt.colorbar(cp)


####### SAVING FIGURE AS AN EPS #####


