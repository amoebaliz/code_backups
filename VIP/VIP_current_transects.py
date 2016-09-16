import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp
import pyroms

def plot_domain(GRID):
    bathym = np.ma.masked_where(GRID.hgrid.mask_rho == 0, GRID.vgrid.h)
    plt.figure(figsize=(16,12))
    plt.pcolor(bathym,vmin=0,vmax=3600)
    plt.ylim (0,bathym.shape[0])         
    plt.xlim (0,bathym.shape[1])
    # plot transect
    for ntrans in range(4):
        a = 2*ntrans
        b = 2*ntrans+2
        plt.plot(tsect_xi[ntrans]*np.ones(2),tsect_eta[a:b],'k',lw=3) 
        plt.plot(tsect_xi[ntrans],tsect_eta[a]-2,'o',ms=20,mfc='w',mew=1,alpha=0.9)
        plt.text(tsect_xi[ntrans],tsect_eta[a]-2,panel_id[ntrans],fontsize=13,ha = 'center',va ='center')
#        plt.text(tsect_xi[ntrans] + 2,tsect_eta[a] + map_offsets[ntrans],panel_id[ntrans],fontsize=13)
def plot_transects(GRID,nvars,ntransects):

    fig=plt.figure(figsize=(14,6))
 
    # setting relative widths of transect plots
    dy_vals = GRID.hgrid.dy
    widths = np.zeros(ntransects)
    for ntrans in range(ntransects):  
        a = 2*ntrans
        widths[ntrans] = np.sum(dy_vals[tsect_eta[a]:tsect_eta[a+1]+1,tsect_xi[ntrans]])
    gs = gsp.GridSpec(nvars,ntransects,width_ratios = widths)
    gs.update(bottom = 0.2, left = 0.1, right = 0.88, top = 0.87, wspace=0.1)
    cbar_ax_1 = fig.add_axes([0.91,0.59,0.015,0.2])
    cbar_ax_2 = fig.add_axes([0.91,0.22,0.015,0.2])
    cmap = plt.cm.bwr

    for nvar in range(nvars):
        ncfile = data_dir + data_files[nvar]
        fid = nc.Dataset(ncfile)
        u_var = np.squeeze(fid.variables['u'][:])

        if nvar == 1:
           bias_frac = ((365.*4+1)/(365*4))
           u_var = np.sqrt(bias_frac*u_var) 

        for ntrans in range(4):
            a = 2*ntrans
            col = tsect_xi[ntrans]
            min_row = tsect_eta[a]
            max_row = tsect_eta[a+1]
            rows = np.arange(min_row,max_row+1)
            X = np.tile(rows,(50,1))
            u_slice, u_Z, u_X, u_Y = pyroms.tools.islice(u_var,col,GRID,'u')

            if ntrans == 0:
               axs = fig.add_subplot(gs[nvar,ntrans])
               if nvar == 0:
                  axs.set_ylabel('Depth (m)', x= -20, y= -.1,fontsize=18)
               axs.set_yticks([max_axis_dep+100,(max_axis_dep+100)/2,0])
               axs.text(min_row+2,15,Panel_ID[nvar],fontsize=14)
            else:
               axs = fig.add_subplot(gs[nvar,ntrans],sharey=axs)
               plt.setp(axs.get_yticklabels(),visible=False)

            C = axs.pcolor(X,u_Z[:,rows],u_slice[:,rows],vmin=vmins[nvar],vmax=vmaxs[nvar],cmap=cmaps[nvar])
            axs.set_ylim((max_axis_dep,0))
            axs.tick_params(axis = 'y',labelsize=13)
            axs.set_xlim(min_row, max_row)
            axs.get_xaxis().set_ticks([]) 
            axs.set_axis_bgcolor((0.6,0.6,0.6))
            axs.yaxis.tick_left()
            #axs.plot([tsect_eta[a],tsect_eta[a+1]],max_axis_dep*np.ones(2),'-k',lw=2,clip_on=False,zorder=99)
            #axs.plot(tsect_eta[a],max_axis_dep,'o',ms=20,mfc='w',mew=1,clip_on=False,zorder=100)
            #axs.text(tsect_eta[a],max_axis_dep, panel_id[ntrans],fontsize=13,ha = 'center',va ='center',zorder=101)

            bar_depth = -800
            bar_end   = 192
            if nvar == 0 and ntrans == 3:
               # Add 10 km scalebar
               bar_len = (10000/widths[ntrans])*(max_row-min_row+1)
               axs.plot([bar_end,bar_end-bar_len],[bar_depth,bar_depth],'-k',lw=2)
               axs.text(bar_end-3*bar_len/4,bar_depth+20,'10 km',fontsize=13)
 
        if nvar == 0:
           fig.colorbar(C,cax=cbar_ax_1,ticks = tick_vals[nvar])
           cbar_ax_1.set_title('Velocity \n (m/s)',y=1.1)
           cbar_ax_1.tick_params(labelsize=13)
        #else:
        #   fig.colorbar(C,cax=cbar_ax_2,ticks = tick_vals[nvar])
        #   cbar_ax_2.set_title('Standard \n Deviation \n (m/s)',y=1.1)
        #   cbar_ax_2.tick_params(labelsize=13)
 
# ------------------------------------------------------ #
VIP = pyroms.grid.get_ROMS_grid('VIP')

data_dir = '/t3/workdir/liz/MODELS/VIP/PostProc/'
data_files = ['VIP_U_JJAS_1999.nc']

#tsect_eta = np.array([70,148,95,152,144,178,140,200])
#tsect_xi = np.array([192,240,280,314])

tsect_eta = np.array([65,149,102,152,145,178,146,199])
tsect_xi = np.array([180,240,280,318])
map_offsets = [-1,-2,-3,-3]
vmins = [-0.5,0]
vmaxs = [0.5,.6]
cmaps = ['PuOr','Purples']
panel_id = ['i','ii','iii','iv']
Panel_ID = ['A','B']
tick_vals = np.array([[-.6,0,.6],[0,.3,.6]])
max_axis_dep = -900

#plot_domain(VIP)
plot_transects(VIP,1,4) 
plt.savefig('a2.png')
plt.show()
