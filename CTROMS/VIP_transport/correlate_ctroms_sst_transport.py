import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# POINTS ALONG VIP TO COMPARE...
itime = range(732,18262+1)
ilons = range(366,426+1)
jlats = range(450,500+1)

# LOAD TRANSPORT TIME SERIES
fid = nc.Dataset('/t3/workdir/liz/scripts/CT_roms/ctroms_u_transport_col395.nc','r')
pos_transport = fid.variables['pos_transport'][itime]
neg_transport = fid.variables['neg_transport'][itime]
Ttrans = pos_transport + neg_transport
# GET TEMPERATURE RECORDS FOR DIFFERENT POINTS
fid2 = nc.Dataset('/t3/workdir/liz/scripts/CT_roms/ctroms_sst.nc')
mask_val = -7
def get_temp(j,i):
    sst = np.array(np.squeeze(fid2.variables['temp'][itime,:,jlats[j],ilons[i]]))
    return sst

def cor_ssvals(j,i):
    sst_loc = get_temp(j,i)
    cor_mat = np.corrcoef(Ttrans,sst_loc)
    T = cor_mat[0,1]
    if np.isnan(T):
       T = mask_val
    #print T
    return T

cor_vals = np.zeros((len(jlats),len(ilons)))
for i in range(len(ilons)):
    for j in range(len(jlats)):
         cor_vals[j,i] = cor_ssvals(j,i)

plot_vals = np.ma.masked_equal(cor_vals,mask_val)

plt.figure()
plt.pcolor(plot_vals, vmin=-1, vmax=1)
plt.title('SST correlation w/ surface transport')
plt.plot([395-ilons[0],395-ilons[0]],[469-jlats[0],473-jlats[0]],'-k',linewidth=2.0)

plt.xlim(0,len(ilons))
plt.ylim(0,len(jlats))
plt.colorbar()
plt.show()
