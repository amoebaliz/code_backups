import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pyroms

# Grid Specific details
GRID = pyroms.grid.get_ROMS_hgrid('CORAL')
rho_mask_in = GRID.mask_rho
angle = GRID.angle_rho
x_vert = GRID.x_vert
y_vert = GRID.y_vert

rho_mask = np.zeros((rho_mask_in.shape[0]+2,rho_mask_in.shape[1]+2))
rho_mask[1:-1,1:-1] = rho_mask_in
rho_mask[0,1:-1] = rho_mask_in[0,:]
rho_mask[-1,1:-1] = rho_mask_in[-1,:]
rho_mask[1:-1,0] = rho_mask_in[:,0]
rho_mask[1:-1,-1] = rho_mask_in[:,-1]

xi_dif = np.diff(x_vert,axis=1)
eta_dif = np.diff(y_vert,axis=0)

# directionality of array = [[N],[E],[S],[W]]
grid_edge = np.array([[np.pi],[np.pi/2],[0],[3*np.pi/2]])
mask_store = np.empty([4,0])
dim_store = np.empty([4,0])
angle_store = np.empty([1,0])

for neta in range(rho_mask_in.shape[0]):
    for nxi in range(rho_mask_in.shape[1]):
        eta_2 = neta + 1
        xi_2 = nxi + 1
        # If point[neta,nxi] is land, check neighbors
        if rho_mask[eta_2,xi_2] == 0: 

           north_mask = rho_mask[eta_2+1,xi_2]
           east_mask = rho_mask[eta_2,xi_2+1]
           south_mask = rho_mask[eta_2-1,xi_2]
           west_mask = rho_mask[eta_2,xi_2-1]

           # If land_ point[neta,nxi] is NOT surrounded by land, store values
           if (north_mask + east_mask + south_mask + west_mask) < 4:
              # Add column to coastline ID array (4xN)
              #print mask_store.shape
              #print np.array([[north_mask],[east_mask],[south_mask],[west_mask]]).shape
              mask_store = np.concatenate((mask_store,np.array([[north_mask],[east_mask],[south_mask],[west_mask]])),axis=1)

              # Add column to cell rotation array (1xN)
              nangle = angle[neta,nxi]
              angle_store = np.concatenate((angle_store,np.array([[nangle]])),axis=1)
 
              # Add column to cell dimension array (4xN)
              north_dim = xi_dif[neta+1,nxi]
              east_dim = eta_dif[neta,nxi+1] 
              south_dim = xi_dif[neta,nxi] 
              west_dim = eta_dif[neta,nxi]
              dim_store = np.concatenate((dim_store,np.array([[north_dim],[east_dim],[south_dim],[west_dim]])),axis=1)

# Calculate land angle
edge_angle = np.tile(grid_edge,(1,angle_store.shape[1]))
grid_angle = np.tile(angle_store,(4,1))
land_angle = edge_angle+grid_angle 

# Compute K value 
K = mask_store*dim_store

# Wind stress components
N = 500

wind_angle = np.linspace(0,2*np.pi,N)
x_stress = np.cos(wind_angle)
y_stress = np.sin(wind_angle)

#Initialize Storage matrix
store_pot = np.zeros(len(wind_angle)) 

# Iterate over 0-360 wind angles
for nt in range(len(x_stress)):

    # Calculate cos(angle dif: wind-land)
    cos_angle = np.cos(wind_angle[nt] - land_angle)

    #Compute summation matrix
    wind_stress_pot = K*cos_angle

    # Only consider upwelling
    upwell = np.ma.masked_where(wind_stress_pot<0,wind_stress_pot)
    upwelling_pot = np.sum(upwell)
    store_pot[nt] = upwelling_pot

# Plot in polar coord 
norm_pot = store_pot/np.max(store_pot)  
ax = plt.subplot(111, projection='polar')
for nt in range(len(norm_pot)):
    ax.plot([0,wind_angle[nt]],[0,norm_pot[nt]],'r')
ax.set_rmax(1.5)
ax.grid(True)
#ax.set_yticklabels([1])
ax.set_title("CT_ROMS coastal upwelling potential", va='bottom')
plt.show()    

