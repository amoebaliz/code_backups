import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection

map = Basemap(resolution='c')
fig = plt.figure()
ax = fig.add_subplot(111)
map.drawcoastlines(linewidth=0)
#map.fillcontinents(color='gray',lake_color='gray')
polys = []

for polygon in map.landpolygons:
     polys.append(polygon.get_coords())                      

lc = PolyCollection(polys, edgecolor='black',
                    facecolor='green', closed=True)
ax.add_collection(lc)
plt.show()






