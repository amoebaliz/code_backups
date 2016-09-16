import scalab_plot as sp

rootdir='/t1/scratch/liz/tmpdir_VIP-LD.HCo13T/scalab/'
list=[]
for nc in 480, 240, 192, 144:
        list.append( rootdir + 'log.0.' + str(nc) )

myexp = {'machine':'triton','experiment':'VIP','compiler':'gfortran-4.4','timesteps':'1152','nbpoints':'0600x0340x050'}
vip = sp.scalability_test(list,myexp)


