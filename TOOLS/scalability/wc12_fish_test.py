import scalab_plot as sp

rootdir='/t1/raphael/tests_ROMS/WC12_SAN/results/scalab_tests/'
list=[]
for nc in 256, 192, 128, 96, 64, 32, 16:
        list.append( rootdir + 'ncores_' + str(nc) + '/log.wc12.' + str(nc) )

myexp = {'machine':'triton','experiment':'WC12-Fish','compiler':'gfortran-4.4','timesteps':'2880','nbpoints':'0184x0179x042'}
wc12 = sp.scalability_test(list,myexp)


