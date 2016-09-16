#!/usr/bin/env python

import scalab_plot as sp

rootdir='./logs_ccs1fish_ys/'
list=[]
for nc in 512, 256, 128, 96, 64, 32, 16:
        list.append( rootdir + 'log.' + str(nc) )

myexp = {'machine':'yellowstone','experiment':'CCS1-Fish','compiler':'ifort','timesteps':'1152','nbpoints':'0182x0482x050'}
ccs1 = sp.scalability_test(list,myexp)


