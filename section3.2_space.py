# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 17:27:10 2020

@author: chens
"""

import numpy as np
from geoist.inversion import pfmodel_ts as its
import pandas as pd
import matplotlib.pyplot as plt
import openpyxl as xl

def checkboard(density, nx, ny, size):
    solution = np.ones([nx,ny])*density
    for i in range(nx):
        inum = int(i/size)
        for j in range(ny):
            jnum = int(j/size)
            if (inum % 2) == 0:
                if (jnum % 2 ) == 0:
                    solution[i,j] = solution[i,j]*-1.0
            else:
                if (jnum % 2 ) != 0:
                    solution[i,j] = solution[i,j]*-1.0
                
    return solution

def addnoise(field, noise):

    sigma = noise*(field2.max()-field2.min())
#    sigma =0.005
    n1 = np.random.normal(0,sigma,len(field2))
    return field + n1

def srfgrd(solution, filename, nx, ny, xmin, xmax, ymin, ymax):
    from geoist.pfm.grdio import grddata
    g1out = grddata()
    g1out.cols = ny
    g1out.rows = nx
    g1out.xmin = xmin
    g1out.ymin = ymin
    g1out.xmax = xmax
    g1out.ymax = ymax    
    g1out.data0 = solution.reshape(nx,ny)
    #g1out.data0 = np.transpose(solution.reshape(nx,ny))
    g1out.export_surfer(filename, False, 'ascii')


if __name__ == '__main__':
    

    
    nyx = (12,10)
    xmin = 100
    xmax = 105
    ymin = 26
    ymax = 32
    

    cell_type = 'tesseroid'
    source_volume = [xmin, xmax, ymin, ymax, -10000, -11000]
    margin = [1,1,1,1]

    smooth_components = ['dx','dy','dxy']
    optimize_weights = ['obs','dx','dy','dxy']
    weights = {'obs':1,'dx':1,'dy':1,'dxy':1} #,'dt':1}

# 0.158166352
# 7.59E-08
# 0.153465488
# 0.540261426 
    
    model = its.InvModelTS(nyx=nyx,
                     source_volume=source_volume,
                          smooth_components=smooth_components,
                          weights=weights,
                          cell_type=cell_type,
                          margin = margin,
                          optimize_weights=optimize_weights,
                          data_dir='./input/section3.2')
    model.load_data(usecols=[0,1,3],header=1,delim_whitespace=False)
    #model.deg2xy()
    model.plot_field(surveys=[0],plot_station=True)
    model.gen_mesh()
    model.gen_kernel()
    sol1 = checkboard(0.001, nyx[0], nyx[1], 1)
    # soln = np.concatenate([sol1, sol1, sol1, sol1])
    model.plot_density(sol1)

    
    field2  = model.forward(sol1)
    

    # field2n = addnoise(field2, 0.03)
    field3 = xl.load_workbook('modnoise05.xlsx')
    field4 = field3.get_sheet_by_name('Sheet1')
    field2n=[item.value for item in list(field4.columns)[0]]
    field2n=np.array(field2n)
    
    model.orig_data['g'] = field2n
    # model.orig_data['g'] = field2
    model.plot_field(field2n,plot_station=False)
    # model.plot_field(field2,plot_station=False)
    
    model.abic_optimize()
    # model.do_linear_solve()
    model.plot_density()
    
    filedinv = model.forward()
    model.plot_field(filedinv,plot_station=False)
    
    res = model.orig_data['g']-filedinv
    model.plot_field(res,plot_station=True)
    

    
    srfgrd(model.solution, './output/section3.2_space/solution-inv.grd', nyx[0], nyx[1],xmin, xmax, ymin, ymax)
    srfgrd(sol1, './output/section3.2_space/solution-model.grd', nyx[0], nyx[1],xmin, xmax, ymin, ymax)
    
    df1 = pd.DataFrame(columns=['lon','lat','mod','modwithnoise','noise','modfy','res'])
    df1['lon'] = model.orig_data['lon']
    df1['lat'] = model.orig_data['lat']
    df1['mod'] = field2
    df1['modwithnoise'] = field2n
    df1['noise'] = field2n-field2
    df1['modfy'] = filedinv
    df1['res'] = filedinv-field2n
    df1.to_csv('./output/section3.2_space/mod2-inv.txt') 
    
   