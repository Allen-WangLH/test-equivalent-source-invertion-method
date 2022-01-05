# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 17:27:10 2020

@author: chens
"""

import numpy as np
from geoist.inversion import pfmodel_ts as its
import pandas as pd
import matplotlib.pyplot as plt

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

def addnoise(field, noise = 0.05):
    sigma = noise*(field.max()-field.min())
    n1 = np.random.normal(0,sigma,len(field))
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
    

    
    nyx = (2,2)
    xmin = 100
    xmax = 105
    ymin = 26
    ymax = 32
    
    cell_type = 'tesseroid'
    # source_volume = [xmin, xmax, ymin, ymax, -12500, -13500]
    source_volume = [xmin, xmax, ymin, ymax, -10000, -11000]
    margin = [1,1,1,1]

    
    smooth_components = ['dx','dy','dtt']
    # optimize_weights = ['obs','dx','dy','dtt']
    optimize_weights = ['obs','dtt']

    weights = {'obs':1,'refer':1,'dx':1,'dy':1,'dtt':1}
 
    model = its.InvModelTS(nyx=nyx,
                     source_volume=source_volume,
                          smooth_components=smooth_components,
                          weights=weights,
                          cell_type=cell_type,
                          margin = margin,
                          optimize_weights=optimize_weights,
                          data_dir='./input/section3.2_time')
    model.load_data(usecols=[0,1,3],header=1,delim_whitespace=False)
    refer_density = np.zeros(model.ns*model.nx*model.ny) #minimun model
    model.set_refer(refer_density)    
    model.gen_mesh()
    model.gen_kernel()
         # aabic[i]=model.calc_abic_quiet()
    model.abic_optimize()
    # model.do_linear_solve()
    model.plot_density()
  
    
    fieldinv = model.forward()
    model.plot_field(fieldinv,plot_station=False)
    
    res = model.orig_data['g']-fieldinv
    model.plot_field(res,plot_station=True)
    
    
    for i in range(40):
        filename = './output/section3.2_time/solution-mdtlu{}.grd'.format(i)
        print(filename,i*4,4*(i+1))
        srfgrd(model.solution[i*4:4*(i+1)], 
                              filename, nyx[0], nyx[1],xmin, xmax, ymin, ymax)

    #srfgrd(sol1, './solution-m21.grd', nyx[0], nyx[1],xmin, xmax, ymin, ymax)
    
    df1 = pd.DataFrame(columns=['lon','lat','field','inv','res'])
    df1['lon'] = model.orig_data['lon']
    df1['lat'] = model.orig_data['lat']
    df1['field'] = model.orig_data['g']
    df1['inv'] = fieldinv
    df1['res'] = res
    df1.to_csv('./output/section3.2_time/invfieldmlu.csv')    
    
   