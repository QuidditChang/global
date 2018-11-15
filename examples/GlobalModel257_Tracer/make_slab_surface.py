#!/usr/bin/env python

import os,sys
import math
import numpy as np

timestep=int(sys.argv[1])

for cap in range(12):
    if cap<10:
        cmd="awk '{if(NR>1 && $3>0.965) print $1,$2,$3,$9}' combined_files/global.comp0%d.%d > all.coor.comp.txt" % (cap,timestep)
        print cmd
        os.system(cmd)
    else:
        cmd="awk '{if(NR>1 && $3>0.965) print $1,$2,$3,$9}' combined_files/global.comp%d.%d > all.coor.comp.txt" % (cap,timestep)
        print cmd
        os.system(cmd)
    
    nx=257
    ny=257
    nz=18
    R0=6371.0
    max_dist=30.0
    r2d=180.0/math.pi
    data=np.loadtxt('all.coor.comp.txt')
    data[:,0]=90.0-data[:,0]*r2d
    data[:,1]=data[:,1]*r2d
    data[:,2]=(1-data[:,2])*R0
    data=data[::-1,:]
    store_data=[]
    for i in range(nx):
        for j in range(ny):
            depths=[]
            for k in range(nz):
                n=k+j*nz+i*ny*nz
                if data[n,3]==1:
                    depths.append(data[n,2])
                if data[n,2]>200.0:
                    break
            if len(depths)==0:
                store_data.append([data[n,0],data[n,1],999])
            elif len(depths)==1:
                store_data.append([data[n,0],data[n,1],data[n,2]])
            else:
                for ii in range(len(depths)-1):
                    if depths[ii+1]-depths[ii]>max_dist:
                        break
                store_data.append([data[n,0],data[n,1],np.mean(depths[:i+1])])
    
    store_data=np.array(store_data)
    if cap<10:
        np.savetxt('weakzone_tracer.0%d.txt' % cap,store_data)
    else:
        np.savetxt('weakzone_tracer.%d.txt' % cap,store_data)

cmd="cat weakzone_tracer.[01][0-9].txt > weakzone_tracer.txt"
print cmd
os.system(cmd)
# convert to grd
res=0.1
cmd="blockmean weakzone_tracer.txt -I%(res)f -di999 -R0/360/-89.9/89.9 -: > weakzone_tracer.blockmean.txt" % vars()
print cmd
os.system(cmd)
cmd="surface weakzone_tracer.blockmean.txt -I%(res)f -R0/360/-89.9/89.9 -T0.1 -di999 -Gweakzone_tracer.grd -:" % vars()
print cmd
os.system(cmd)
# plot out grd
cmd="makecpt -Crainbow -T-10/200/10 -Z > depth.cpt"
print cmd
os.system(cmd)
cmd="grdimage weakzone_tracer.grd -Cdepth.cpt -JN180/5.0i -B30 -X3 -Y15 -K -P > weakzone.ps"
print cmd
os.system(cmd)
cmd="psscale -Cdepth.cpt -D3/-0.5/8/0.5h -Ba50g25 -X0.0 -Y-3.0 -O -K -P >> weakzone.ps"
print cmd
os.system(cmd) 
                

cmd="rm all.comp.txt all.coor.comp.txt weakzone_tracer.blockmean.txt depth.cpt"
print cmd
os.system(cmd)
