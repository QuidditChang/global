#!/usr/bin/env python

import os,sys,math
import numpy as np
from scipy import interpolate
import matplotlib.mlab as mlab

#=========================================================================================================================
def my_interp(all_data,xmin,xmax,ymin,ymax,res):
    print xmin, ymin
    f1=np.logical_and(all_data[:,0]>xmin-20*res,all_data[:,0]<xmax+20*res)
    f2=np.logical_and(all_data[:,1]>ymin-20*res,all_data[:,1]<ymax+20*res)
    fl=np.logical_and(f1,f2)
    x = np.arange(xmin+res/2, xmax, res)
    y = np.arange(ymin+res/2, ymax, res)
    xx, yy = np.meshgrid(x,y)
    zz = mlab.griddata(all_data[fl,0], all_data[fl,1], all_data[fl,2], xx, yy, interp='linear')
    data=np.array([xx.flatten(),yy.flatten(),zz.flatten()]).T
    np.savetxt('interpolated.%f.%f.txt' % (xmin,ymin),data)

    return
#=========================================================================================================================
def my_interp2(all_data,xmin,xmax,ymin,ymax,res):
    f1=np.logical_and(all_data[:,0]>xmin-10*res,all_data[:,0]<xmax+10*res)
    f2=np.logical_and(all_data[:,1]>ymin-10*res,all_data[:,1]<ymax+10*res)
    fl=np.logical_and(f1,f2)
    f = interpolate.interp2d(all_data[fl,0], all_data[fl,1], all_data[fl,2], kind='linear')
    x = np.arange(xmin+res/2, xmax, res)
    y = np.arange(ymin+res/2, ymax, res)
    xx, yy = np.meshgrid(x,y)
    zz = f(x,y)
    data=np.array([xx.flatten(),yy.flatten(),zz.flatten()]).T
    np.savetxt('interpolated.%f.%f.txt' % (xmin,ymin),data)
    
    return
#=========================================================================================================================
def main():
    timestep=int(sys.argv[1])
    
    model_dir="/home1/02505/hujs/scratch/GlobalModel257_Tracer_plume_80Ma"
    for cap in range(12):
        if cap<10:
            cmd="awk '{if(NR>1 && $3>0.965) print $1,$2,$3,$9}' %s/combined_files/global.comp0%d.%d > all.coor.comp.txt" % (model_dir,cap,timestep)
            print cmd
            os.system(cmd)
        else:
            cmd="awk '{if(NR>1 && $3>0.965) print $1,$2,$3,$9}' %s/combined_files/global.comp%d.%d > all.coor.comp.txt" % (model_dir,cap,timestep)
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
        for i in range(2,nx-2):
            for j in range(2,ny-2):
                depths=[]
                for k in range(nz):
                    n=k+j*nz+i*ny*nz
                    if data[n,3]==1:
                        depths.append(data[n,2])
                    if data[n,2]>200.0:
                        break
                if len(depths)==0:
                    store_data.append([data[n,0],data[n,1],float('nan')])
                elif len(depths)==1:
                    store_data.append([data[n,0],data[n,1],data[n,2]])
                else:
                    for ii in range(len(depths)-1):
                        if depths[ii+1]-depths[ii]>max_dist:
                            break
                    store_data.append([data[n,0],data[n,1],np.mean(depths[:i+1])])
        
        store_data=np.array(store_data)
        if cap==0:
            all_data=store_data
        else:
            all_data=np.append(all_data,store_data,axis=0)
        if cap<10:
            np.savetxt('weakzone_tracer.0%d.txt' % cap,store_data)
        else:
            np.savetxt('weakzone_tracer.%d.txt' % cap,store_data)

    res=0.1
    incr=60
    for x in range(-90,90,incr):
        for y in range(0,360,incr):
            my_interp(all_data,x,x+incr,y,y+incr,res)

    cmd="cat interpolated.*000000.txt > weakzone_tracer.txt"
    print cmd
    os.system(cmd)
    # convert to grd
    res=0.125
    cmd="blockmean weakzone_tracer.txt -I%(res)f -C -R0/360/-89.9/89.9 -: > weakzone_tracer.blockmean.txt" % vars()
    print cmd
    os.system(cmd)
    cmd="surface weakzone_tracer.blockmean.txt -I%(res)f -R0/360/-89.9/89.9 -T0.5 -Gweakzone_tracer.%(timestep)d.grd -:" % vars()
    print cmd
    os.system(cmd)
    cmd="xyz2grd weakzone_tracer.blockmean.txt -I%(res)f -R0/360/-89.9/89.9 -Gweakzone_tracer.mask.%(timestep)d.grd -:" % vars()
    print cmd
    os.system(cmd)
    cmd="grdmath weakzone_tracer.%(timestep)d.grd weakzone_tracer.mask.%(timestep)d.grd OR = weakzone_tracer.%(timestep)d.grd" % vars()
    print cmd
    os.system(cmd)
    # plot out grd
    cmd="makecpt -Crainbow -T-10/200/10 -Z > depth.cpt"
    print cmd
    os.system(cmd)
    cmd="grdimage weakzone_tracer.%(timestep)d.grd -Cdepth.cpt -JN180/5.0i -B30 -X3 -Y15 -K -P > weakzone.%(timestep)d.ps" % vars()
    print cmd
    os.system(cmd)
    cmd="psscale -Cdepth.cpt -D3/-0.5/8/0.5h -Ba50g25 -X0.0 -Y-3.0 -O -K -P >> weakzone.%(timestep)d.ps" % vars()
    print cmd
    os.system(cmd) 
                    
    
    cmd="rm all.comp.txt all.coor.comp.txt interpolated.*.txt weakzone_tracer*.txt depth.cpt"
    print cmd
    os.system(cmd)
#=========================================================================================================================
if __name__=='__main__':
    main()
