#!/usr/bin/env python

import os,sys
import math
import time
import Core_GMT
import numpy as np

timestep=int(sys.argv[1])

slab_surface_grd="weakzone_tracer.grd"
fname="topology_subduction_boundaries_%.2fMa.xy" % timestep
fname="dense_reconstructed_%dMa.xy" % timestep

# get parameters for the grd file
dict={}
dict['grid']=slab_surface_grd
Core_GMT.grdinfo(dict)
long_min = float(dict['west'])
long_max = float(dict['east'])
lat_min = float(dict['south'])
lat_max = float(dict['north'])
dx = float(dict['dx'])
#dx = float("%.2f" % dx)
dy = float(dict['dy'])
#dy = float("%.2f" % dy)

trench_clean_f = fname + ".clean"
cmd="awk '!seen[$0]++' " + fname + " > " + trench_clean_f
print cmd
#os.system(cmd)

#trench_f=open(trench_clean_f,'r')
trench_f=open(fname,'r')
lines_rawdata=trench_f.readlines()
trench_f.close()
lines=[]
for line in lines_rawdata:
    if line[0]!='>' or line.find('name')!=-1:
        lines.append(line)

i=0
trenches=[]
trench_names=[]
length=[]
while i<len(lines):
    if lines[i][0]=='>':
        if i==len(lines)-1:
            break
        trench_pts=[]
        trench_pts.append(lines[i])
        trench_name=lines[i].split("#")[1].split(":")[1].strip()
        trench_name=trench_name.replace(" ","_").replace("(","").replace(")","").replace("/","")
        trench_name="trench."+trench_name+".xy"
        #print trench_name
        i=i+1
        while i<len(lines) and lines[i][0]!='>':
            trench_pts.append(lines[i])
            i=i+1
        if trench_name not in trench_names:
            trench_names.append(trench_name)
        else:
            trench_name=trench_name.rstrip(".xy")+time.ctime().split()[3].replace(":","-")+".xy"
            trench_names.append(trench_name)
        trenches.append(trench_pts)
        length.append(len(trench_pts))

indexes=np.argsort(length)[::-1]
for index in indexes:
    trench_pts=trenches[index]
    trench_name=trench_names[index]
    print trench_name

    # output the pts along the trench if the length > 1
    if len(trench_pts)>1:
        # the output file for the trench
        trench_f=open(trench_name,'w')
        for pts in trench_pts:
            trench_f.write("%s" % pts)
        trench_f.close()

        # compute the continental side dist through the 
        # assimilation depth function
        cmd="./mk_assim %s" % trench_name
        print cmd
        os.system(cmd)

        # grid the output data file
        res=0.1
        cmd="blockmean asm_depth.xy -I%(res)f -C -R%(long_min)f/%(long_max)f/%(lat_min)f/%(lat_max)f > asm_depth.blockmean.txt" % vars()
        print cmd
        os.system(cmd)
        cmd="surface asm_depth.blockmean.txt -I%(res)f -R%(long_min)f/%(long_max)f/%(lat_min)f/%(lat_max)f -T0.1 -Gasm_depth.grd" % vars()
        print cmd
        os.system(cmd)
        # clip the values > 0
        cmd="grdclip asm_depth.grd -Gasm_depth.clip.grd -Sa0.05/NAN -V"
        print cmd
        os.system(cmd)
        # resample the asm_depth grd file
        cmd="grdsample asm_depth.clip.grd -I%(dx)f/%(dy)f -Gasm_depth.clip.resample.grd" % vars()
        print cmd
        os.system(cmd)

        cmd="awk '{if($3<0) print $1,$2,-$3; else print $1,$2,$3}' asm_depth.xy | blockmean -I%(res)f -C -R%(long_min)f/%(long_max)f/%(lat_min)f/%(lat_max)f > asm_depth.blockmean.txt" % vars()
        print cmd
        os.system(cmd)
        cmd="surface asm_depth.blockmean.txt -I%(res)f -R%(long_min)f/%(long_max)f/%(lat_min)f/%(lat_max)f -T0.1 -Gasm_depth2.grd" % vars()
        print cmd
        os.system(cmd)
        # clip the values > 0
        cmd="grdclip asm_depth2.grd -Gasm_depth2.clip.grd -Sa0.8/NAN -V"
        print cmd
        os.system(cmd)
        # resample the asm_depth grd file
        cmd="grdsample asm_depth2.clip.grd -I%(dx)f/%(dy)f -Gasm_depth2.clip.resample.grd" % vars()
        print cmd
        os.system(cmd)

        # mask the surface slab grd file
        trench_grd_name=trench_name.replace(".xy",".grd")
        cmd="grdmath %(slab_surface_grd)s asm_depth.clip.resample.grd OR asm_depth2.clip.resample.grd OR = %(trench_grd_name)s" % vars()
        print cmd
        os.system(cmd)

        # plot out grd
        cmd="makecpt -Crainbow -T-10/200/10 -Z > depth.cpt"
        print cmd
        os.system(cmd)
        trench_ps_name=trench_name.replace(".xy",".ps")
        cmd="grdimage %(trench_grd_name)s -Cdepth.cpt -JN180/5.0i -B30 -X3 -Y15 -K -P > %(trench_ps_name)s" % vars()
        print cmd
        os.system(cmd)
        cmd="psxy %(fname)s -J -m -R -fg -Wwhite -K -O -P >> %(trench_ps_name)s" % vars()
        print cmd
        os.system(cmd)
        cmd="psxy %(trench_name)s -J -m -R -fg -Worange -K -O -P >> %(trench_ps_name)s" % vars()
        print cmd
        os.system(cmd)
        cmd="psscale -Cdepth.cpt -D3/-0.5/8/0.5h -Ba50g25 -X0.0 -Y-3.0 -O -K -P >> %(trench_ps_name)s" % vars()
        print cmd
        os.system(cmd)

cmd="rm trench*.xy asm_depth*"
print cmd
os.system(cmd)
