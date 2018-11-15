#!/usr/bin/env python

import sys

age=int(sys.argv[1])
print age
nx=257
ny=257
nz=65

for time in range(age,age+100,100):
    for cap in range(12):
	if cap<10:
            fname="global.comp0%d.%d" % (cap,time)
	else:
            fname="global.comp%d.%d" % (cap,time)
        infile=open(fname,"r");
        lines=infile.readlines();
        infile.close()
        
	if cap<10:
            outname="global.0%d.%d.txt" % (cap,time)
	else:
            outname="global.%d.%d.txt" % (cap,time)
        outfile=open(outname,"w")
        
        for j in range(ny):
            for i in range(nx):
                for k in range(nz):
        	    if i%2==0 and j%2==0:
    		    #if i%3==0:
    		    #if i%4==0 and j%4==0:
        		outfile.write(lines[j*nx*nz+i*nz+k+1])
        
        outfile.close()
