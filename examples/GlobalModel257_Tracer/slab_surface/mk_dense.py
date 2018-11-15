#!/usr/bin/python
import os, string,math,sys

age = int(sys.argv[1])
sz_dir='/net/beno/data1/jiashun/weakzone_grd/slab_surface/'
r2d=180.0/math.pi
d2r=math.pi/180.0
r=6371.0
input_file_name="%sreconstructed_%dMa.xy" % (sz_dir, age)
input_file=open(input_file_name,'r')
output_file_name="dense_reconstructed_%dMa.xy" % age 
output_file=open(output_file_name,'w')
min_deg=0.1
max_deg=50
#record the last point of a trench
temp_lat=0.0
temp_lon=0.0
#record the number of lines
num_line=0
#determine where it is the first line that has '>' in a bunch of lines with '>'
firstline=0
while 1:
        line=input_file.readline()
        if(line):
            c1 = line[0]
            if c1 == '>':
                m=0 
		#only for fixing the bug of deleting the last point of a trench
		if firstline==0:
		    firstline=1
		else:
		    firstline=2
		#only for fixing the bug of deleting the last point of a trench
		if num_line!=0 and firstline==1:
		    output_file.write("%f %f\n" % (temp_lon,temp_lat))
        #       print 'make new line'
                #NEW.write("%s\n" % gmt_char)
                output_file.write("%s" % line)
            if c1 != '#' and c1 != '>':
		firstline=0
                m+=1
                lon, lat =line.split()
                flat=float(lat)
                flon=float(lon)
                if flon < 0.0:
                    flon=flon+360.0
                if flon > 360.0:
                    flon=flon-360.0
                if m>1:
                    dist=0.0
                    if (flon-flon_1)<-200.0:
                        flon=flon+360
                    if (flon-flon_1)>200:
                        flon=flon-360
                    if (flon_1 > 0.1 and flon > 0.1) or (flon_1 < 359.9 and flon < 359.9) :
                        dist=math.sqrt((flat_1-flat)**2 + (flon_1-flon)**2)
                    if dist > min_deg:
#and dist<max_deg:
                        num_new=int(dist/min_deg)
                        dlon=(flon-flon_1)/num_new
                        dlat=(flat-flat_1)/num_new
                        nn=0
                        flon_new=flon_1
                        flat_new=flat_1
                        while nn<num_new:
                            output_file.write("%f %f\n" %(flon_new,flat_new))
                            flon_new=flon_new+dlon
                            flat_new=flat_new+dlat
                            nn+=1
                    else:
			#here is a bug that has been commented out
                        #output_file.write("%f %f\n" %(flon,flat))
			output_file.write("%f %f\n" %(flon_1,flat_1))
                flat_1=flat
                flon_1=flon
		#only for fixing the bug of deleting the last point of a trench
		temp_lat=flat
		temp_lon=flon
		#only for fixing the bug of deleting the last point of a trench
		num_line=num_line+1
        else:
                break
input_file.close()
output_file.close()
