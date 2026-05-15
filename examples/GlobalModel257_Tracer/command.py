#!/usr/bin/env python

import sys, string, os, math

cmd="mkdir temp"
print "cmd=",cmd
#os.system(cmd)

i=660
while i<=660:
	cmd1="autocombine.py localhost pid*.cfg %d" %(i*1)
	print "cmd1=",cmd1
       	os.system(cmd1)

	cmd="tail -n +2 -q ~/scratch/GlobalModel257_Tracer_plume_80Ma_output/DATA/*/global.tracer*.%d > tracer_combined.%d" % (i,i*1)
        print "cmd=",cmd
        #os.system(cmd)

        #cmd="cp sub.cap00.%d tracer_combined.%d ~/work/2D_Lab/" % (i*1,i*1)
       # cmd1="cp sub.cap00.%d /u/sciteam/zhou2/2D_plot/sub.cap00.%d" % (i*1,i*1)
 #       cmd2="cp tracer_combined.%d ~/work/2D_plot/tracer_combined.%d" % (i*1,i*1)
       # print "cmd1=",cmd1
  #      print "cmd2=",cmd2
       # os.system(cmd1)
   #     os.system(cmd2)

	i+=200

cmd="cp DATA/0/sub.time /u/sciteam/zhou2/2D_plot"  
print "cmd=",cmd
#os.system(cmd)
cmd="cp pid* /u/sciteam/zhou2/2D_plot"
print "cmd=",cmd
#os.system(cmd)
