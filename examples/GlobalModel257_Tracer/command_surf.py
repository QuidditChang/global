#!/usr/bin/env python

import os

timestep=0
dest="/u/sciteam/hu1/scratch/GlobalModel257_Tracer/DATA"
cmd="cp %s/*/global.coord.* %s" % (dest, dest)
os.system(cmd)

cmd="cp %s/*/global.surf.*.%d %s" % (dest, timestep, dest)
os.system(cmd)

cmd="~/AssimDepth_CitcomS-3.0.3/visual/batchsurf.py localhost DATA global %d 257 257 129 12 8 8 4" % timestep
print cmd
os.system(cmd)

cmd="rm %s/global.coord.*" % dest
os.system(cmd)

cmd="rm %s/global.surf.*.%d" % (dest, timestep)
os.system(cmd)
