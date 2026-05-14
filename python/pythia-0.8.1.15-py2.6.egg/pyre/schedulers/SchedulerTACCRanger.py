#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                      California Institute of Technology
#                        (C) 2008  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from SchedulerSGE import SchedulerSGE
import os, sys


class SchedulerTACCRanger(SchedulerSGE):
    
    
    name = "tacc-ranger"
    

    import pyre.inventory as pyre
    
#    command      = pyre.str("command", default="qsub")
    command      = pyre.str("command", default="sbatch")
    tpn          = pyre.int("tpn", default=16,
                            validator=pyre.choice([1,2,4,8,12,15,16,32,34]))
    tpn.meta['tip'] = 'Task per node'
    qsubOptions  = pyre.list("qsub-options")

    
    def schedule(self, job):
        from math import ceil
        # round up to multiple of 16
        nodes = ceil(job.nodes / float(self.tpn))

	# Lijun add
	self.nodes = nodes
	# end of add

        self.cores = int(nodes * self.tpn)

        SchedulerSGE.schedule(self, job)
        
# end of file 
