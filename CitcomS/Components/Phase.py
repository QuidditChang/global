#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#<LicenseText>
#
# CitcomS.py by Eh Tan, Eun-seo Choi, and Pururav Thoutireddy.
# Copyright (C) 2002-2005, California Institute of Technology.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#</LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from CitcomComponent import CitcomComponent

class Phase(CitcomComponent):


    def __init__(self, name="phase", facility="phase"):
        CitcomComponent.__init__(self, name, facility)
        return



    def setProperties(self, stream):
        from CitcomSLib import Phase_set_properties
        Phase_set_properties(self.all_variables, self.inventory, stream)
        return


    class Inventory(CitcomComponent.Inventory):


        import pyre.inventory


        phase_depth = pyre.inventory.array(
            "phase_depth", default=[410.0/6371.0, 520.0/6371.0, 660.0/6371.0])
        phase_delta_rho = pyre.inventory.array(
            "phase_delta_rho", default=[0.0, 0.0, 0.0])
        phase_delta_s = pyre.inventory.array(
            "phase_delta_s", default=[0.0, 0.0, 0.0])
        # Deprecated compatibility input. Runtime derives this from delta rho.
        phase_Ra = pyre.inventory.array(
            "phase_Ra", default=[0.0, 0.0, 0.0])
        phase_width = pyre.inventory.array(
            "phase_width", default=[0.0058, 0.0058, 0.0058])
        phase_clapeyron = pyre.inventory.array(
            "phase_clapeyron", default=[0.0235, 0.0, -0.0235])
        phase_transT = pyre.inventory.array(
            "phase_transT", default=[0.78, 0.78, 0.78])


# version
__id__ = "$Id: Phase.py 4957 2006-10-12 14:48:43Z leif $"

# End of file
