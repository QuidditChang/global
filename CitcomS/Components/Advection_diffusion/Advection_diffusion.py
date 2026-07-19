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

from CitcomS.Components.CitcomComponent import CitcomComponent


class Advection_diffusion(CitcomComponent):


    def __init__(self, name, facility):
        CitcomComponent.__init__(self, name, facility)
        return



    def setProperties(self, stream):
        from CitcomSLib import Advection_diffusion_set_properties
        Advection_diffusion_set_properties(self.all_variables, self.inventory, stream)
        return



    def run(self,dt):
        self._solve(dt)
        return



    def setup(self):
        from CitcomSLib import set_convection_defaults
        set_convection_defaults(self.all_variables)
	return


    def launch(self):
        from CitcomSLib import PG_timestep_init
        PG_timestep_init(self.all_variables)
        return

    #def fini(self):
	#return



    def _solve(self,dt):
        from CitcomSLib import PG_timestep_solve
        PG_timestep_solve(self.all_variables, dt)
	return



    def stable_timestep(self):
        from CitcomSLib import stable_timestep
        dt = stable_timestep(self.all_variables)
        return dt



    class Inventory(CitcomComponent.Inventory):

        import pyre.inventory as prop

        ADV = prop.bool("ADV", default=True)
        filter_temp = prop.bool("filter_temp", default=False)
        monitor_max_T = prop.bool("monitor_max_T", default=True)

        fixed_timestep = prop.float("fixed_timestep", default=0.0)
        finetunedt = prop.float("finetunedt", default=0.9)
        adv_gamma = prop.float("adv_gamma", default=0.5)
        adv_sub_iterations = prop.int("adv_sub_iterations", default=2)

        # inputdiffusivity is a deprecated cfg alias retained by setProperties.
        inputdiffusivity = prop.float("inputdiffusivity", default=-1.0)
        reference_conductivity = prop.float("reference_conductivity", default=-1.0)
        kd_mantle_thickness_km = prop.float("kd_mantle_thickness_km", default=2890.0)
        kd_transition_depth_km = prop.float("kd_transition_depth_km", default=660.0)
        kd_upper_prefactor = prop.float("kd_upper_prefactor", default=3.0)
        kd_upper_linear = prop.float("kd_upper_linear", default=15.66)
        kd_upper_quadratic = prop.float("kd_upper_quadratic", default=-16.38)
        kd_lower_prefactor = prop.float("kd_lower_prefactor", default=5.33)
        kd_lower_linear = prop.float("kd_lower_linear", default=4.98)
        kd_lower_quadratic = prop.float("kd_lower_quadratic", default=-0.81)
        kT_exponent = prop.float("kT_exponent", default=0.0)
        kC_ratio = prop.float("kC_ratio", default=0.8)




# version
__id__ = "$Id: Advection_diffusion.py 8088 2007-10-05 19:57:48Z tan2 $"

# End of file
