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


class Incompressible(CitcomComponent):


    def __init__(self, name, facility):
        CitcomComponent.__init__(self, name, facility)

        return



    def run(self):
        from CitcomSLib import general_stokes_solver
        general_stokes_solver(self.all_variables)
	return



    def setup(self):
        from CitcomSLib import set_cg_defaults, set_mg_defaults, set_mg_el_defaults
        if self.inventory.Solver == "cgrad":
            set_cg_defaults(self.all_variables)
        elif self.inventory.Solver == "multigrid":
            set_mg_defaults(self.all_variables)
        elif self.inventory.Solver == "multigrid-el":
            set_mg_el_defaults(self.all_variables)
        return



    def launch(self):
        from CitcomSLib import general_stokes_solver_setup
        general_stokes_solver_setup(self.all_variables)
        return



    #def fini(self):
	#return



    def setProperties(self, stream):
        from CitcomSLib import Incompressible_set_properties
        Incompressible_set_properties(self.all_variables, self.inventory, stream)
        return



    class Inventory(CitcomComponent.Inventory):

        import pyre.inventory as prop


        Solver = prop.str("Solver", default="cgrad",
                 validator=prop.choice(["cgrad",
                                        "multigrid",
                                        "multigrid-el"]))
        node_assemble = prop.bool("node_assemble", default=True)
        precond = prop.bool("precond", default=True)

        accuracy = prop.float("accuracy", default=1.0e-6)
        tole_compressibility = prop.float("tole_compressibility", default=1.0e-7)
        mg_cycle = prop.int("mg_cycle", default=1)
        down_heavy = prop.int("down_heavy", default=3)
        up_heavy = prop.int("up_heavy", default=3)

        vlowstep = prop.int("vlowstep", default=1000)
        vhighstep = prop.int("vhighstep", default=3)
        piterations = prop.int("piterations", default=1000)

        aug_lagr = prop.bool("aug_lagr", default=True)
        aug_number = prop.float("aug_number", default=2.0e3)

        compressible_formulation = prop.str("compressible_formulation",
                                            default="tala",
                                            validator=prop.choice(["tala",
                                                                   "ala"]))
        ala_schur_symmetry_check = prop.bool("ala_schur_symmetry_check",
                                             default=False)
        ala_schur_symmetry_tolerance = prop.float(
            "ala_schur_symmetry_tolerance", default=1.0e-3)
        ala_augmented_lagrangian_gamma = prop.float(
            "ala_augmented_lagrangian_gamma", default=0.0)
        ala_beta_element_source = prop.str(
            "ala_beta_element_source", default="supplied_average",
            validator=prop.choice(["supplied_average", "density_log_secant",
                                   "interval"]))
        ala_inner_accuracy_max = prop.float("ala_inner_accuracy_max",
                                            default=1.0e-4)
        ala_inner_accuracy_factor = prop.float("ala_inner_accuracy_factor",
                                               default=1.0e-2)
        ala_pcg_restart_interval = prop.int("ala_pcg_restart_interval",
                                            default=20)
        ala_outer_solver = prop.str("ala_outer_solver", default="pcg",
                                    validator=prop.choice(
                                        ["pcg", "fgmres",
                                         "coupled_fgmres"]))
        ala_coupled_initial_velocity_relative_tolerance = prop.float(
            "ala_coupled_initial_velocity_relative_tolerance", default=0.0)
        ala_coupled_inner_relative_tolerance = prop.float(
            "ala_coupled_inner_relative_tolerance", default=1.0e-2)
        ala_coupled_inner_max_cycles = prop.int(
            "ala_coupled_inner_max_cycles", default=200)
        ala_coupled_inner_progress_interval = prop.int(
            "ala_coupled_inner_progress_interval", default=20)
        ala_coupled_defect_corrections = prop.int(
            "ala_coupled_defect_corrections", default=0)
        ala_coupled_factor2_coarse_correction = prop.bool(
            "ala_coupled_factor2_coarse_correction", default=False)
        ala_coupled_multilevel_audit_only = prop.bool(
            "ala_coupled_multilevel_audit_only", default=False)
        ala_coupled_first_preconditioner_audit_only = prop.bool(
            "ala_coupled_first_preconditioner_audit_only", default=False)
        ala_coupled_debug_stop_iteration = prop.int(
            "ala_coupled_debug_stop_iteration", default=0)
        ala_coupled_element_vanka = prop.bool(
            "ala_coupled_element_vanka", default=False)
        ala_coupled_multilevel_vcycle = prop.bool(
            "ala_coupled_multilevel_vcycle", default=False)
        ala_coupled_multilevel_coarse_sweeps = prop.int(
            "ala_coupled_multilevel_coarse_sweeps", default=2)
        ala_coupled_multilevel_coarse_weight = prop.float(
            "ala_coupled_multilevel_coarse_weight", default=1.0)
        ala_viscosity_spectrum_diagnostics = prop.bool(
            "ala_viscosity_spectrum_diagnostics", default=False)
        ala_viscosity_spectrum_interval = prop.int(
            "ala_viscosity_spectrum_interval", default=1)
        ala_coupled_shallow_vanka_layers = prop.int(
            "ala_coupled_shallow_vanka_layers", default=0)
        ala_coupled_shallow_vanka_core_layers = prop.int(
            "ala_coupled_shallow_vanka_core_layers", default=-1)
        ala_coupled_shallow_vanka_band_sweeps = prop.int(
            "ala_coupled_shallow_vanka_band_sweeps", default=0)
        ala_coupled_shallow_vanka_sweeps = prop.int(
            "ala_coupled_shallow_vanka_sweeps", default=0)
        ala_coupled_shallow_vanka_warm_sweeps = prop.int(
            "ala_coupled_shallow_vanka_warm_sweeps", default=-1)
        ala_unaugmented_momentum_tolerance = prop.float(
            "ala_unaugmented_momentum_tolerance", default=0.0)
        ala_feasibility_audit = prop.bool("ala_feasibility_audit",
                                          default=False)
        ala_feasibility_window = prop.int("ala_feasibility_window",
                                           default=20)
        ala_feasibility_min_reduction = prop.float(
            "ala_feasibility_min_reduction", default=0.02)
        ala_hybrid_convergence = prop.bool("ala_hybrid_convergence",
                                           default=False)
        ala_div_v_tolerance = prop.float("ala_div_v_tolerance",
                                         default=1.0e-7)
        ala_update_tolerance = prop.float("ala_update_tolerance",
                                          default=1.0e-3)
        ala_consecutive_steps = prop.int("ala_consecutive_steps", default=3)
        ala_depth_diagnostics = prop.bool("ala_depth_diagnostics",
                                           default=False)
        ala_depth_diagnostic_interval = prop.int(
            "ala_depth_diagnostic_interval", default=5)
        ala_depth_diagnostic_bins = prop.int("ala_depth_diagnostic_bins",
                                              default=8)
        ala_beta_causal_diagnostics = prop.bool(
            "ala_beta_causal_diagnostics", default=False)
        ala_coarse_residual_diagnostics = prop.bool(
            "ala_coarse_residual_diagnostics", default=False)
        ala_coarse_residual_interval = prop.int(
            "ala_coarse_residual_interval", default=5)
        ala_coarse_residual_levels = prop.int(
            "ala_coarse_residual_levels", default=3)
        ala_two_level_preconditioner = prop.bool(
            "ala_two_level_preconditioner", default=False)
        ala_two_level_offset = prop.int("ala_two_level_offset", default=2)
        ala_two_level_coarse_iterations = prop.int(
            "ala_two_level_coarse_iterations", default=12)
        ala_two_level_coarse_damping = prop.float(
            "ala_two_level_coarse_damping", default=0.03)
        ala_two_level_coarse_solver = prop.str(
            "ala_two_level_coarse_solver", default="chebyshev",
            validator=prop.choice(
                ["jacobi", "chebyshev", "scaled_diagonal"]))
        ala_two_level_coarse_eigenvalue_min = prop.float(
            "ala_two_level_coarse_eigenvalue_min", default=0.01)
        ala_two_level_coarse_eigenvalue_max = prop.float(
            "ala_two_level_coarse_eigenvalue_max", default=27.0)
        ala_two_level_coarse_weight = prop.float(
            "ala_two_level_coarse_weight", default=0.05)
        ala_two_level_velocity_solver = prop.str(
            "ala_two_level_velocity_solver", default="chebyshev",
            validator=prop.choice(["diagonal", "chebyshev"]))
        ala_two_level_velocity_iterations = prop.int(
            "ala_two_level_velocity_iterations", default=16)
        ala_two_level_velocity_eigenvalue_min = prop.float(
            "ala_two_level_velocity_eigenvalue_min", default=0.01)
        ala_two_level_velocity_eigenvalue_max = prop.float(
            "ala_two_level_velocity_eigenvalue_max", default=4.0)
        ala_pressure_multigrid = prop.bool(
            "ala_pressure_multigrid", default=False)
        ala_pressure_multigrid_galerkin = prop.bool(
            "ala_pressure_multigrid_galerkin", default=False)
        ala_pressure_bpi_weight = prop.float(
            "ala_pressure_bpi_weight", default=1.0)
        ala_pressure_factor2_coarse_action_scale = prop.float(
            "ala_pressure_factor2_coarse_action_scale", default=1.0)
        ala_pressure_multigrid_min_level = prop.int(
            "ala_pressure_multigrid_min_level", default=0)
        ala_pressure_multigrid_pre_smooth = prop.int(
            "ala_pressure_multigrid_pre_smooth", default=2)
        ala_pressure_multigrid_post_smooth = prop.int(
            "ala_pressure_multigrid_post_smooth", default=2)
        ala_pressure_multigrid_coarse_iterations = prop.int(
            "ala_pressure_multigrid_coarse_iterations", default=20)
        ala_pressure_multigrid_damping = prop.float(
            "ala_pressure_multigrid_damping", default=0.04)
        ala_pressure_multigrid_weight = prop.float(
            "ala_pressure_multigrid_weight", default=0.5)
        ala_global_coarse_preconditioner = prop.bool(
            "ala_global_coarse_preconditioner", default=False)
        ala_global_coarse_weight = prop.float(
            "ala_global_coarse_weight", default=1.0)
        ala_global_coarse_regularization = prop.float(
            "ala_global_coarse_regularization", default=1.0e-10)
        ala_shallow_patch_preconditioner = prop.bool(
            "ala_shallow_patch_preconditioner", default=False)
        ala_shallow_patch_depth_km = prop.float(
            "ala_shallow_patch_depth_km", default=410.0)
        ala_shallow_patch_weight = prop.float(
            "ala_shallow_patch_weight", default=0.25)
        ala_shallow_patch_regularization = prop.float(
            "ala_shallow_patch_regularization", default=1.0e-3)
        ala_shallow_patch_horizontal_elements = prop.int(
            "ala_shallow_patch_horizontal_elements", default=4)
        ala_shallow_patch_horizontal_stride = prop.int(
            "ala_shallow_patch_horizontal_stride", default=2)
        ala_shallow_patch_mpi_overlap = prop.int(
            "ala_shallow_patch_mpi_overlap", default=2)
        ala_shallow_patch_velocity_solver = prop.str(
            "ala_shallow_patch_velocity_solver", default="diagonal",
            validator=prop.choice(
                ["diagonal", "node_block", "element_vanka"]))
        ala_geneo_preconditioner = prop.bool(
            "ala_geneo_preconditioner", default=False)
        ala_geneo_basis_type = prop.str(
            "ala_geneo_basis_type", default="spectral",
            validator=prop.choice(["spectral", "radial_partition"]))
        ala_geneo_eigenvalue_threshold = prop.float(
            "ala_geneo_eigenvalue_threshold", default=0.20)
        ala_geneo_min_modes_per_rank = prop.int(
            "ala_geneo_min_modes_per_rank", default=1)
        ala_geneo_max_modes_per_rank = prop.int(
            "ala_geneo_max_modes_per_rank", default=2)
        ala_geneo_horizontal_bins = prop.int(
            "ala_geneo_horizontal_bins", default=4)
        ala_geneo_radial_bins = prop.int(
            "ala_geneo_radial_bins", default=2)
        ala_geneo_rank_group_x = prop.int(
            "ala_geneo_rank_group_x", default=1)
        ala_geneo_rank_group_y = prop.int(
            "ala_geneo_rank_group_y", default=1)
        ala_geneo_max_global_modes = prop.int(
            "ala_geneo_max_global_modes", default=400)
        ala_geneo_weight = prop.float("ala_geneo_weight", default=1.0)
        ala_geneo_regularization = prop.float(
            "ala_geneo_regularization", default=1.0e-8)
        ala_radial_line_preconditioner = prop.bool(
            "ala_radial_line_preconditioner", default=False)
        ala_element_vanka_smoother = prop.bool(
            "ala_element_vanka_smoother", default=False)
        ala_element_vanka_damping = prop.float(
            "ala_element_vanka_damping", default=0.8)
        ala_element_vanka_pressure_damping = prop.float(
            "ala_element_vanka_pressure_damping", default=1.0)
        ala_element_vanka_regularization = prop.float(
            "ala_element_vanka_regularization", default=1.0e-8)
        ala_element_vanka_external_diagonal_weight = prop.float(
            "ala_element_vanka_external_diagonal_weight", default=1.0)
        ala_element_vanka_galerkin_schur = prop.bool(
            "ala_element_vanka_galerkin_schur", default=False)
        ala_element_vanka_rebuild_interval = prop.int(
            "ala_element_vanka_rebuild_interval", default=1)
        uzawa = prop.str("uzawa", default="cg",
                         validator=prop.choice(["cg", "bicg", "ala_cg"]))
        compress_iter_maxstep = prop.int("compress_iter_maxstep", default=100)
        relative_err_accuracy = prop.float("relative_err_accuracy", default=0.001)
        remove_rigid_rotation = prop.bool("remove_rigid_rotation", default=True)

# version
__id__ = "$Id: Incompressible.py 8095 2007-10-10 20:11:00Z tan2 $"

# End of file
