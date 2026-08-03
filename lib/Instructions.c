/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/* Set up the finite element problem to suit: returns with all memory */
/* allocated, temperature, viscosity, node locations and how to use */
/* them all established. 8.29.92 or 29.8.92 depending on your nationality*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <sys/stat.h>
#include <sys/errno.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include "element_definitions.h"
#include "global_defs.h"

#include "citcom_init.h"
#include "initial_temperature.h"
#include "lith_age.h"
#include "material_properties.h"
#include "output.h"
#include "output_h5.h"
#include "parallel_related.h"
#include "parsing.h"
#include "phase_change.h"
#include "interuption.h"

void parallel_process_termination();
void allocate_common_vars(struct All_variables*);
void allocate_velocity_vars(struct All_variables*);
void check_bc_consistency(struct All_variables*);
void construct_elt_gs(struct All_variables*);
void construct_elt_cs(struct All_variables*);
void construct_id(struct All_variables*);
void construct_ien(struct All_variables*);
void construct_lm(struct All_variables*);
void construct_masks(struct All_variables*);
void construct_shape_functions(struct All_variables*);
void construct_sub_element(struct All_variables*);
void construct_surf_det (struct All_variables*);
void construct_bdry_det (struct All_variables*);
void construct_surface (struct All_variables*);
void get_initial_elapsed_time(struct All_variables*);
void lith_age_init(struct All_variables *E);
void mass_matrix(struct All_variables*);
void output_init(struct All_variables*);
void set_elapsed_time(struct All_variables*);
void set_sphere_harmonics (struct All_variables*);
void set_starting_age(struct All_variables*);
void tracer_initial_settings(struct All_variables*);
void tracer_input(struct All_variables*);
void viscosity_input(struct All_variables*);
void get_vtk_filename(char *,int,struct All_variables *,int);
void myerror(struct All_variables *,char *);
static void expand_str(char *, size_t, const char *, const char *);
void open_qfiles(struct All_variables *) ;


void initial_mesh_solver_setup(struct All_variables *E)
{

    E->monitor.cpu_time_at_last_cycle =
        E->monitor.cpu_time_at_start = CPU_time0();

    output_init(E);
    (E->problem_derived_values)(E);   /* call this before global_derived_  */
    (E->solver.global_derived_values)(E);

    (E->solver.parallel_processor_setup)(E);   /* get # of proc in x,y,z */
    (E->solver.parallel_domain_decomp0)(E);  /* get local nel, nno, elx, nox et al */
    allocate_common_vars(E);
    (E->problem_allocate_vars)(E);
    (E->solver_allocate_vars)(E);

           /* logical domain */
    construct_ien(E);
    construct_surface(E);
    (E->solver.construct_boundary)(E);
    (E->solver.parallel_domain_boundary_nodes)(E);

           /* physical domain */
    (E->solver.node_locations)(E);
    allocate_velocity_vars(E);
   // fprintf(stderr,"After allocate_velo_vars\n");

    get_initial_elapsed_time(E);  /* Set elapsed time */
    set_starting_age(E);  /* set the starting age to elapsed time, if desired */
    set_elapsed_time(E);         /* reset to elapsed time to zero, if desired */


    /* open the heatflow files here because we need to know about loc_me */
    if(E->output.write_q_files)
      open_qfiles(E);
    else{
      E->output.fpqt = E->output.fpqb = NULL;
    }


   // fprintf(stderr,"Before lith_age_init\n");
    if(E->control.lith_age)
        lith_age_init(E);
   // fprintf(stderr,"After lith_age_init\n");

    (E->problem_boundary_conds)(E);

    check_bc_consistency(E);

    construct_masks(E);		/* order is important here */
    construct_id(E);
    construct_lm(E);

    (E->solver.parallel_communication_routs_v)(E);
    (E->solver.parallel_communication_routs_s)(E);

    reference_state(E);

    construct_sub_element(E);
    construct_shape_functions(E);
    construct_elt_gs(E);
    if(E->control.inv_gruneisen != 0)
        construct_elt_cs(E);

    /* this matrix results from spherical geometry */
    /* construct_c3x3matrix(E); */

    mass_matrix(E);
   // fprintf(stderr,"After mass_matrix\n");
    construct_surf_det (E);
    construct_bdry_det (E);

    set_sphere_harmonics (E);

    if(E->control.tracer) {
	tracer_initial_settings(E);
	//fprintf(stderr,"After tracer_initial_settings\n");
	(E->problem_tracer_setup)(E);
	//fprintf(stderr,"After tracer_setup\n");
    }

}


void read_instructions(struct All_variables *E, char *filename)
{
    void read_initial_settings();
    void global_default_values();

    void setup_parser();
    void shutdown_parser();

    /* =====================================================
       Global interuption handling routine defined once here
       =====================================================  */

    set_signal();

    /* ==================================================
       Initialize from the command line
       from startup files. (See Parsing.c).
       ==================================================  */

    setup_parser(E,filename);

    global_default_values(E);
    read_initial_settings(E);

    shutdown_parser(E);

    return;
}


/* This function is replaced by CitcomS.Solver._setup() */
void initial_setup(struct All_variables *E)
{
    void general_stokes_solver_setup();
    void initial_mesh_solver_setup();

    initial_mesh_solver_setup(E);

    general_stokes_solver_setup(E);

    (E->next_buoyancy_field_init)(E);


    if (E->parallel.me==0) fprintf(stderr,"time=%f\n",
                                   CPU_time0()-E->monitor.cpu_time_at_start);

    return;
}


void initialize_material(struct All_variables *E)
{
    void construct_mat_group();
    void read_mat_from_file();

    if(E->control.mat_control)
        read_mat_from_file(E);
    else
        construct_mat_group(E);
}


/* This function is replaced by CitcomS.Components.IC.launch()*/
void initial_conditions(struct All_variables *E)
{
    void initialize_tracers();
    void init_composition();
    void common_initial_fields();

    initialize_material(E);

    if (E->control.tracer==1) {
        initialize_tracers(E);

        if (E->composition.on)
            init_composition(E);
    }

    (E->problem_initial_fields)(E);   /* temperature/chemistry/melting etc */
    common_initial_fields(E);  /* velocity/pressure/viscosity (viscosity must be done LAST) */

    return;
}


void read_initial_settings(struct All_variables *E)
{
  void set_convection_defaults();
  void set_cg_defaults();
  void set_mg_defaults();
  float tmp;
  int m=E->parallel.me;

  /* first the problem type (defines subsequent behaviour) */

  input_string("Problem",E->control.PROBLEM_TYPE,"convection",m);
  if ( strcmp(E->control.PROBLEM_TYPE,"convection") == 0)  {
    E->control.CONVECTION = 1;
    set_convection_defaults(E);
  }

  else if ( strcmp(E->control.PROBLEM_TYPE,"convection-chemical") == 0) {
    E->control.CONVECTION = 1;
    set_convection_defaults(E);
  }

  else {
    fprintf(E->fp,"Unable to determine problem type, assuming convection ... \n");
    E->control.CONVECTION = 1;
    set_convection_defaults(E);
  }

  input_string("Geometry",E->control.GEOMETRY,"sphere",m);
  if ( strcmp(E->control.GEOMETRY,"cart2d") == 0)
    { E->control.CART2D = 1;
    (E->solver.set_2dc_defaults)(E);}
  else if ( strcmp(E->control.GEOMETRY,"axi") == 0)
    { E->control.AXI = 1;
    }
  else if ( strcmp(E->control.GEOMETRY,"cart2pt5d") == 0)
    { E->control.CART2pt5D = 1;
    (E->solver.set_2pt5dc_defaults)(E);}
  else if ( strcmp(E->control.GEOMETRY,"cart3d") == 0)
    { E->control.CART3D = 1;
    (E->solver.set_3dc_defaults)(E);}
  else if ( strcmp(E->control.GEOMETRY,"sphere") == 0)
    {
      (E->solver.set_3dsphere_defaults)(E);}
  else
    { fprintf(E->fp,"Unable to determine geometry, assuming cartesian 2d ... \n");
    E->control.CART2D = 1;
    (E->solver.set_2dc_defaults)(E); }

  input_string("Solver",E->control.SOLVER_TYPE,"cgrad",m);
  if ( strcmp(E->control.SOLVER_TYPE,"cgrad") == 0)
    { E->control.CONJ_GRAD = 1;
    set_cg_defaults(E);}
  else if ( strcmp(E->control.SOLVER_TYPE,"multigrid") == 0)
    { E->control.NMULTIGRID = 1;
    set_mg_defaults(E);}
  else if ( strcmp(E->control.SOLVER_TYPE,"multigrid-el") == 0)
    { E->control.EMULTIGRID = 1;
    set_mg_defaults(E);}
  else
    { if (E->parallel.me==0) fprintf(stderr,"Unable to determine how to solve, specify Solver=VALID_OPTION \n");
    parallel_process_termination();
    }


  /* admin */

  input_string("Spacing",E->control.NODE_SPACING,"regular",m);
  if ( strcmp(E->control.NODE_SPACING,"regular") == 0)
    E->control.GRID_TYPE = 1;
  else if ( strcmp(E->control.NODE_SPACING,"bound_lyr") == 0)
    E->control.GRID_TYPE = 2;
  else if ( strcmp(E->control.NODE_SPACING,"region") == 0)
    E->control.GRID_TYPE = 3;
  else if ( strcmp(E->control.NODE_SPACING,"ortho_files") == 0)
    E->control.GRID_TYPE = 4;
  else
    {  E->control.GRID_TYPE = 1; }

  /* Information on which files to print, which variables of the flow to calculate and print.
     Default is no information recorded (apart from special things for given applications.
  */

  input_string("datadir",E->control.data_dir,".",m);
  input_string("datafile",E->control.data_prefix,"initialize",m);
  input_string("datadir_old",E->control.data_dir_old,".",m);
  input_string("datafile_old",E->control.data_prefix_old,"initialize",m);
  input_string("logfile",E->control.log_template,"datafile",m);

  input_int("mgunitx",&(E->mesh.mgunitx),"1",m);
  input_int("mgunitz",&(E->mesh.mgunitz),"1",m);
  input_int("mgunity",&(E->mesh.mgunity),"1",m);
  input_int("levels",&(E->mesh.levels),"0",m);

  input_int("coor",&(E->control.coor),"0",m);
  if(E->control.coor == 2){
    /*
       refinement in two layers
    */
    /* number of refinement layers */
    E->control.coor_refine[0] = 0.10; /* bottom 10% */
    E->control.coor_refine[1] = 0.15; /* get 15% of the nodes */
    E->control.coor_refine[2] = 0.10; /* top 10% */
    E->control.coor_refine[3] = 0.20; /* get 20% of the nodes */
    input_float_vector("coor_refine",4,E->control.coor_refine,m);
  }else if(E->control.coor == 3){
    /*

    refinement CitcomCU style, by reading in layers, e.g.

	r_grid_layers=3		# minus 1 is number of layers with uniform grid in r
	rr=0.5,0.75,1.0 	#    starting and ending r coodinates
	nr=1,37,97		#    starting and ending node in r direction

    */
    input_int("r_grid_layers", &(E->control.rlayers), "1",m);
    if(E->control.rlayers > 20)
      myerror(E,"number of rlayers out of bounds (20) for coor = 3");
    /* layers radii */
    input_float_vector("rr", E->control.rlayers, (E->control.rrlayer),m);
    /* associated node numbers */
    input_int_vector("nr", E->control.rlayers, (E->control.nrlayer),m);
   }

  input_string("coor_file",E->control.coor_file,"",m);

  input_int("nprocx",&(E->parallel.nprocx),"1",m);
  input_int("nprocy",&(E->parallel.nprocy),"1",m);
  input_int("nprocz",&(E->parallel.nprocz),"1",m);
  input_int("nproc_surf",&(E->parallel.nprocxy),"1",m);


  input_boolean("node_assemble",&(E->control.NASSEMBLE),"off",m);
  /* general mesh structure */

  input_boolean("verbose",&(E->control.verbose),"off",m);
  input_boolean("see_convergence",&(E->control.print_convergence),"off",m);

  input_int("stokes_flow_only",&(E->control.stokes),"0",m);

  /* restart from checkpoint file */
  input_boolean("restart",&(E->control.restart),"off",m);
  input_int("post_p",&(E->control.post_p),"0",m);
  input_int("solution_cycles_init",&(E->monitor.solution_cycles_init),"0",m);

  /* for layers    */
  /*

  these boundaries are a little wacko


  */
  input_float("z_cmb",&(E->viscosity.zcmb),"0.45",m); /* does this ever get used? */
  input_float("z_lmantle",&(E->viscosity.zlm),"0.45",m);
  input_float("z_410",&(E->viscosity.z410),"0.225",m); /* 0.06434, more like it */
  input_float("z_lith",&(E->viscosity.zlith),"0.225",m); /* 0.0157, more like it */

  /*  the start age and initial subduction history   */
  input_float("start_age",&(E->control.start_age),"0.0",m);
  input_int("reset_startage",&(E->control.reset_startage),"0",m);
  input_int("zero_elapsed_time",&(E->control.zero_elapsed_time),"0",m);

  input_int("output_ll_max",&(E->output.llmax),"1",m);

  input_int("topvbc",&(E->mesh.topvbc),"0",m);
  input_int("botvbc",&(E->mesh.botvbc),"0",m);

  input_float("topvbxval",&(E->control.VBXtopval),"0.0",m);
  input_float("botvbxval",&(E->control.VBXbotval),"0.0",m);
  input_float("topvbyval",&(E->control.VBYtopval),"0.0",m);
  input_float("botvbyval",&(E->control.VBYbotval),"0.0",m);


  input_float("T_interior_max_for_exit",&(E->monitor.T_interior_max_for_exit),"1.5",m);

  input_int("pseudo_free_surf",&(E->control.pseudo_free_surf),"0",m);

  input_int("toptbc",&(E->mesh.toptbc),"1",m);
  input_int("bottbc",&(E->mesh.bottbc),"1",m);
  input_float("toptbcval",&(E->control.TBCtopval),"0.0",m);
  input_float("bottbcval",&(E->control.TBCbotval),"1.0",m);

  input_boolean("side_sbcs",&(E->control.side_sbcs),"off",m);

  input_int("file_vbcs",&(E->control.vbcs_file),"0",m);
  input_string("vel_bound_file",E->control.velocity_boundary_file,"",m);

  input_int("reference_state",&(E->refstate.choice),"1",m);
  if(E->refstate.choice == 0) {
      input_string("refstate_file",E->refstate.filename,"refstate.dat",m);
      input_string("ala_beta_interval_file",
                   E->refstate.beta_interval_filename,
                   "interval_ALA_strict.txt",m);
  }

  input_int("mat_control",&(E->control.mat_control),"0",m);
  input_string("mat_file",E->control.mat_file,"",m);

  input_int("nodex",&(E->mesh.nox),"essential",m);
  input_int("nodez",&(E->mesh.noz),"essential",m);
  input_int("nodey",&(E->mesh.noy),"essential",m);

  input_boolean("aug_lagr",&(E->control.augmented_Lagr),"off",m);
  input_double("aug_number",&(E->control.augmented),"0.0",m);

  input_boolean("remove_rigid_rotation",&(E->control.remove_rigid_rotation),"on",m);

  input_float("tole_compressibility",&(E->control.tole_comp),"0.0",m);

  input_int("storage_spacing",&(E->control.record_every),"10",m);
  input_int("checkpointFrequency",&(E->control.checkpoint_frequency),"100",m);
  input_int("cpu_limits_in_seconds",&(E->control.record_all_until),"5",m);
  input_int("write_q_files",&(E->output.write_q_files),"0",m);/* write additional
								 heat flux files? */
  if(E->output.write_q_files){	/* make sure those get written at
				   least as often as velocities */
    E->output.write_q_files = min(E->output.write_q_files,E->control.record_every);
  }

  input_int("cmbhf_CBF_freq",&(E->output.cmbhf_CBF_freq),"0",m); /* legacy C-side CBF heat flux
                                                                      frequency (0 = disabled) */
  input_boolean("cbf_output_shflux",&(E->output.cbf_output_shflux),"on",m); /* write surface/top
                                                                                shflux_CBF files */
  input_boolean("cbf_output_bhflux",&(E->output.cbf_output_bhflux),"on",m); /* write bottom/CMB
                                                                                bhflux_CBF files */
  input_boolean("cbf_use_advection",&(E->output.cbf_use_advection),"on",m); /* include u.gradT
                                                                                in CBF. off is an
                                                                                approximate diagnostic,
                                                                                not equation-consistent */


  input_boolean("precond",&(E->control.precondition),"off",m);
  input_int("mg_cycle",&(E->control.mg_cycle),"2,0,nomax",m);
  input_int("down_heavy",&(E->control.down_heavy),"1,0,nomax",m);
  input_int("up_heavy",&(E->control.up_heavy),"1,0,nomax",m);
  input_double("accuracy",&(E->control.accuracy),"1.0e-4,0.0,1.0",m);

  input_int("vhighstep",&(E->control.v_steps_high),"1,0,nomax",m);
  input_int("vlowstep",&(E->control.v_steps_low),"250,0,nomax",m);
  input_int("piterations",&(E->control.p_iterations),"100,0,nomax",m);

  input_float("rayleigh",&(E->control.Atemp),"essential",m);
  input_float("dissipation_number",&(E->control.disptn_number),"0.0",m);
  input_float("gruneisen",&(tmp),"0.0",m);
  /* special case: if tmp==0, set gruneisen as inf */
  if(tmp != 0)
      E->control.inv_gruneisen = 1/tmp;
  else
      E->control.inv_gruneisen = 0;

  input_string("compressible_formulation",
               E->control.compressible_formulation,"tala",m);
  if(strcmp(E->control.compressible_formulation, "tala") == 0) {
      E->control.ala_pressure_buoyancy = 0;
  }
  else if(strcmp(E->control.compressible_formulation, "ala") == 0) {
      E->control.ala_pressure_buoyancy = 1;
      if(E->control.inv_gruneisen == 0)
          myerror(E, "compressible_formulation=ala requires gruneisen != 0");
  }
  else {
      myerror(E, "compressible_formulation must be tala or ala");
  }

  input_boolean("ala_schur_symmetry_check",
                &(E->control.ala_schur_symmetry_check),"off",m);
  input_double("ala_schur_symmetry_tolerance",
               &(E->control.ala_schur_symmetry_tolerance),"1.0e-3",m);
  input_double("ala_augmented_lagrangian_gamma",
               &(E->control.ala_augmented_lagrangian_gamma),"0.0",m);
  input_string("ala_beta_element_source",
               E->control.ala_beta_element_source,"supplied_average",m);
  input_double("ala_inner_accuracy_max",
               &(E->control.ala_inner_accuracy_max),"1.0e-4",m);
  input_double("ala_inner_accuracy_factor",
               &(E->control.ala_inner_accuracy_factor),"1.0e-2",m);
  input_int("ala_pcg_restart_interval",
            &(E->control.ala_pcg_restart_interval),"20",m);
  input_string("ala_outer_solver",E->control.ala_outer_solver,"pcg",m);
  input_double("ala_unaugmented_momentum_tolerance",
               &(E->control.ala_unaugmented_momentum_tolerance),"0.0",m);
  input_boolean("ala_feasibility_audit",
                &(E->control.ala_feasibility_audit),"off",m);
  input_int("ala_feasibility_window",
            &(E->control.ala_feasibility_window),"20",m);
  input_double("ala_feasibility_min_reduction",
               &(E->control.ala_feasibility_min_reduction),"0.02",m);
  input_boolean("ala_hybrid_convergence",
                &(E->control.ala_hybrid_convergence),"off",m);
  input_double("ala_div_v_tolerance",
               &(E->control.ala_div_v_tolerance),"1.0e-7",m);
  input_double("ala_update_tolerance",
               &(E->control.ala_update_tolerance),"1.0e-3",m);
  input_int("ala_consecutive_steps",
            &(E->control.ala_consecutive_steps),"3",m);
  input_boolean("ala_depth_diagnostics",
                &(E->control.ala_depth_diagnostics),"off",m);
  input_int("ala_depth_diagnostic_interval",
            &(E->control.ala_depth_diagnostic_interval),"5",m);
  input_int("ala_depth_diagnostic_bins",
            &(E->control.ala_depth_diagnostic_bins),"8",m);
  input_boolean("ala_beta_causal_diagnostics",
                &(E->control.ala_beta_causal_diagnostics),"off",m);
  input_boolean("ala_coarse_residual_diagnostics",
                &(E->control.ala_coarse_residual_diagnostics),"off",m);
  input_int("ala_coarse_residual_interval",
            &(E->control.ala_coarse_residual_interval),"5",m);
  input_int("ala_coarse_residual_levels",
            &(E->control.ala_coarse_residual_levels),"3",m);
  input_boolean("ala_two_level_preconditioner",
                &(E->control.ala_two_level_preconditioner),"off",m);
  input_int("ala_two_level_offset",
            &(E->control.ala_two_level_offset),"2",m);
  input_int("ala_two_level_coarse_iterations",
            &(E->control.ala_two_level_coarse_iterations),"12",m);
  input_double("ala_two_level_coarse_damping",
               &(E->control.ala_two_level_coarse_damping),"0.03",m);
  input_string("ala_two_level_coarse_solver",
               E->control.ala_two_level_coarse_solver,"chebyshev",m);
  input_double("ala_two_level_coarse_eigenvalue_min",
               &(E->control.ala_two_level_coarse_eigenvalue_min),"0.01",m);
  input_double("ala_two_level_coarse_eigenvalue_max",
               &(E->control.ala_two_level_coarse_eigenvalue_max),"27.0",m);
  input_double("ala_two_level_coarse_weight",
               &(E->control.ala_two_level_coarse_weight),"0.05",m);
  input_string("ala_two_level_velocity_solver",
               E->control.ala_two_level_velocity_solver,"chebyshev",m);
  input_int("ala_two_level_velocity_iterations",
            &(E->control.ala_two_level_velocity_iterations),"16",m);
  input_double("ala_two_level_velocity_eigenvalue_min",
               &(E->control.ala_two_level_velocity_eigenvalue_min),"0.01",m);
  input_double("ala_two_level_velocity_eigenvalue_max",
               &(E->control.ala_two_level_velocity_eigenvalue_max),"4.0",m);
  input_boolean("ala_pressure_multigrid",
                &(E->control.ala_pressure_multigrid),"off",m);
  input_int("ala_pressure_multigrid_min_level",
            &(E->control.ala_pressure_multigrid_min_level),"0",m);
  input_int("ala_pressure_multigrid_pre_smooth",
            &(E->control.ala_pressure_multigrid_pre_smooth),"2",m);
  input_int("ala_pressure_multigrid_post_smooth",
            &(E->control.ala_pressure_multigrid_post_smooth),"2",m);
  input_int("ala_pressure_multigrid_coarse_iterations",
            &(E->control.ala_pressure_multigrid_coarse_iterations),"20",m);
  input_double("ala_pressure_multigrid_damping",
               &(E->control.ala_pressure_multigrid_damping),"0.04",m);
  input_double("ala_pressure_multigrid_weight",
               &(E->control.ala_pressure_multigrid_weight),"0.5",m);
  input_boolean("ala_global_coarse_preconditioner",
                &(E->control.ala_global_coarse_preconditioner),"off",m);
  input_double("ala_global_coarse_weight",
               &(E->control.ala_global_coarse_weight),"1.0",m);
  input_double("ala_global_coarse_regularization",
               &(E->control.ala_global_coarse_regularization),"1.0e-10",m);
  input_boolean("ala_shallow_patch_preconditioner",
                &(E->control.ala_shallow_patch_preconditioner),"off",m);
  input_double("ala_shallow_patch_depth_km",
               &(E->control.ala_shallow_patch_depth_km),"410.0",m);
  input_double("ala_shallow_patch_weight",
               &(E->control.ala_shallow_patch_weight),"0.25",m);
  input_double("ala_shallow_patch_regularization",
               &(E->control.ala_shallow_patch_regularization),"1.0e-3",m);
  input_int("ala_shallow_patch_horizontal_elements",
            &(E->control.ala_shallow_patch_horizontal_elements),"4",m);
  input_int("ala_shallow_patch_horizontal_stride",
            &(E->control.ala_shallow_patch_horizontal_stride),"2",m);
  input_int("ala_shallow_patch_mpi_overlap",
            &(E->control.ala_shallow_patch_mpi_overlap),"2",m);
  input_string("ala_shallow_patch_velocity_solver",
               E->control.ala_shallow_patch_velocity_solver,"diagonal",m);
  input_boolean("ala_geneo_preconditioner",
                &(E->control.ala_geneo_preconditioner),"off",m);
  input_string("ala_geneo_basis_type",
               E->control.ala_geneo_basis_type,"spectral",m);
  input_double("ala_geneo_eigenvalue_threshold",
               &(E->control.ala_geneo_eigenvalue_threshold),"0.20",m);
  input_int("ala_geneo_min_modes_per_rank",
            &(E->control.ala_geneo_min_modes_per_rank),"1",m);
  input_int("ala_geneo_max_modes_per_rank",
            &(E->control.ala_geneo_max_modes_per_rank),"2",m);
  input_int("ala_geneo_horizontal_bins",
            &(E->control.ala_geneo_horizontal_bins),"4",m);
  input_int("ala_geneo_radial_bins",
            &(E->control.ala_geneo_radial_bins),"2",m);
  input_int("ala_geneo_rank_group_x",
            &(E->control.ala_geneo_rank_group_x),"1",m);
  input_int("ala_geneo_rank_group_y",
            &(E->control.ala_geneo_rank_group_y),"1",m);
  input_int("ala_geneo_max_global_modes",
            &(E->control.ala_geneo_max_global_modes),"400",m);
  input_double("ala_geneo_weight",
               &(E->control.ala_geneo_weight),"1.0",m);
  input_double("ala_geneo_regularization",
               &(E->control.ala_geneo_regularization),"1.0e-8",m);
  input_boolean("ala_radial_line_preconditioner",
                &(E->control.ala_radial_line_preconditioner),"off",m);
  input_boolean("ala_element_vanka_smoother",
                &(E->control.ala_element_vanka_smoother),"off",m);
  input_double("ala_element_vanka_damping",
               &(E->control.ala_element_vanka_damping),"0.8",m);
  input_double("ala_element_vanka_regularization",
               &(E->control.ala_element_vanka_regularization),"1.0e-8",m);

  if(E->control.ala_schur_symmetry_tolerance <= 0.0)
      myerror(E, "ala_schur_symmetry_tolerance must be positive");
  if(E->control.ala_augmented_lagrangian_gamma < 0.0)
      myerror(E, "ala_augmented_lagrangian_gamma must be nonnegative");
  if(E->control.ala_augmented_lagrangian_gamma > 0.0 &&
     !E->control.ala_pressure_buoyancy)
      myerror(E, "ala_augmented_lagrangian_gamma requires "
              "compressible_formulation=ala");
  if(E->control.ala_augmented_lagrangian_gamma > 0.0 &&
     E->control.augmented_Lagr)
      myerror(E, "ala_augmented_lagrangian_gamma and aug_lagr are "
              "mutually exclusive");
  if(strcmp(E->control.ala_beta_element_source,"supplied_average") != 0 &&
     strcmp(E->control.ala_beta_element_source,"density_log_secant") != 0 &&
     strcmp(E->control.ala_beta_element_source,"interval") != 0)
      myerror(E, "ala_beta_element_source must be supplied_average, "
              "density_log_secant, or interval");
  if(strcmp(E->control.ala_beta_element_source,"interval") == 0 &&
     (!E->control.ala_pressure_buoyancy || E->refstate.choice != 0))
      myerror(E, "ala_beta_element_source=interval requires "
              "compressible_formulation=ala and reference_state=0");
  if(E->control.ala_inner_accuracy_max <= 0.0)
      myerror(E, "ala_inner_accuracy_max must be positive");
  if(E->control.ala_inner_accuracy_factor <= 0.0)
      myerror(E, "ala_inner_accuracy_factor must be positive");
  if(E->control.ala_pcg_restart_interval < 1)
      myerror(E, "ala_pcg_restart_interval must be at least one");
  if(strcmp(E->control.ala_outer_solver,"pcg") != 0 &&
     strcmp(E->control.ala_outer_solver,"fgmres") != 0 &&
     strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)
      myerror(E, "ala_outer_solver must be pcg, fgmres, or coupled_fgmres");
  if(E->control.ala_unaugmented_momentum_tolerance < 0.0)
      myerror(E, "ala_unaugmented_momentum_tolerance must be nonnegative");
  if(E->control.ala_unaugmented_momentum_tolerance > 0.0 &&
     strcmp(E->control.ala_outer_solver,"fgmres") != 0 &&
     strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)
      myerror(E, "ala_unaugmented_momentum_tolerance requires "
              "ala_outer_solver=fgmres or coupled_fgmres");
  if(E->control.ala_feasibility_window < 1)
      myerror(E, "ala_feasibility_window must be at least one");
  if(E->control.ala_feasibility_min_reduction < 0.0 ||
     E->control.ala_feasibility_min_reduction >= 1.0)
      myerror(E, "ala_feasibility_min_reduction must be in [0,1)");
  if(E->control.ala_div_v_tolerance <= 0.0)
      myerror(E, "ala_div_v_tolerance must be positive");
  if(E->control.ala_update_tolerance <= 0.0)
      myerror(E, "ala_update_tolerance must be positive");
  if(E->control.ala_consecutive_steps < 1)
      myerror(E, "ala_consecutive_steps must be at least one");
  if(E->control.ala_depth_diagnostic_interval < 1)
      myerror(E, "ala_depth_diagnostic_interval must be at least one");
  if(E->control.ala_depth_diagnostic_bins < 1 ||
     E->control.ala_depth_diagnostic_bins > 128)
      myerror(E, "ala_depth_diagnostic_bins must be between 1 and 128");
  if(E->control.ala_coarse_residual_interval < 1)
      myerror(E, "ala_coarse_residual_interval must be at least one");
  if(E->control.ala_coarse_residual_levels < 1 ||
     E->control.ala_coarse_residual_levels > 10)
      myerror(E, "ala_coarse_residual_levels must be between 1 and 10");
  if(E->control.ala_two_level_offset < 1 ||
     E->control.ala_two_level_offset > 10)
      myerror(E, "ala_two_level_offset must be between 1 and 10");
  if(E->control.ala_two_level_coarse_iterations < 1)
      myerror(E, "ala_two_level_coarse_iterations must be at least one");
  if(E->control.ala_two_level_coarse_damping <= 0.0 ||
     E->control.ala_two_level_coarse_damping >= 2.0/27.0)
      myerror(E, "ala_two_level_coarse_damping must be in (0,2/27)");
  if(strcmp(E->control.ala_two_level_coarse_solver,"jacobi") != 0 &&
     strcmp(E->control.ala_two_level_coarse_solver,"chebyshev") != 0)
      myerror(E, "ala_two_level_coarse_solver must be jacobi or chebyshev");
  if(E->control.ala_two_level_coarse_eigenvalue_min <= 0.0 ||
     E->control.ala_two_level_coarse_eigenvalue_max <=
       E->control.ala_two_level_coarse_eigenvalue_min)
      myerror(E, "ALA two-level Chebyshev eigenvalue interval is invalid");
  if(E->control.ala_two_level_coarse_weight <= 0.0 ||
     E->control.ala_two_level_coarse_weight > 1.0)
      myerror(E, "ala_two_level_coarse_weight must be in (0,1]");
  if(strcmp(E->control.ala_two_level_velocity_solver,"diagonal") != 0 &&
     strcmp(E->control.ala_two_level_velocity_solver,"chebyshev") != 0)
      myerror(E, "ala_two_level_velocity_solver must be diagonal or chebyshev");
  if(E->control.ala_two_level_velocity_iterations < 1)
      myerror(E, "ala_two_level_velocity_iterations must be at least one");
  if(E->control.ala_two_level_velocity_eigenvalue_min <= 0.0 ||
     E->control.ala_two_level_velocity_eigenvalue_max <=
       E->control.ala_two_level_velocity_eigenvalue_min)
      myerror(E, "ALA two-level velocity Chebyshev eigenvalue interval is invalid");
  if(E->control.ala_pressure_multigrid &&
     (E->control.ala_pressure_multigrid_min_level < E->mesh.gridmin ||
      E->control.ala_pressure_multigrid_min_level > E->mesh.levmax))
      myerror(E, "ala_pressure_multigrid_min_level is outside the mesh hierarchy");
  if(E->control.ala_pressure_multigrid &&
     (E->control.ala_pressure_multigrid_pre_smooth < 1 ||
      E->control.ala_pressure_multigrid_post_smooth < 1))
      myerror(E, "ALA pressure multigrid smoothing counts must be positive");
  if(E->control.ala_pressure_multigrid &&
     E->control.ala_pressure_multigrid_coarse_iterations < 1)
      myerror(E, "ALA pressure multigrid coarse iterations must be positive");
  if(E->control.ala_pressure_multigrid &&
     (E->control.ala_pressure_multigrid_damping <= 0.0 ||
      E->control.ala_pressure_multigrid_damping >= 2.0/27.0))
      myerror(E, "ala_pressure_multigrid_damping must be in (0,2/27)");
  if(E->control.ala_pressure_multigrid &&
     (E->control.ala_pressure_multigrid_weight <= 0.0 ||
      E->control.ala_pressure_multigrid_weight > 1.0))
      myerror(E, "ala_pressure_multigrid_weight must be in (0,1]");
  if(E->control.ala_global_coarse_weight <= 0.0 ||
     E->control.ala_global_coarse_weight > 1.0)
      myerror(E, "ala_global_coarse_weight must be in (0,1]");
  if(E->control.ala_global_coarse_regularization < 0.0 ||
     E->control.ala_global_coarse_regularization > 1.0e-4)
      myerror(E, "ala_global_coarse_regularization must be in [0,1e-4]");
  if(E->control.ala_shallow_patch_depth_km <= 0.0)
      myerror(E, "ala_shallow_patch_depth_km must be positive");
  if(E->control.ala_shallow_patch_weight <= 0.0 ||
     E->control.ala_shallow_patch_weight > 1.0)
      myerror(E, "ala_shallow_patch_weight must be in (0,1]");
  if(E->control.ala_shallow_patch_regularization < 0.0 ||
     E->control.ala_shallow_patch_regularization > 0.1)
      myerror(E, "ala_shallow_patch_regularization must be in [0,0.1]");
  if(E->control.ala_shallow_patch_horizontal_elements < 2 ||
     E->control.ala_shallow_patch_horizontal_elements > 8)
      myerror(E, "ala_shallow_patch_horizontal_elements must be in [2,8]");
  if(E->control.ala_shallow_patch_horizontal_stride < 1 ||
     E->control.ala_shallow_patch_horizontal_stride >
       E->control.ala_shallow_patch_horizontal_elements)
      myerror(E, "ala_shallow_patch_horizontal_stride must be in "
              "[1,ala_shallow_patch_horizontal_elements]");
  if(E->control.ala_shallow_patch_mpi_overlap < 1 ||
     E->control.ala_shallow_patch_mpi_overlap > 4)
      myerror(E, "ala_shallow_patch_mpi_overlap must be in [1,4]");
  if(2*E->control.ala_shallow_patch_mpi_overlap >
     E->control.ala_shallow_patch_horizontal_elements)
      myerror(E, "twice ala_shallow_patch_mpi_overlap must not exceed "
              "ala_shallow_patch_horizontal_elements");
  if(strcmp(E->control.ala_shallow_patch_velocity_solver,"diagonal") != 0 &&
     strcmp(E->control.ala_shallow_patch_velocity_solver,"node_block") != 0)
      myerror(E, "ala_shallow_patch_velocity_solver must be diagonal or "
              "node_block");
  if(E->control.ala_geneo_eigenvalue_threshold <= 0.0)
      myerror(E, "ala_geneo_eigenvalue_threshold must be positive");
  if(strcmp(E->control.ala_geneo_basis_type,"spectral") != 0 &&
     strcmp(E->control.ala_geneo_basis_type,"radial_partition") != 0)
      myerror(E, "ala_geneo_basis_type must be spectral or "
              "radial_partition");
  if(strcmp(E->control.ala_geneo_basis_type,"radial_partition") == 0 &&
     (E->control.ala_geneo_min_modes_per_rank
        != E->control.ala_geneo_radial_bins ||
      E->control.ala_geneo_max_modes_per_rank
        != E->control.ala_geneo_radial_bins))
      myerror(E, "radial_partition requires GenEO min/max modes equal "
              "ala_geneo_radial_bins");
  if(strcmp(E->control.ala_geneo_basis_type,"radial_partition") == 0 &&
     E->control.ala_geneo_rank_group_x == 1 &&
     E->control.ala_geneo_rank_group_y == 1)
      myerror(E, "radial_partition requires a cross-rank GenEO group");
  if(E->control.ala_geneo_min_modes_per_rank < 1 ||
     E->control.ala_geneo_max_modes_per_rank <
       E->control.ala_geneo_min_modes_per_rank ||
     E->control.ala_geneo_max_modes_per_rank > 8)
      myerror(E, "ALA GenEO modes per rank must satisfy 1 <= min <= max <= 8");
  if(E->control.ala_geneo_horizontal_bins < 2 ||
     E->control.ala_geneo_horizontal_bins > 8 ||
     E->control.ala_geneo_radial_bins < 1 ||
     E->control.ala_geneo_radial_bins > 4)
      myerror(E, "ALA GenEO bins require horizontal in [2,8] and radial in [1,4]");
  if(E->control.ala_geneo_rank_group_x < 1 ||
     E->control.ala_geneo_rank_group_x > E->parallel.nprocx ||
     E->control.ala_geneo_rank_group_y < 1 ||
     E->control.ala_geneo_rank_group_y > E->parallel.nprocy)
      myerror(E, "ALA GenEO rank groups must fit the horizontal processor grid");
  if(E->control.ala_geneo_horizontal_bins *
       E->control.ala_geneo_rank_group_x *
       E->control.ala_geneo_horizontal_bins *
       E->control.ala_geneo_rank_group_y *
       E->control.ala_geneo_radial_bins > 256)
      myerror(E, "ALA GenEO cross-rank aggregate exceeds 256 spectral bins");
  if(E->control.ala_geneo_max_global_modes < 1 ||
     E->control.ala_geneo_max_global_modes > 4096)
      myerror(E, "ala_geneo_max_global_modes must be between 1 and 4096");
  if(E->control.ala_geneo_weight <= 0.0 ||
     E->control.ala_geneo_weight > 1.0)
      myerror(E, "ala_geneo_weight must be in (0,1]");
  if(E->control.ala_geneo_regularization < 0.0 ||
     E->control.ala_geneo_regularization > 1.0e-3)
      myerror(E, "ala_geneo_regularization must be in [0,1e-3]");
  if(E->control.ala_element_vanka_damping <= 0.0 ||
     E->control.ala_element_vanka_damping > 1.0)
      myerror(E, "ala_element_vanka_damping must be in (0,1]");
  if(E->control.ala_element_vanka_regularization < 0.0 ||
     E->control.ala_element_vanka_regularization > 0.1)
      myerror(E, "ala_element_vanka_regularization must be in [0,0.1]");
  if(E->control.ala_element_vanka_smoother &&
     (strcmp(E->control.SOLVER_TYPE,"multigrid") != 0 ||
      !E->control.ala_pressure_buoyancy ||
      E->control.ala_augmented_lagrangian_gamma <= 0.0))
      myerror(E, "ala_element_vanka_smoother requires multigrid, "
              "compressible_formulation=ala, and positive gamma");

  if(E->control.inv_gruneisen != 0) {
      /* "cg" is the legacy split compressible solver.  Strict ALA may use
       * BiCGStab or its dedicated complete-operator PCG experiment. */
      input_string("uzawa",E->control.uzawa,"cg",m);
      if(strcmp(E->control.uzawa, "cg") == 0) {
          if(E->control.ala_pressure_buoyancy)
              myerror(E,
                      "strict ALA requires uzawa=bicg or uzawa=ala_cg");
          /* more convergence parameters for "cg" */
          input_int("compress_iter_maxstep",&(E->control.compress_iter_maxstep),"100",m);
          input_float("relative_err_accuracy",&(E->control.relative_err_accuracy),"0.001",m);
      }
      else if(strcmp(E->control.uzawa, "bicg") == 0) {
          input_int("compress_iter_maxstep",&(E->control.compress_iter_maxstep),"100",m);
          input_float("relative_err_accuracy",&(E->control.relative_err_accuracy),"0.001",m);
      }
      else if(strcmp(E->control.uzawa, "ala_cg") == 0) {
          if(!E->control.ala_pressure_buoyancy)
              myerror(E,
                      "uzawa=ala_cg requires compressible_formulation=ala");
          input_int("compress_iter_maxstep",&(E->control.compress_iter_maxstep),"100",m);
          input_float("relative_err_accuracy",&(E->control.relative_err_accuracy),"0.001",m);
      }
      else
          myerror(E, "Error: unknown Uzawa iteration\n");
  }
  if((strcmp(E->control.ala_outer_solver,"fgmres") == 0 ||
      strcmp(E->control.ala_outer_solver,"coupled_fgmres") == 0) &&
     strcmp(E->control.uzawa,"ala_cg") != 0)
      myerror(E, "ALA FGMRES outer solvers require uzawa=ala_cg");
  if(E->control.ala_radial_line_preconditioner &&
     (!E->control.precondition ||
      !E->control.ala_pressure_buoyancy ||
      strcmp(E->control.uzawa,"ala_cg") != 0))
      myerror(E,
              "ala_radial_line_preconditioner requires precond=on, "
              "compressible_formulation=ala, and uzawa=ala_cg");
  if(E->control.ala_two_level_preconditioner &&
     (!E->control.precondition ||
      !E->control.ala_pressure_buoyancy ||
      strcmp(E->control.uzawa,"ala_cg") != 0))
      myerror(E,
              "ala_two_level_preconditioner requires precond=on, "
              "compressible_formulation=ala, and uzawa=ala_cg");
  if(E->control.ala_two_level_preconditioner &&
     E->control.ala_radial_line_preconditioner)
      myerror(E, "ALA two-level and radial-line preconditioners are "
              "mutually exclusive");
  if(E->control.ala_global_coarse_preconditioner &&
     !E->control.ala_two_level_preconditioner)
      myerror(E, "ala_global_coarse_preconditioner requires "
              "ala_two_level_preconditioner=on");
  if(E->control.ala_pressure_multigrid &&
     (!E->control.precondition ||
      !E->control.ala_pressure_buoyancy ||
      E->control.ala_augmented_lagrangian_gamma <= 0.0 ||
      (strcmp(E->control.ala_outer_solver,"fgmres") != 0 &&
       strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)))
      myerror(E, "ala_pressure_multigrid requires precond=on, strict ALA, "
              "positive gamma, and an FGMRES outer solver");
  if(E->control.ala_pressure_multigrid &&
     E->control.ala_two_level_preconditioner)
      myerror(E, "ALA pressure multigrid and legacy two-level "
              "preconditioners are mutually exclusive");
  if(E->control.ala_shallow_patch_preconditioner &&
     (!E->control.precondition ||
      !E->control.ala_pressure_buoyancy ||
      strcmp(E->control.uzawa,"ala_cg") != 0))
      myerror(E,
              "ala_shallow_patch_preconditioner requires precond=on, "
              "compressible_formulation=ala, and uzawa=ala_cg");
  if(E->control.ala_shallow_patch_preconditioner &&
     E->control.ala_radial_line_preconditioner)
      myerror(E, "ALA shallow-patch and radial-line preconditioners are "
              "mutually exclusive");
  if(E->control.ala_geneo_preconditioner &&
     (!E->control.precondition ||
      !E->control.ala_pressure_buoyancy ||
      strcmp(E->control.uzawa,"ala_cg") != 0 ||
      !E->control.ala_shallow_patch_preconditioner))
      myerror(E, "ala_geneo_preconditioner requires strict ALA PCG and "
              "ala_shallow_patch_preconditioner=on");
  if(E->control.ala_geneo_preconditioner &&
     E->control.ala_two_level_preconditioner)
      myerror(E, "ALA GenEO and geometric two-level preconditioners are "
              "mutually exclusive");
  if(E->control.ala_feasibility_audit &&
     (!E->control.ala_pressure_buoyancy ||
      strcmp(E->control.uzawa,"ala_cg") != 0))
      myerror(E, "ala_feasibility_audit requires "
              "compressible_formulation=ala and uzawa=ala_cg");

  /* Deprecated compatibility input.  Explicit Ttop/Tbottom are read below. */
  input_float("surfaceT",&(E->control.surface_temp),"-1.0",m);
  /*input_float("adiabaticT0",&(E->control.adiabaticT0),"0.4",m);*/
  input_float("Q0",&(E->control.Q0),"0.0",m);
  /* Q0_enriched gets read in Tracer_setup.c */

  /* data section */
  input_float("gravacc",&(E->data.grav_acc),"9.81",m);
  input_float("thermexp",&(E->data.therm_exp),"3.0e-5",m);
  input_float("cp",&(E->data.Cp),"1200.0",m);
  input_float("thermdiff",&(E->data.therm_diff),"1.0e-6",m);
  input_float("density",&(E->data.density),"3340.0",m);
  input_float("k0",&(E->data.k0),"0.0",m);
  input_float("ks",&(E->data.ks),"0.0",m);
  input_float("rho0",&(E->data.rho0),"0.0",m);
  input_float("Cp0",&(E->data.Cp0),"0.0",m);
  input_float("g0",&(E->data.g0),"0.0",m);
  input_float("alpha0",&(E->data.alpha0),"0.0",m);
  input_float("Ttop",&(E->data.Ttop),"-1.0",m);
  input_float("Tbottom",&(E->data.Tbottom),"-1.0",m);
  input_float("deltaT",&(E->data.ref_temperature),"-1.0",m);
  input_float("density_above",&(E->data.density_above),"1030.0",m);
  input_float("density_below",&(E->data.density_below),"6600.0",m);
  input_float("refvisc",&(E->data.ref_viscosity),"1.0e21",m);

  input_float("radius",&tmp,"6371e3.0",m);
  E->data.radius_km = tmp / 1e3;

  /* The explicit reference system is authoritative when supplied.  A zero
     value selects the historical density/cp/thermdiff input triplet. */
  if(E->data.rho0 <= 0.0)
      E->data.rho0 = E->data.density;
  if(E->data.Cp0 <= 0.0)
      E->data.Cp0 = E->data.Cp;
  if(E->data.g0 <= 0.0)
      E->data.g0 = E->data.grav_acc;
  if(E->data.alpha0 <= 0.0)
      E->data.alpha0 = E->data.therm_exp;
  if(E->data.k0 <= 0.0)
      E->data.k0 = E->data.therm_diff * E->data.rho0 * E->data.Cp0;
  /* Legacy files did not distinguish the model normalization from k0. */
  if(E->data.ks <= 0.0)
      E->data.ks = E->data.k0;

  if(!isfinite(E->data.ks) || !isfinite(E->data.k0) ||
     !isfinite(E->data.rho0) || !isfinite(E->data.Cp0) ||
     !isfinite(E->data.g0) || !isfinite(E->data.alpha0) ||
     E->data.ks <= 0.0 || E->data.k0 <= 0.0 ||
     E->data.rho0 <= 0.0 || E->data.Cp0 <= 0.0 ||
     E->data.g0 <= 0.0 || E->data.alpha0 <= 0.0)
      myerror(E, "ks, k0, rho0, Cp0, g0, and alpha0 must be finite and positive");

  E->data.kappa0 = E->data.k0 / (E->data.rho0 * E->data.Cp0);
  E->control.reference_conductivity = E->data.ks / E->data.k0;

  /* Preserve old call sites and old input files without retaining two
     independent thermal scales. */
  E->data.therm_cond = E->data.k0;
  E->data.therm_diff = E->data.kappa0;
  E->data.density = E->data.rho0;
  E->data.Cp = E->data.Cp0;
  E->data.grav_acc = E->data.g0;
  E->data.therm_exp = E->data.alpha0;

  if((E->data.Ttop >= 0.0) != (E->data.Tbottom >= 0.0))
      myerror(E, "Ttop and Tbottom must be specified together");
  if(E->data.Ttop >= 0.0) {
      if(!isfinite(E->data.Ttop) || !isfinite(E->data.Tbottom) ||
         E->data.Tbottom <= E->data.Ttop)
          myerror(E, "Ttop and Tbottom must be finite with Tbottom > Ttop >= 0");
      if(E->parallel.me == 0 && E->data.ref_temperature > 0.0)
          fprintf(stderr, "WARNING: deltaT is deprecated and ignored; using Tbottom-Ttop\n");
      if(E->parallel.me == 0 && E->control.surface_temp >= 0.0)
          fprintf(stderr, "WARNING: surfaceT is deprecated and ignored; using Ttop/DeltaT\n");
      E->data.ref_temperature = E->data.Tbottom - E->data.Ttop;
      E->control.surface_temp = E->data.Ttop / E->data.ref_temperature;
  }
  else {
      if(E->parallel.me == 0)
          fprintf(stderr, "WARNING: Ttop/Tbottom absent; deltaT/surfaceT compatibility mode is deprecated\n");
      if(E->control.surface_temp < 0.0)
          E->control.surface_temp = 0.1;
  }

  if(E->data.ref_temperature > 0.0)
      E->control.Atemp = E->data.rho0 * E->data.g0 * E->data.alpha0
        * E->data.ref_temperature
        * (E->data.radius_km * E->data.radius_km * E->data.radius_km * 1e9)
        / (E->data.ref_viscosity * E->data.kappa0);
  else
      E->data.ref_temperature = E->control.Atemp * E->data.kappa0
        * E->data.ref_viscosity
        / (E->data.rho0 * E->data.g0 * E->data.alpha0)
        / (E->data.radius_km * E->data.radius_km * E->data.radius_km * 1e9);

  if(E->data.Ttop < 0.0) {
      E->data.Ttop = E->control.surface_temp * E->data.ref_temperature;
      E->data.Tbottom = E->data.Ttop + E->data.ref_temperature;
  }
  if(E->parallel.me == 0)
      fprintf(stderr,
              "Temperature reference:\nTtop = %.9g K\nTbottom = %.9g K\nDeltaT = %.9g K\n",
              E->data.Ttop, E->data.Tbottom, E->data.ref_temperature);

  output_common_input(E);
  h5input_params(E);
  phase_change_input(E);
  lith_age_input(E);

  tic_input(E);
  tracer_input(E);

  viscosity_input(E);		/* moved the viscosity input behind
				   the tracer input */

  (E->problem_settings)(E);


  return;
}


/* ===================================
   Functions which set up details
   common to all problems follow ...
   ===================================  */

void allocate_common_vars(E)
     struct All_variables *E;

{
    void set_up_nonmg_aliases();
    int m,n,snel,nsf,elx,ely,nox,noy,noz,nno,nel,npno;
    int k,i,j,d,l,nno_l,npno_l,nozl,nnov_l,nxyz;

    m=0;
    n=1;

 for (j=1;j<=E->sphere.caps_per_proc;j++)  {

  npno = E->lmesh.npno;
  nel  = E->lmesh.nel;
  nno  = E->lmesh.nno;
  nsf  = E->lmesh.nsf;
  noz  = E->lmesh.noz;
  nox  = E->lmesh.nox;
  noy  = E->lmesh.noy;
  elx  = E->lmesh.elx;
  ely  = E->lmesh.ely;

  E->P[j]        = (double *) malloc((npno+1)*sizeof(double));
  E->T[j]        = (double *) malloc((nno+1)*sizeof(double));
  E->NP[j]       = (float *) malloc((nno+1)*sizeof(float));
  E->edot[j]     = (float *) malloc((nno+1)*sizeof(float));
  E->buoyancy[j] = (double *) malloc((nno+1)*sizeof(double));

  E->gstress[j] = (float *) malloc((6*nno+1)*sizeof(float));
  E->stress[j]   = (float *) malloc((12*nsf+1)*sizeof(float));

  for(i=1;i<=E->mesh.nsd;i++)
      E->sphere.cap[j].TB[i] = (float *)  malloc((nno+1)*sizeof(float));

  E->age[j]      = (float *)malloc((nsf+2)*sizeof(float));

  E->slice.tpg[j]      = (float *)malloc((nsf+2)*sizeof(float));
  E->slice.tpgb[j]     = (float *)malloc((nsf+2)*sizeof(float));
  E->slice.divg[j]     = (float *)malloc((nsf+2)*sizeof(float));
  E->slice.vort[j]     = (float *)malloc((nsf+2)*sizeof(float));
  E->slice.shflux[j]      = (float *)malloc((nsf+2)*sizeof(float));
  E->slice.bhflux[j]      = (float *)malloc((nsf+2)*sizeof(float));
  E->slice.shflux_CBF[j]  = (float *)malloc((nsf+2)*sizeof(float));
  E->slice.bhflux_CBF[j]  = (float *)malloc((nsf+2)*sizeof(float));
  /*  if(E->mesh.topvbc==2 && E->control.pseudo_free_surf) */
  E->slice.freesurf[j]    = (float *)malloc((nsf+2)*sizeof(float));

  E->mat[j] = (int *) malloc((nel+2)*sizeof(int));
  E->VIP[j] = (float *) malloc((nel+2)*sizeof(float));

  E->heating_adi[j]    = (double *) malloc((nel+1)*sizeof(double));
  E->heating_adi_base[j] = (double *) malloc((nel+1)*sizeof(double));
  E->heating_phase_adi[j] = (double *) malloc((nel+1)*sizeof(double));
  E->heating_visc[j]   = (double *) malloc((nel+1)*sizeof(double));
  E->heating_latent[j] = (double *) malloc((nel+1)*sizeof(double));
  E->heating_internal[j] = (double *) malloc((nel+1)*sizeof(double));
  E->heating_assim[j] = (double *) malloc((nel+1)*sizeof(double));
  E->assim_delta_T[j] = (double *) malloc((nno+1)*sizeof(double));

  /* lump mass matrix for the energy eqn */
  E->TMass[j] = (double *) malloc((nno+1)*sizeof(double));

  nxyz = max(nox*noz,nox*noy);
  nxyz = 2*max(nxyz,noz*noy);

  E->sien[j]         = (struct SIEN *) malloc((nxyz+2)*sizeof(struct SIEN));
  E->surf_element[j] = (int *) malloc((nxyz+2)*sizeof(int));
  E->surf_node[j]    = (int *) malloc((nsf+2)*sizeof(int));

  }         /* end for cap j  */

  /* density field */
  E->rho      = (double *) malloc((nno+1)*sizeof(double));

  /* horizontal average */
  E->Have.T         = (float *)malloc((E->lmesh.noz+2)*sizeof(float));
  E->Have.V[1]      = (float *)malloc((E->lmesh.noz+2)*sizeof(float));
  E->Have.V[2]      = (float *)malloc((E->lmesh.noz+2)*sizeof(float));

 for(i=E->mesh.levmin;i<=E->mesh.levmax;i++) {
  E->sphere.R[i] = (double *)  malloc((E->lmesh.NOZ[i]+1)*sizeof(double));
  for (j=1;j<=E->sphere.caps_per_proc;j++)  {
    nno  = E->lmesh.NNO[i];
    npno = E->lmesh.NPNO[i];
    nel  = E->lmesh.NEL[i];
    nox = E->lmesh.NOX[i];
    noz = E->lmesh.NOZ[i];
    noy = E->lmesh.NOY[i];
    elx = E->lmesh.ELX[i];
    ely = E->lmesh.ELY[i];
    snel=E->lmesh.SNEL[i];

    for(d=1;d<=E->mesh.nsd;d++)   {
      E->X[i][j][d]  = (double *)  malloc((nno+1)*sizeof(double));
      E->SX[i][j][d]  = (double *)  malloc((nno+1)*sizeof(double));
      }

    for(d=0;d<=3;d++)
      E->SinCos[i][j][d]  = (double *)  malloc((nno+1)*sizeof(double));

    E->IEN[i][j] = (struct IEN *)   malloc((nel+2)*sizeof(struct IEN));
    E->EL[i][j]  = (struct SUBEL *) malloc((nel+2)*sizeof(struct SUBEL));
    E->sphere.area1[i][j] = (double *) malloc((snel+1)*sizeof(double));
    for (k=1;k<=4;k++)
      E->sphere.angle1[i][j][k] = (double *) malloc((snel+1)*sizeof(double));

    E->MASS[i][j]     = (float *) malloc((nno+1)*sizeof(float));
    E->ECO[i][j] = (struct COORD *) malloc((nno+2)*sizeof(struct COORD));

    E->TWW[i][j] = (struct FNODE *)   malloc((nel+2)*sizeof(struct FNODE));

    for(d=1;d<=E->mesh.nsd;d++)
      for(l=1;l<=E->lmesh.NNO[i];l++)  {
        E->SX[i][j][d][l] = 0.0;
        E->X[i][j][d][l] = 0.0;
        }

    }
  }

 for(i=0;i<=E->output.llmax;i++)
  E->sphere.hindex[i] = (int *) malloc((E->output.llmax+3)
				       *sizeof(int));

 for(i=E->mesh.gridmin;i<=E->mesh.gridmax;i++)
  for (j=1;j<=E->sphere.caps_per_proc;j++)  {

    nno  = E->lmesh.NNO[i];
    npno = E->lmesh.NPNO[i];
    nel  = E->lmesh.NEL[i];
    nox = E->lmesh.NOX[i];
    noz = E->lmesh.NOZ[i];
    noy = E->lmesh.NOY[i];
    elx = E->lmesh.ELX[i];
    ely = E->lmesh.ELY[i];

    nxyz = elx*ely;
    E->CC[i][j] =(struct CC *)  malloc((1)*sizeof(struct CC));
    E->CCX[i][j]=(struct CCX *)  malloc((1)*sizeof(struct CCX));
    /* Test */
    E->ELEMENT[i][j] = (unsigned int *) malloc ((nel+1)*sizeof(unsigned int));

    for (k=1;k<=nel;k++)
       E->ELEMENT[i][j][k] = 0;
    /*ccccc*/

    E->elt_del[i][j] = (struct EG *) malloc((nel+1)*sizeof(struct EG));

    if(E->control.inv_gruneisen != 0)
        E->elt_c[i][j] = (struct EC *) malloc((nel+1)*sizeof(struct EC));

    E->EVI[i][j] = (float *) malloc((nel+1)*vpoints[E->mesh.nsd]*sizeof(float));
    E->BPI[i][j] = (double *) malloc((npno+1)*sizeof(double));
    E->ALA_BPI_line_diag[i][j] =
        (double *) malloc((npno+1)*sizeof(double));
    E->ALA_BPI_line_lower[i][j] =
        (double *) malloc((npno+1)*sizeof(double));
    E->ALA_BPI_line_valid[i][j] =
        (unsigned char *) malloc((elx*ely+1)*sizeof(unsigned char));

    E->ID[i][j]  = (struct ID *)    malloc((nno+1)*sizeof(struct ID));
    E->VI[i][j]  = (float *)        malloc((nno+1)*sizeof(float));
    E->NODE[i][j] = (unsigned int *)malloc((nno+1)*sizeof(unsigned int));

    /* Lijun add: a time-dependent viscosity component */
    E->oldEVI[i][j] = (float *) malloc((nel+1)*vpoints[E->mesh.nsd]*sizeof(float));
    E->oldVI[i][j]  = (float *)        malloc((nno+1)*sizeof(float));
    E->long_t[j] = (float *) malloc((nno+1)*sizeof(float));
    /* end of add */

    nxyz = max(nox*noz,nox*noy);
    nxyz = 2*max(nxyz,noz*noy);
    nozl = max(noy,nox*2);



    E->parallel.EXCHANGE_sNODE[i][j] = (struct PASS *) malloc((nozl+2)*sizeof(struct PASS));
    E->parallel.NODE[i][j]   = (struct BOUND *) malloc((nxyz+2)*sizeof(struct BOUND));
    E->parallel.EXCHANGE_NODE[i][j]= (struct PASS *) malloc((nxyz+2)*sizeof(struct PASS));
    E->parallel.EXCHANGE_ID[i][j] = (struct PASS *) malloc((nxyz*E->mesh.nsd+3)*sizeof(struct PASS));

    for(l=1;l<=E->lmesh.NNO[i];l++)  {
      E->NODE[i][j][l] = (INTX | INTY | INTZ);  /* and any others ... */
      E->VI[i][j][l] = 1.0;
      }


    }         /* end for cap and i & j  */


 for (j=1;j<=E->sphere.caps_per_proc;j++)  {

  for(k=1;k<=E->mesh.nsd;k++)
    for(i=1;i<=E->lmesh.nno;i++)
      E->sphere.cap[j].TB[k][i] = 0.0;

  for(i=1;i<=E->lmesh.nno;i++) {
     E->T[j][i] = 0.0;
     E->assim_delta_T[j][i] = 0.0;
  }

  for(i=1;i<=E->lmesh.nel;i++)   {
      E->mat[j][i]=1;
      E->VIP[j][i]=1.0;

      E->heating_adi[j][i] = 0;
      E->heating_adi_base[j][i] = 0;
      E->heating_phase_adi[j][i] = 0;
      E->heating_visc[j][i] = 0;
      E->heating_latent[j][i] = 1.0;
      E->heating_internal[j][i] = 0;
      E->heating_assim[j][i] = 0;
  }

  for(i=1;i<=E->lmesh.npno;i++)
      E->P[j][i] = 0.0;

  mat_prop_allocate(E);
  phase_change_allocate(E);
  set_up_nonmg_aliases(E,j);

  }         /* end for cap j  */

  if (strcmp(E->output.format, "hdf5") == 0)
      h5output_allocate_memory(E);

  return;
  }

/*  =========================================================  */

void allocate_velocity_vars(E)
     struct All_variables *E;

{
    int m,n,i,j,k,l;
 m=0;
 n=1;
  for (j=1;j<=E->sphere.caps_per_proc;j++)   {
    E->lmesh.nnov = E->lmesh.nno;
    E->lmesh.neq = E->lmesh.nnov * E->mesh.nsd;

    E->temp[j] = (double *) malloc((E->lmesh.neq+1)*sizeof(double));
    E->temp1[j] = (double *) malloc(E->lmesh.neq*sizeof(double));
    E->F[j] = (double *) malloc(E->lmesh.neq*sizeof(double));
    E->U[j] = (double *) malloc((E->lmesh.neq+1)*sizeof(double));
    E->u1[j] = (double *) malloc((E->lmesh.neq+1)*sizeof(double));


    for(i=1;i<=E->mesh.nsd;i++) {
      E->sphere.cap[j].V[i] = (float *) malloc((E->lmesh.nnov+1)*sizeof(float));
      E->sphere.cap[j].VB[i] = (float *)malloc((E->lmesh.nnov+1)*sizeof(float));
      E->sphere.cap[j].Vprev[i] = (float *) malloc((E->lmesh.nnov+1)*sizeof(float));
    }

    for(i=0;i<E->lmesh.neq;i++)
      E->U[j][i] = E->temp[j][i] = E->temp1[j][i] = 0.0;

    /* the following lines are causing a memory outflow */
    /*if(E->control.tracer == 1)  {
      for(i=1;i<=E->mesh.nsd;i++)     {
        E->GV[j][i]=(float*) malloc(((E->lmesh.nno+1)*E->parallel.nproc+1)*sizeof(float));
        E->GV1[j][i]=(float*) malloc(((E->lmesh.nno+1)*E->parallel.nproc+1)*sizeof(float));
        E->V[j][i]=(float*) malloc((E->lmesh.nno+1)*sizeof(float));

        for(k=0;k<(E->lmesh.nno+1)*E->parallel.nproc;k++)   {
          E->GV[j][i][k]=0.0;
          E->GV1[j][i][k]=0.0;
        }
      }
    }*/

    for(k=1;k<=E->mesh.nsd;k++)
      for(i=1;i<=E->lmesh.nnov;i++)
        E->sphere.cap[j].VB[k][i] = 0.0;

  }       /* end for cap j */

  for(l=E->mesh.gridmin;l<=E->mesh.gridmax;l++)
    for (j=1;j<=E->sphere.caps_per_proc;j++)   {
      E->lmesh.NEQ[l] = E->lmesh.NNOV[l] * E->mesh.nsd;

      E->BI[l][j] = (double *) malloc((E->lmesh.NEQ[l])*sizeof(double));
      E->ALA_velocity_BI[l][j] =
          (double *) malloc((E->lmesh.NEQ[l])*sizeof(double));
      if(E->control.ala_element_vanka_smoother) {
        E->ALA_vanka_overlap_BI[l][j] =
            (double *) malloc((E->lmesh.NEQ[l])*sizeof(double));
        E->ALA_vanka_chol[l][j] = (higher_precision *)malloc(
            (E->lmesh.NEL[l]+1)*ALA_VANKA_CHOL_SIZE
            *sizeof(higher_precision));
        E->ALA_vanka_valid[l][j] = (unsigned char *)malloc(
            (E->lmesh.NEL[l]+1)*sizeof(unsigned char));
        if(E->ALA_vanka_overlap_BI[l][j]==NULL ||
           E->ALA_vanka_chol[l][j]==NULL ||
           E->ALA_vanka_valid[l][j]==NULL)
          myerror(E,"Unable to allocate full ALA element-Vanka cache");
      }
      else {
        E->ALA_vanka_overlap_BI[l][j] = NULL;
        E->ALA_vanka_chol[l][j] = NULL;
        E->ALA_vanka_valid[l][j] = NULL;
      }
      k = (E->lmesh.NOX[l]*E->lmesh.NOZ[l]+E->lmesh.NOX[l]*E->lmesh.NOY[l]+
          E->lmesh.NOY[l]*E->lmesh.NOZ[l])*6;
      E->zero_resid[l][j] = (int *) malloc((k+2)*sizeof(int));
      E->parallel.Skip_id[l][j] = (int *) malloc((k+2)*sizeof(int));

      for(i=0;i<E->lmesh.NEQ[l];i++) {
         E->BI[l][j][i]=0.0;
         E->ALA_velocity_BI[l][j][i]=0.0;
         if(E->control.ala_element_vanka_smoother) {
           E->ALA_vanka_overlap_BI[l][j][i]=0.0;
         }
         }

      }   /* end for j & l */

  return;
 }


/*  =========================================================  */

void global_default_values(E)
     struct All_variables *E;
{

  /* FIRST: values which are not changed routinely by the user */

  E->control.v_steps_low = 10;
  E->control.v_steps_upper = 1;
  E->control.accuracy = 1.0e-6;
  E->control.vaccuracy = 1.0e-8;
  E->control.verbose=0; /* debugging/profiles */

  /* SECOND: values for which an obvious default setting is useful */

    E->control.stokes=0;
    E->control.restart=0;
    E->control.CONVECTION = 0;
    E->control.CART2D = 0;
    E->control.CART3D = 0;
    E->control.CART2pt5D = 0;
    E->control.AXI = 0;
    E->control.CONJ_GRAD = 0;
    E->control.NMULTIGRID = 0;
    E->control.EMULTIGRID = 0;
    E->control.augmented_Lagr = 0;
    E->control.augmented = 0.0;
    E->control.ala_augmented_lagrangian_gamma = 0.0;
    E->control.ala_schur_symmetry_check = 0;
    E->control.ala_schur_symmetry_tolerance = 1.0e-3;
    strcpy(E->control.ala_beta_element_source,"supplied_average");
    E->control.ala_inner_accuracy_max = 1.0e-4;
    E->control.ala_inner_accuracy_factor = 1.0e-2;
    E->control.ala_pcg_restart_interval = 20;
    strcpy(E->control.ala_outer_solver, "pcg");
    E->control.ala_unaugmented_momentum_tolerance = 0.0;
    E->control.ala_feasibility_audit = 0;
    E->control.ala_feasibility_window = 20;
    E->control.ala_feasibility_min_reduction = 0.02;
    E->control.ala_hybrid_convergence = 0;
    E->control.ala_div_v_tolerance = 1.0e-7;
    E->control.ala_update_tolerance = 1.0e-3;
    E->control.ala_consecutive_steps = 3;
    E->control.ala_depth_diagnostics = 0;
    E->control.ala_depth_diagnostic_interval = 5;
    E->control.ala_depth_diagnostic_bins = 8;
    E->control.ala_beta_causal_diagnostics = 0;
    E->control.ala_coarse_residual_diagnostics = 0;
    E->control.ala_coarse_residual_interval = 5;
    E->control.ala_coarse_residual_levels = 3;
    E->control.ala_two_level_preconditioner = 0;
    E->control.ala_two_level_offset = 2;
    E->control.ala_two_level_coarse_iterations = 12;
    E->control.ala_two_level_coarse_damping = 0.03;
    strcpy(E->control.ala_two_level_coarse_solver,"chebyshev");
    E->control.ala_two_level_coarse_eigenvalue_min = 0.01;
    E->control.ala_two_level_coarse_eigenvalue_max = 27.0;
    E->control.ala_two_level_coarse_weight = 0.05;
    strcpy(E->control.ala_two_level_velocity_solver,"chebyshev");
    E->control.ala_two_level_velocity_iterations = 16;
    E->control.ala_two_level_velocity_eigenvalue_min = 0.01;
    E->control.ala_two_level_velocity_eigenvalue_max = 4.0;
    E->control.ala_pressure_multigrid = 0;
    E->control.ala_pressure_multigrid_min_level = 0;
    E->control.ala_pressure_multigrid_pre_smooth = 2;
    E->control.ala_pressure_multigrid_post_smooth = 2;
    E->control.ala_pressure_multigrid_coarse_iterations = 20;
    E->control.ala_pressure_multigrid_damping = 0.04;
    E->control.ala_pressure_multigrid_weight = 0.5;
    E->control.ala_global_coarse_preconditioner = 0;
    E->control.ala_global_coarse_weight = 1.0;
    E->control.ala_global_coarse_regularization = 1.0e-10;
    E->control.ala_shallow_patch_preconditioner = 0;
    E->control.ala_shallow_patch_depth_km = 410.0;
    E->control.ala_shallow_patch_weight = 0.25;
    E->control.ala_shallow_patch_regularization = 1.0e-3;
    E->control.ala_shallow_patch_horizontal_elements = 4;
    E->control.ala_shallow_patch_horizontal_stride = 2;
    E->control.ala_shallow_patch_mpi_overlap = 2;
    strcpy(E->control.ala_shallow_patch_velocity_solver,"diagonal");
    E->control.ala_geneo_preconditioner = 0;
    strcpy(E->control.ala_geneo_basis_type,"spectral");
    E->control.ala_geneo_eigenvalue_threshold = 0.20;
    E->control.ala_geneo_min_modes_per_rank = 1;
    E->control.ala_geneo_max_modes_per_rank = 2;
    E->control.ala_geneo_horizontal_bins = 4;
    E->control.ala_geneo_radial_bins = 2;
    E->control.ala_geneo_rank_group_x = 1;
    E->control.ala_geneo_rank_group_y = 1;
    E->control.ala_geneo_max_global_modes = 400;
    E->control.ala_geneo_weight = 1.0;
    E->control.ala_geneo_regularization = 1.0e-8;
    E->control.ala_radial_line_preconditioner = 0;
    E->control.ala_element_vanka_smoother = 0;
    E->control.ala_element_vanka_damping = 0.8;
    E->control.ala_element_vanka_regularization = 1.0e-8;
    strcpy(E->control.log_template,"datafile");

    E->control.GRID_TYPE=1;

    E->trace.fpt = NULL;
    E->control.tracer = 0;
    E->composition.on = 0;

  E->parallel.nprocx=1; E->parallel.nprocz=1; E->parallel.nprocy=1;

  E->mesh.levmax=0;
  E->mesh.levmin=0;
  E->mesh.gridmax=0;
  E->mesh.gridmin=0;
  E->mesh.noz = 1;    E->mesh.nzs = 1;  E->lmesh.noz = 1;    E->lmesh.nzs = 1;
  E->mesh.noy = 1;    E->mesh.nys = 1;  E->lmesh.noy = 1;    E->lmesh.nys = 1;
  E->mesh.nox = 1;    E->mesh.nxs = 1;  E->lmesh.nox = 1;    E->lmesh.nxs = 1;

  E->sphere.ro = 1.0;
  E->sphere.ri = 0.5;

  E->control.precondition = 0;  /* for larger visc contrasts turn this back on  */

  E->mesh.toptbc = 1; /* fixed t */
  E->mesh.bottbc = 1;
  E->mesh.topvbc = 0; /* stress */
  E->mesh.botvbc = 0;
  E->control.VBXtopval=0.0;
  E->control.VBYtopval=0.0;
  E->control.VBXbotval=0.0;
  E->control.VBYbotval=0.0;

  E->data.radius_km = 6370.0; /* Earth, whole mantle defaults */
  E->data.grav_acc = 9.81;
  E->data.therm_diff = 1.0e-6;
  E->data.therm_exp = 3.e-5;
  E->data.density = 3300.0;
  E->data.ref_viscosity=1.e21;
  E->data.density_above = 1000.0;    /* sea water */
  E->data.density_below = 6600.0;    /* sea water */

  E->data.Cp = 1200.0;
  E->data.ks = 0.0;
  E->data.rho0 = E->data.density;
  E->data.Cp0 = E->data.Cp;
  E->data.g0 = E->data.grav_acc;
  E->data.alpha0 = E->data.therm_exp;
  E->data.k0 = E->data.therm_diff * E->data.rho0 * E->data.Cp0;
  E->data.kappa0 = E->data.therm_diff;
  E->data.therm_cond = E->data.k0;
  E->control.reference_conductivity = 1.0;
  E->control.requested_reference_conductivity = -1.0;
  E->data.res_density = 3300.0;  /* density when X = ... */
  E->data.res_density_X = 0.3;
  E->data.melt_density = 2800.0;
  E->data.permeability = 3.0e-10;
  E->data.gas_const = 8.3;
  E->data.surf_heat_flux = 4.4e-2;
  E->data.grav_const = 6.673e-11;
  E->data.youngs_mod = 1.0e11;
  E->data.Te = 0.0;
  E->data.T_sol0 = 1373.0;      /* Dave's values 1991 (for the earth) */
  E->data.Tsurf = 273.0;
  E->data.dTsol_dz = 3.4e-3 ;
  E->data.dTsol_dF = 440.0;
  E->data.dT_dz = 0.48e-3;
  E->data.delta_S = 250.0;
  E->data.ref_temperature = 2 * 1350.0; /* fixed temperature ... delta T */
  E->data.Ttop = 0.1 * E->data.ref_temperature;
  E->data.Tbottom = E->data.Ttop + E->data.ref_temperature;
  E->control.surface_temp = 0.1;

  /* THIRD: you forgot and then went home, let's see if we can help out */

    sprintf(E->control.data_prefix,"citcom.tmp.%d",getpid());

    E->control.NASSEMBLE = 0;

    E->monitor.elapsed_time=0.0;

    E->control.record_all_until = 10000000;

  return;  }


/* =============================================================
   ============================================================= */

void check_bc_consistency(E)
     struct All_variables *E;

{ int i,j,lev;

  for (j=1;j<=E->sphere.caps_per_proc;j++)  {
    for(i=1;i<=E->lmesh.nno;i++)    {
      if ((E->node[j][i] & VBX) && (E->node[j][i] & SBX))
        printf("Inconsistent x velocity bc at %d\n",i);
      if ((E->node[j][i] & VBZ) && (E->node[j][i] & SBZ))
        printf("Inconsistent z velocity bc at %d\n",i);
      if ((E->node[j][i] & VBY) && (E->node[j][i] & SBY))
        printf("Inconsistent y velocity bc at %d\n",i);
      if ((E->node[j][i] & TBX) && (E->node[j][i] & FBX))
        printf("Inconsistent x temperature bc at %d\n",i);
      if ((E->node[j][i] & TBZ) && (E->node[j][i] & FBZ))
        printf("Inconsistent z temperature bc at %d\n",i);
      if ((E->node[j][i] & TBY) && (E->node[j][i] & FBY))
        printf("Inconsistent y temperature bc at %d\n",i);
      }
    }          /* end for j */

  for(lev=E->mesh.gridmin;lev<=E->mesh.gridmax;lev++)
    for (j=1;j<=E->sphere.caps_per_proc;j++)  {
      for(i=1;i<=E->lmesh.NNO[lev];i++)        {
        if ((E->NODE[lev][j][i] & VBX) && (E->NODE[lev][j][i]  & SBX))
          printf("Inconsistent x velocity bc at %d,%d\n",lev,i);
        if ((E->NODE[lev][j][i] & VBZ) && (E->NODE[lev][j][i]  & SBZ))
          printf("Inconsistent z velocity bc at %d,%d\n",lev,i);
        if ((E->NODE[lev][j][i] & VBY) && (E->NODE[lev][j][i]  & SBY))
          printf("Inconsistent y velocity bc at %d,%d\n",lev,i);
        /* Tbc's not applicable below top level */
        }

    }   /* end for  j and lev */

  return;

}

void set_up_nonmg_aliases(E,j)
     struct All_variables *E;
     int j;

{ /* Aliases for functions only interested in the highest mg level */

  int i;

  E->eco[j] = E->ECO[E->mesh.levmax][j];
  E->ien[j] = E->IEN[E->mesh.levmax][j];
  E->id[j] = E->ID[E->mesh.levmax][j];
  E->Vi[j] = E->VI[E->mesh.levmax][j];
  E->EVi[j] = E->EVI[E->mesh.levmax][j];
  E->node[j] = E->NODE[E->mesh.levmax][j];
  E->cc[j] = E->CC[E->mesh.levmax][j];
  E->ccx[j] = E->CCX[E->mesh.levmax][j];
  E->Mass[j] = E->MASS[E->mesh.levmax][j];
  E->element[j] = E->ELEMENT[E->mesh.levmax][j];

  for (i=1;i<=E->mesh.nsd;i++)    {
    E->x[j][i] = E->X[E->mesh.levmax][j][i];
    E->sx[j][i] = E->SX[E->mesh.levmax][j][i];
    }

  return; }

void report(E,string)
     struct All_variables *E;
     char * string;
{ if(E->control.verbose && E->parallel.me==0)
    { fprintf(stderr,"%s\n",string);
      fflush(stderr);
    }
  return;
}

void record(E,string)
     struct All_variables *E;
     char * string;
{ if(E->control.verbose && E->fp)
    { fprintf(E->fp,"%s\n",string);
      fflush(E->fp);
    }

  return;
}



/* =============================================================
   Initialize values which are not problem dependent.
   NOTE: viscosity may be a function of all previous
   input fields (temperature, pressure, velocity, chemistry) and
   so is always to be done last.
   ============================================================= */


/* This function is replaced by CitcomS.Components.IC.launch()*/
void common_initial_fields(E)
    struct All_variables *E;
{
    void initial_pressure();
    void initial_velocity();
    /*void read_viscosity_option();*/
    void initial_viscosity();

    report(E,"Initialize pressure field");
    initial_pressure(E);
    report(E,"Initialize velocity field");
    initial_velocity(E);
    report(E,"Initialize viscosity field");
    /*get_viscosity_option(E);*/
    initial_viscosity(E);

    return;

}

/* ========================================== */

void initial_pressure(E)
     struct All_variables *E;
{
    int i,m;

  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=E->lmesh.npno;i++)
      E->P[m][i]=0.0;

  return;
}

void initial_velocity(E)
     struct All_variables *E;
{
    int i,m;

  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=E->lmesh.nnov;i++)   {
        E->sphere.cap[m].V[1][i]=0.0;
        E->sphere.cap[m].V[2][i]=0.0;
        E->sphere.cap[m].V[3][i]=0.0;
        }

    return;
}



static void log_timestamp(struct All_variables *E, char *timestamp,
                          size_t timestamp_size)
{
  time_t now;
  struct tm *utc;

  memset(timestamp,0,timestamp_size);
  if(E->parallel.me == 0) {
    now=time(NULL);
    utc=gmtime(&now);
    if(utc == NULL ||
       strftime(timestamp,timestamp_size,"%Y%m%dT%H%M%SZ",utc) == 0)
      snprintf(timestamp,timestamp_size,"unknown-time");
  }
  MPI_Bcast(timestamp,(int)timestamp_size,MPI_CHAR,0,E->parallel.world);
}


static void open_log(struct All_variables *E)
{
  char logfile[512];
  char basename[256];
  char pressure_iterations[32];
  char timestamp[32];
  int written;

  E->fp = NULL;
  if(strcmp(E->control.log_template,"datafile") == 0) {
    if (strcmp(E->output.format, "ascii-gz") == 0)
      snprintf(logfile,sizeof(logfile),"%s/log",E->control.data_dir);
    else
      snprintf(logfile,sizeof(logfile),"%s.log",E->control.data_file);
  }
  else {
    snprintf(basename,sizeof(basename),"%s",E->control.log_template);
    snprintf(pressure_iterations,sizeof(pressure_iterations),"%d",
             E->control.p_iterations);
    log_timestamp(E,timestamp,sizeof(timestamp));
    expand_str(basename,sizeof(basename),"%P",pressure_iterations);
    expand_str(basename,sizeof(basename),"%T",timestamp);
    written = snprintf(logfile,sizeof(logfile),"%s/%s.log",E->control.data_dir,
                       basename);
    if(written < 0 || (size_t)written >= sizeof(logfile)) {
      fprintf(stderr,"error: expanded logfile path is too long\n");
      parallel_process_termination();
    }
  }

  if (E->control.restart || E->control.post_p)
      /* append the log file if restart */
      E->fp = output_open(logfile, "a");
  else
      E->fp = output_open(logfile, "w");

  fprintf(E->fp,"RUN_LOG file=%s template=%s piterations=%d\n",
          logfile,E->control.log_template,E->control.p_iterations);
  fflush(E->fp);
  if(E->parallel.me == 0)
    fprintf(stderr,"RUN_LOG file=%s piterations=%d\n",logfile,
            E->control.p_iterations);

  return;
}


static void open_time(struct All_variables *E)
{
  char timeoutput[255];

  E->fptime = NULL;
  if (E->parallel.me == 0) {
  if (strcmp(E->output.format, "ascii-gz") == 0)
    sprintf(timeoutput,"%s/time", E->control.data_dir);
  else
    sprintf(timeoutput,"%s.time", E->control.data_file);

  if (E->control.restart || E->control.post_p)
      /* append the time file if restart */
      E->fptime = output_open(timeoutput, "a");
  else
      E->fptime = output_open(timeoutput, "w");
  }

  return;
}


static void open_info(struct All_variables *E)
{
  char output_file[255];

  E->fp_out = NULL;
  if (E->control.verbose) {
  if (strcmp(E->output.format, "ascii-gz") == 0)
    sprintf(output_file,"%s/info.%d", E->control.data_dir, E->parallel.me);
  else
    sprintf(output_file,"%s.info.%d", E->control.data_file, E->parallel.me);
  E->fp_out = output_open(output_file, "w");
  }

  return;
}

void open_qfiles(struct All_variables *E) /* additional heat
					     flux output */
{
  char output_file[255];

  /* only one CPU will write to those */
  if((E->parallel.me_loc[3] == E->parallel.nprocz-1) &&
     (E->parallel.me==E->parallel.nprocz-1)){
    /* top heat flux and other stat quantities */
    if (strcmp(E->output.format, "ascii-gz") == 0)
      sprintf(output_file,"%s/qt.dat", E->control.data_dir);
    else
      sprintf(output_file,"%s.qt.dat", E->control.data_file);
    E->output.fpqt = output_open(output_file, "w");
  }else{
    E->output.fpqt = NULL;
  }
  if (E->parallel.me_loc[3] == 0)    {
    /* bottom heat flux and other stat quantities */
    if (strcmp(E->output.format, "ascii-gz") == 0)
      sprintf(output_file,"%s/qb.dat", E->control.data_dir);
    else
      sprintf(output_file,"%s.qb.dat", E->control.data_file);
    E->output.fpqb = output_open(output_file, "w");
  }else{
    E->output.fpqb = NULL;
  }


  return;
}


static void output_parse_optional(struct  All_variables *E)
{
    char* strip(char*);

    int pos, len;
    char *prev, *next;

    len = strlen(E->output.optional);
    /* fprintf(stderr, "### length of optional is %d\n", len); */
    pos = 0;
    next = E->output.optional;

    E->output.connectivity = 0;
    E->output.stress = 0;
    E->output.pressure = 0;
    E->output.surf = 0;
    E->output.botm = 0;
    E->output.geoid = 0;
    E->output.horiz_avg = 0;
    E->output.tracer = 0;
    E->output.comp_el = 0;
    E->output.comp_nd = 0;
    E->output.k = 0;
    E->output.kd = 0;
    E->output.kT = 0;
    E->output.kC = 0;
    E->output.k_total = 0;
    E->output.kappa_eff = 0;
    E->output.rho_ref = 0;
    E->output.Cp = 0;
    E->output.heating = 0;

    while(1) {
        /* get next field */
        prev = strsep(&next, ",");

        /* break if no more field */
        if(prev == NULL) break;

        /* skip if empty */
        if(prev[0] == '\0') continue;

        /* strip off leading and trailing whitespaces */
        prev = strip(prev);

        /* skip empty field */
        if (strlen(prev) == 0) continue;

        /* fprintf(stderr, "### %s: %s\n", prev, next); */
        if(strcmp(prev, "connectivity")==0)
            E->output.connectivity = 1;
        else if(strcmp(prev, "stress")==0)
            E->output.stress = 1;
        else if(strcmp(prev, "pressure")==0)
            E->output.pressure = 1;
        else if(strcmp(prev, "surf")==0)
            E->output.surf = 1;
        else if(strcmp(prev, "botm")==0)
            E->output.botm = 1;
        else if(strcmp(prev, "geoid")==0)
	    if (E->parallel.nprocxy != 12) {
		fprintf(stderr, "Warning: geoid calculation only works in full version. Disabled\n");
	    }
	    else {
		/* geoid calculation requires surface and CMB topo */
		/* make sure the topos are available!              */
		E->output.geoid  = 1;
	    }
        else if(strcmp(prev, "horiz_avg")==0)
            E->output.horiz_avg = 1;
        else if(strcmp(prev, "tracer")==0)
            E->output.tracer = 1;
        else if(strcmp(prev, "comp_el")==0)
            E->output.comp_el = 1;
        else if(strcmp(prev, "comp_nd")==0)
            E->output.comp_nd = 1;
        else if(strcmp(prev, "k")==0)
            E->output.k = 1;
        else if(strcmp(prev, "kd")==0)
            E->output.kd = 1;
        else if(strcmp(prev, "kT")==0)
            E->output.kT = 1;
        else if(strcmp(prev, "kC")==0)
            E->output.kC = 1;
        else if(strcmp(prev, "k_total")==0)
            E->output.k_total = 1;
        else if(strcmp(prev, "kappa_eff")==0)
            E->output.kappa_eff = 1;
        else if(strcmp(prev, "rho_ref")==0)
            E->output.rho_ref = 1;
        else if(strcmp(prev, "Cp")==0 || strcmp(prev, "cp")==0)
            E->output.Cp = 1;
        else if(strcmp(prev, "heating")==0)
            E->output.heating = 1;
        else
            if(E->parallel.me == 0)
                fprintf(stderr, "Warning: unknown field for output_optional: %s\n", prev);

    }

    return;
}

/* check whether E->control.data_file contains a path seperator */
static void chk_prefix(struct  All_variables *E)
{
  char *found;

  found = strchr(E->control.data_prefix, '/');
  if (found) {
      fprintf(stderr, "error in input parameter: datafile='%s' contains '/'\n", E->control.data_file);
      parallel_process_termination();
  }

  if(E->control.log_template[0] == '\0' ||
     strchr(E->control.log_template,'/') != NULL) {
      fprintf(stderr,
              "error in input parameter: logfile='%s' must be a nonempty basename template\n",
              E->control.log_template);
      parallel_process_termination();
  }

  if (E->control.restart || E->control.post_p ||
      (E->convection.tic_method == -1) ||
      (E->control.tracer && (E->trace.ic_method == 2))) {
      found = strchr(E->control.data_prefix_old, '/');
      if (found) {
	  fprintf(stderr, "error in input parameter: datafile_old='%s' contains '/'\n", E->control.data_file);
	  parallel_process_termination();
      }
  }
}


/* search src and substitue the 1st occurance of target by value */
static void expand_str(char *src, size_t max_size,
		       const char *target, const char *value)
{
    char *pos, *end, *new_end;
    size_t end_len, value_len;

    /* is target a substring of src? */
    pos = strstr(src, target);
    if (pos != NULL) {
        value_len = strlen(value);

	/* the end part of the original string... */
	end = pos + strlen(target);
        /* ...and where it is going */
        new_end = pos + value_len;
        end_len = strlen(end);
        if (new_end + end_len >= src + max_size) {
            /* too long */
            return;
        }

	/* move the end part of the original string */
        memmove(new_end, end, end_len + 1); /* incl. null byte */

        /* insert the value */
        memcpy(pos, value, value_len);
    }
}

static void expand_datadir(struct All_variables *E, char *datadir)
{
    char *found, *err;
    char tmp[150];
    int diff;
    FILE *pipe;
    const char str1[] = "%HOSTNAME";
    const char str2[] = "%RANK";
    const char str3[] = "%DATADIR";
    const char str3_prog[] = "citcoms_datadir";

    /* expand str1 by machine's hostname */
    found = strstr(datadir, str1);
    if (found) {
	gethostname(tmp, 100);
	expand_str(datadir, 150, str1, tmp);
    }

    /* expand str2 by MPI rank */
    found = strstr(datadir, str2);
    if (found) {
	sprintf(tmp, "%d", E->parallel.me);
	expand_str(datadir, 150, str2, tmp);
    }

    /* expand str3 by the result of the external program */
    diff = strcmp(datadir, str3);
    if (!diff) {
	pipe = popen(str3_prog, "r");
	err = fgets(tmp, 150, pipe);
	pclose(stdout);
	if (err != NULL)
	    sscanf(tmp, " %s", datadir);
	else {
	    fprintf(stderr, "Cannot get datadir from command '%s'\n", str3_prog);
	    parallel_process_termination();
	}
    }
}


void mkdatadir(const char *dir)
{
  int err;

  err = mkdir(dir, 0755);
  if (err && errno != EEXIST) {
      /* if error occured and the directory is not exisitng */
      fprintf(stderr, "Cannot make new directory '%s'\n", dir);
      parallel_process_termination();
  }
}


void output_init(struct  All_variables *E)
{
    chk_prefix(E);
    expand_datadir(E, E->control.data_dir);
    mkdatadir(E->control.data_dir);
    snprintf(E->control.data_file, 200, "%s/%s", E->control.data_dir,
	     E->control.data_prefix);

    if (E->control.restart || E->control.post_p ||
        (E->convection.tic_method == -1) ||
        (E->control.tracer && (E->trace.ic_method == 2))) {
	expand_datadir(E, E->control.data_dir_old);
	snprintf(E->control.old_P_file, 200, "%s/%s", E->control.data_dir_old,
		 E->control.data_prefix_old);
    }

    open_log(E);
    open_time(E);
    open_info(E);

    if (strcmp(E->output.format, "ascii") == 0) {
        E->problem_output = output;
    }
    else if (strcmp(E->output.format, "hdf5") == 0)
        E->problem_output = h5output;
#ifdef USE_GZDIR
    else if (strcmp(E->output.format, "ascii-gz") == 0)
        E->problem_output = gzdir_output;
    else {
        /* indicate error here */
        if (E->parallel.me == 0) {
            fprintf(stderr, "wrong output_format, must be 'ascii', 'hdf5', or 'ascii-gz'\n");
            fprintf(E->fp, "wrong output_format, must be  'ascii', 'hdf5', or 'ascii-gz'\n");
        }
        parallel_process_termination(E);
    }
#else
    else {
        /* indicate error here */
        if (E->parallel.me == 0) {
            fprintf(stderr, "wrong output_format, must be 'ascii' or 'hdf5' (USE_GZDIR undefined)\n");
            fprintf(E->fp, "wrong output_format, must be 'ascii' or 'hdf5' (USE_GZDIR undefined)\n");
        }
        parallel_process_termination(E);
    }
#endif

    output_parse_optional(E);
    profile_parse_optional(E);
}



void output_finalize(struct  All_variables *E)
{
  char message[255],files[255];
  if (E->fp)
    fclose(E->fp);

  if (E->fptime)
    fclose(E->fptime);

  if (E->fp_out)
    fclose(E->fp_out);

  if (E->trace.fpt)
    fclose(E->trace.fpt);

  if(E->output.fpqt)
    fclose(E->output.fpqt);
  if(E->output.fpqb)
    fclose(E->output.fpqb);

  mat_prop_free(E);


#ifdef USE_GZDIR
  /*
     remove VTK geo file in case we used that for IO
  */
  if((E->output.gzdir.vtk_io != 0) &&
     (strcmp(E->output.format, "ascii-gz") == 0)){
    if((E->output.gzdir.vtk_io == 3)||(E->parallel.me == 0)){
      /* delete the geo files */
      get_vtk_filename(files,1,E,0);
      remove(files);
      if(E->parallel.me == 0){
	/* close the log */
	fclose(E->output.gzdir.vtk_fp);
      }
    }
  }
#endif
}


char* strip(char *input)
{
    int end;
    char *str;
    end = strlen(input) - 1;
    str = input;

    /* trim trailing whitespace */
    while (isspace(str[end]))
        end--;

    str[++end] = 0;

    /* trim leading whitespace */
    while(isspace(*str))
        str++;

    return str;
}
