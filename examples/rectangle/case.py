#!/usr/bin/env python3

import json
import math

# Dynamic Viscosity
Mu1 = 0.0000184
#Mu2 = 0.01
rho1 = 0.2199
gam_a = 1.4


# Patch Design
D = 0.1

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    # For these computations, the cylinder is placed at the (0,0,0)
    # domain origin. 
    # axial direction
    'x_domain%beg'                 : -5*D,
    'x_domain%end'                 : 5*D,
    'y_domain%beg'                 : -2.5*D,
    'y_domain%end'                 : 2.5*D,
    'm'                            : 80,
    'n'                            : 40,
    'p'                            : 0,
    'dt'                           : 2.0E-06,
    't_step_start'                 : 0,
    't_step_stop'                  : 5000,
    't_step_save'                  : 10,
    # ==========================================================================
    
    # Simulation Algorithm Parameters ==========================================
    # Only one patches are necessary, the air tube
    'num_patches'                  : 1,
    # Use the 5 equation model
    'model_eqns'                   : 2,
    # 6 equations model does not need the K \div(u) term                    
    'alt_soundspeed'               : 'F',
    # One fluids: air
    'num_fluids'                   : 1,
    # Advect both volume fractions
    # No need to ensure the volume fractions sum to unity at the end of each
    # time step
    # Correct errors when computing speed of sound
    # Use TVD RK3 for time marching
    'time_stepper'                 : 3,
    # Use WENO5
    'weno_order'                   : 3,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',
    'weno_avg'                     : 'F',
    'avg_state'                    : 2,
    # Use the mapped WENO weights to maintain monotinicity
    'mapped_weno'                  : 'F',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    # Use the HLLC  Riemann solver
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    # We use reflective boundary conditions at octant edges and 
    # non-reflective boundary conditions at the domain edges
    'bc_x%beg'                     : -3,
    'bc_x%end'                     : -3,
    'bc_y%beg'                     : -3,
    'bc_y%end'                     : -3,
    'ib'		                   : 'T',
    'num_ibs'                      :  1,
    'perturb_flow'                 : 'T',
    'perturb_flow_fluid'           :  1,
    'perturb_flow_mag'             : 0.01,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    # Export primitive variables in double precision with parallel
    # I/O to minimize I/O computational time during large simulations
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    #'omega_wrt(1)'                :'T',  
    #'omega_wrt(2)'                :'T',  
    'omega_wrt(3)'                 :'T',  
    'fd_order'                     : 2, 
    # ==========================================================================

    # Patch: Middle =====================================
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 0, 
    'patch_icpp(1)%y_centroid'     : 0,
    'patch_icpp(1)%length_x'       : 1000*D,
    'patch_icpp(1)%length_y'       : 1000*D,
    #Specify the patch primitive variables 9
    'patch_icpp(1)%vel(1)'         : 527.2E+00,
    'patch_icpp(1)%vel(2)'         : 0.0E+00,
    'patch_icpp(1)%pres'           : 10918.2549,
    'patch_icpp(1)%alpha_rho(1)'   : (1.0)*rho1,
    'patch_icpp(1)%alpha(1)'       : 1.0,
    # =========================================================================
    # 'patch_ib(1)%geometry'          : 3,
    # 'patch_ib(1)%x_centroid'        : 0,
    # 'patch_ib(1)%y_centroid'        : 0,
    # 'patch_ib(1)%length_x'          : D,
    # 'patch_ib(1)%length_y'          : D,
    # 'patch_ib(1)%slip'              : 'F',

    'patch_ib(1)%geometry'                     : 5,
    'patch_ib(1)%model%filepath'               : '/Users/anshgupta1234/Desktop/Coding/MFC-copy/examples/rectangle/Rec_IBM.stl',
    'patch_ib(1)%model%translate(1)'           : -0.05,
    'patch_ib(1)%model%translate(2)'           : -0.05,
    'patch_ib(1)%model%spc'                    : 100,
    'patch_ib(1)%model%threshold'              : 0.05,
    'patch_ib(1)%slip'                         : 'F',
    # # ==========================================================================
   
    # Fluids Physical Parameters ===============================================
    # Use the same stiffness as the air bubble
    'fluid_pp(1)%gamma'            : 1.E+00/(gam_a-1.E+00),  # 2.50(Not 1.40)
    'fluid_pp(1)%pi_inf'           : 0,
    'fluid_pp(1)%Re(1)'            : 7535533.2,
    # ==========================================================================
}))

