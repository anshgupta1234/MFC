import json
import math

Mu = 1.84E-05
gam_a = 1.4
Ma = 0.5
Re = 40

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    # For these computations, the cylinder is placed at the (0,0,0)
    # domain origin. 
    # axial direction
    'x_domain%beg'                 : -6,
    'x_domain%end'                 : 12,
    # r direction
    'y_domain%beg'                 : -5,
    'y_domain%end'                 : 5,
    'cyl_coord'                    : 'F',
    'm'                            : 2100,
    'n'                            : 700,
    'p'                            : 0,
    'dt'                           : 3.5E-4,
    't_step_start'                 : 0,
    't_step_stop'                  : 44000,  #3000
    't_step_save'                  : 4000,  #10
    # ==========================================================================
    
    # Simulation Algorithm Parameters ==========================================
    # Only one patches are necessary, the air tube
    'num_patches'                  : 1,
    # Use the 5 equation model
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    # One fluids: air
    'num_fluids'                   : 1,
    # time step
    'mpp_lim'                      : 'F',
    # Correct errors when computing speed of sound
    'mixture_err'                  : 'T',
    # Use TVD RK3 for time marching
    'time_stepper'                 : 3,
    # Use WENO5
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'T',
    'weno_avg'                     : 'T',
    'avg_state'                    : 2,
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    # Grid stretching
    'stretch_x'                    : 'T',
    'stretch_y'                    : 'T',
    'x_a'                          : -0.5,
    'x_b'                          : 0.5,
    'a_x'                          : 5,
    'y_a'                          : -0.5,
    'y_b'                          : 0.5,
    'a_y'                          : 3,
    # We use ghost-cell 
    'bc_x%beg'                     : -17,
    'bc_x%end'                     : -3,
    'bc_y%beg'                     : -3,
    'bc_y%end'                     : -3,
    # Set IB to True and add 1 patch
    'ib'                           : 'T',
    'num_ibs'                      : 1,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'E_wrt'                        :'T',
    'c_wrt'                        :'T',
    'parallel_io'                  :'T',
    # ==========================================================================

    # Patch: Constant Tube filled with air =====================================
    # Specify the cylindrical air tube grid geometry
    'patch_icpp(1)%geometry'       : 3,
    'patch_icpp(1)%x_centroid'     : 5,
    # Uniform medium density, centroid is at the center of the domain
    'patch_icpp(1)%y_centroid'     : 0,
    'patch_icpp(1)%length_x'       : 60,
    'patch_icpp(1)%length_y'       : 30,  
    # Specify the patch primitive variables 
    'patch_icpp(1)%vel(1)'         : Ma,
    'patch_icpp(1)%vel(2)'         : 0.0E+00,
    'patch_icpp(1)%pres'           : 1/gam_a,
    'patch_icpp(1)%alpha_rho(1)'   : 1.E+00,
    'patch_icpp(1)%alpha(1)'       : 1.E+00,
    # # ========================================================================

    # Patch: Cylinder Immersed Boundary ========================================
    'patch_ib(1)%geometry'       : 2,
    'patch_ib(1)%x_centroid'     : 0,
    'patch_ib(1)%y_centroid'     : 0,
    'patch_ib(1)%radius'         : 0.5,
    'patch_ib(1)%slip'           : 'F',
    # # ========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(gam_a-1.E+00),  # 2.50(Not 1.40)
    'fluid_pp(1)%pi_inf'           : 0,
    'fluid_pp(1)%Re(1)'            : Re/Ma,
    # ==========================================================================
}))
