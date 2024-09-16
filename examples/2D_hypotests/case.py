import math
import json

# Numerical setup
Nx      = 301  # Number of grid points in x
Ny      = 301  # Number of grid points in y
dx      = 1./(1.*(Nx+1))  # Grid spacing in x
dy      = 1./(1.*(Ny+1))  # Grid spacing in y

Tend    = 64E-06  # End time
Nt      = 3000  # Number of time steps
mydt    = Tend/(1.*Nt)  # Time step size

# Configuring case dictionary
print(json.dumps({
                    # Logistics ================================================                 
                    'run_time_info'                : 'F',                       
                    # ==========================================================
                                                                                
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : 0.E+00,  # x start                    
                    'x_domain%end'                 : 1.E+00,  # x end                    
                    'y_domain%beg'                 : 0.E+00,  # y start
                    'y_domain%end'                 : 1.E+00,  # y end 
                    'm'                            : Nx,  # Number of grid points in x direction                       
                    'n'                            : Ny,  # Number of grid points in y direction                       
                    'p'                            : 0,  # Number of grid points in z (for 3D, change this)                      
                    'dt'                           : 5e-7,  # Time step size                      
                    't_step_start'                 : 0,  # Start time                         
                    't_step_stop'                  : Nt,  # End time                         
                    't_step_save'                  : 100,  # Save frequency
		    # ==========================================================
                                                                                
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 2,  # Two patches                        
                    'model_eqns'                   : 2,  # Number of model equations                       
                    'alt_soundspeed'               : 'F',                      
                    'num_fluids'                   : 1,                        
		        'mpp_lim'                      : 'F',                      
		    # ' mixture_err'                  : 'F',                      
		    'time_stepper'                 : 3,                                               
                    'weno_order'                   : 5,                        
                    'weno_eps'                     : 1.E-16,
		    'weno_Re_flux'                 : 'F',  
    		    'weno_avg'                     : 'F',
                    'mapped_weno'                  : 'F',                     
                    'null_weights'                 : 'F',                      
                    'mp_weno'                      : 'F',                      
		    'riemann_solver'               : 1,                        
                    'wave_speeds'                  : 1,                        
                    'avg_state'                    : 2,                        
                    'bc_x%beg'                     : -3,                       
                    'bc_x%end'                     : -3,                       
                    'bc_y%beg'                     : -3,  # Boundary conditions for y direction
                    'bc_y%end'                     : -3,                       
                    # ==========================================================

                    # Turning on Hypoelasticity ================================
                    'hypoelasticity'               : 'T',                      
                    # ==========================================================
                                                                               
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        
                    'precision'                    : 2,                        
                    'prim_vars_wrt'                :'T',                       
		    'parallel_io'                  :'T',                       
		    # ==========================================================
                                                                                
		    # Patch 1 (background flow) ===================
                    'patch_icpp(1)%geometry'       : 3,  # 2D geometry                     
                    'patch_icpp(1)%x_centroid'     : 0.5,  # x-center                    
                    'patch_icpp(1)%y_centroid'     : 0.5,  # y-center
                    'patch_icpp(1)%length_x'       : 1.0,   # x-length                    
                    'patch_icpp(1)%length_y'       : 1.0,   # y-length
                    'patch_icpp(1)%vel(1)'         : '300*(y-0.5)',   # x-velocity  
                    'patch_icpp(1)%vel(2)'         : 0.0,   # y-velocity   
                    'patch_icpp(1)%pres'           : 1.E+5,  # Pressure                    
                    'patch_icpp(1)%alpha_rho(1)'   : 1050,  # Density                    
                    'patch_icpp(1)%alpha(1)'       : 1.,                
                    'patch_icpp(1)%tau_e(1)'       : 0.0,                
                    # ==========================================================

                    # Patch 2 (hypo material in the center) ================
                    'patch_icpp(2)%geometry'       : 22,  # 2D geometry                     
                    'patch_icpp(2)%x_centroid'     : 0.3,  # x-center                    
                    'patch_icpp(2)%y_centroid'     : 0.4,  # y-center                    
                    'patch_icpp(2)%c'              : 0.4,   # x-length                   
                    'patch_icpp(2)%t'              : 0.10,
                    'patch_icpp(2)%p'              : 0.4,
                    'patch_icpp(2)%m'              : 0.02,              
                    'patch_icpp(2)%theta'          : -45,              
                    'patch_icpp(2)%vel(1)'         : '300*(y-0.5)',   # x-velocity  
                    'patch_icpp(2)%vel(2)'         : 0.0,   # y-velocity   
                    'patch_icpp(2)%pres'           : 1.E+5,  # Pressure                 
                    'patch_icpp(2)%alpha_rho(1)'   : 1000,  # Density                 
                    'patch_icpp(2)%alpha(1)'       : 1.,                
                    'patch_icpp(2)%tau_e(1)'       : 0.0, 
                    'patch_icpp(2)%alter_patch(1)' : 'T',

                    # Fluids Physical Parameters ===============================
                    'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),   
                    'fluid_pp(1)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00 - 1.E+00), 
                    'fluid_pp(1)%G'                : 10E+09,                       
	            # ==========================================================
}))
# ==============================================================================