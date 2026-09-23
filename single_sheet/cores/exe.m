clc
clear all
close all hidden

%% delete
delete( '*.asv')
delete( '*.log')

%% path
add_pathes

%% parameter
param_setting

%% version check
version_check


%% multi threading
maxNumCompThreads( core_num);


%% exe

%%[0] Definition of shape functions
generate_shape_function;

%%[1] Mesh generation
%%[1-0] FEM mesh generation
generate_elements;
%%[1-1] Vortex lattice method grid generation
generate_panel;

%%[1-2] Matrix generation
generate_matrices;

%%[1-3] Disturbance input
generate_Qf_time_mat;




%%[2] Time integration [-]

time_m = 0:d_t:End_Time;
initial_values;


%%[3] solve modes
solve_mode;


%%[4] Speed check
if speed_check == 1
    time_m = time_m(1:10);
    profile on;
end

tic
i_time = 1;
i_time_cnt = 1;
i_wake_time = 1;
time_fluid = 0;
d_t_wake = d_t*dt_wake_per_dt;
time_wake_m = 0;
fluid_compute_flag = 1;
time = 0;
measure_time_struct = 0;
measure_time_fluid = 0;
while time <= time_m(end) || ~fluid_compute_flag
    
    time = i_time*d_t;
    
    disp( [ 'Time = ', num2str( time, '%0.4f'), ' [-]'])

    %%[5] Structural solve
    measure_time_struct_tmp = toc;
    if flag_fluid_bench
        %%[*] Structure: rigid plate
        rigid_structure;   
    else
        %%[*] Structure: elastic plate
        solve_structure;   
    end
    measure_time_struct = measure_time_struct + (toc - measure_time_struct_tmp);
    
    
    %%[6] Fluid solve
    if mod( i_time, dt_wake_per_dt) == 1       
        drawnow                                         	%% Keep the UI from freezing during the run
        
        if fluid_compute_flag
                        
            %% Update the value held from the previous step
            old_Qf_p_global = Qf_p_global;
            old_Qf_p_mat_global = Qf_p_mat_global;
            old_Qf_p_mat0_global = Qf_p_mat0_global;
            old_Qf_p_lift2_mat_global = Qf_p_lift2_mat_global;
        
            measure_time_fluid_tmp = toc;
            solve_fluid;    
            measure_time_fluid = measure_time_fluid + (toc - measure_time_fluid_tmp);
            
            i_wake_time = i_wake_time + 1;   
        else
            
            %% Update the value held from the previous step
            Qf_p_global_a = Qf_p_global;
            Qf_p_mat_global_a = Qf_p_mat_global;
            Qf_p_mat0_global_a = Qf_p_mat0_global;
            Qf_p_lift2_mat_global_a = Qf_p_lift2_mat_global;
                           
            time_fluid = time;                                  %% Fluid parameter update time [-] (for interpolating the fluid force in time)
        end
                
        %%%
        %%% Analytical parameters
        %%%
        disp( [ 'Ma = ', num2str( Ma), ' [-], Ua = ', num2str( Ua), ' [-], theta_a = ', num2str( theta_a), ' [-], C_theta_a = ', num2str( C_theta_a), ' [-], J_a = ', num2str( J_a), ' [-]'])
        %%%
        %%% FSI coupling scheme
        %%%
        if coupling_flag 
            coupling_str = 'Strong coupling';
        else
            coupling_str = 'Weak coupling';
        end
        disp( coupling_str)
    end
    
    
    %%[7] Energy balance evaluation
    solve_energy;
    
    
    %%[8] Intermediate save of the analysis data
    if mod( i_time, 500) == 0 && ~fluid_compute_flag
        save ./save/NUM_DATA -v7.3
    end
    
    %%[9] Iteration
    if mod( i_time, dt_wake_per_dt) == 1      
        
        if fluid_compute_flag
            
            i_time = i_time - i_time_cnt;            
            fluid_compute_flag = 0;
        else
            
            i_time_cnt = 0;
            fluid_compute_flag = 1;
        end
    end
    
    i_time = i_time + 1;
    i_time_cnt = i_time_cnt + 1;
end
measure_time_all = toc;



%% save 

%%[*] Speed check
if speed_check == 1
    profile viewer;
else
    save ./save/NUM_DATA -v7.3
end