%% [0] Assign the initial structural parameters

h_X_vec = zeros(2*N_q_all,length( time_m));
N_time = length( time_m);
N_fluid_time = ceil( N_time/dt_wake_per_dt);

idx_r = reshape( [1:N_qi:N_q_all; 2:N_qi:N_q_all], 1, []);                      %% Assign the nodal coordinates of the initial shape [-]
h_X_vec(idx_r,1) = reshape( coordinates.',[],1);
idx_dx_r = reshape( [4:N_qi:N_q_all; 5:N_qi:N_q_all; 6:N_qi:N_q_all], 1, []); 	%% Assign the x-direction slope of the initial shape (dx_r = [1 0 0]^T [-])
h_X_vec(idx_dx_r,1) = repmat( [1 0 0].', [ N_node 1]);
idx_dy_r = reshape( [7:N_qi:N_q_all; 8:N_qi:N_q_all; 9:N_qi:N_q_all], 1, []); 	%% Assign the y-direction slope of the initial shape (dy_r = [0 1 0]^T [-])
h_X_vec(idx_dy_r,1) = repmat( [0 1 0].', [ N_node 1]);

idx_r = reshape( 3:N_qi:N_q_all, 1, []);                                        %% Assign the initial nodal coordinate velocity (Z direction) [-]
dt_rz0 = dt_rz_end*( coordinates(:,1)/Length).^2;                               %% dt_rz ∝ X^2 ∈ [0,1]
h_X_vec(idx_r + N_q_all,1) = reshape( dt_rz0.',[],1);                      



%% [1] Assign the initial fluid parameters
h_r_panel_vec = zeros(4*N_element,3,N_fluid_time);                            	%% Store the panel node points [-]
h_rcol_vec = zeros(N_element,3,N_fluid_time);                                   %% Store the collocation points [-]
h_n_vec = zeros(N_element,3,N_fluid_time);                                      %% Unit normal vector [-]
h_dt_n_vec = zeros(N_element,3,N_fluid_time);                                 	%% Unit normal vector [-]
h_Amat = zeros(N_element,N_element,N_fluid_time);                               %%  AΓ=b [-]
h_dt_Amat = zeros(N_element,N_element,N_fluid_time);                            %%  AΓ=b [-]

h_V_surf = zeros(N_element,3,N_fluid_time);                                     %% Surface velocity vector [-]
h_V_wake_end = zeros(Ny+1,3,N_fluid_time);                                      %% Trailing-edge velocity vector [-]


h_r_end = zeros(Ny+1,3,N_fluid_time);                                           %% Trailing-edge nodal coordinate vector [-]
h_dp_add = zeros(N_element,N_fluid_time);                                       %% Added-mass contribution (Γ_n+1 - Γ_n)/dt [-]
h_dp_add_estimate = zeros(N_element,N_fluid_time);                            	%% Added-mass contribution [-]
h_dp_lift = zeros(N_element,N_fluid_time);                                      %% Steady fluid force contribution [-]
h_dp_vec = zeros(N_element,N_fluid_time);                                       %% Fluid force [-]
h_Qf_p_global_t = zeros(N_q_all,N_time);                                        %% Fluid force matrix [-]

old_Qf_p_global = zeros(N_q_all,1);
old_Qf_p_mat_global = zeros(N_q_all);
old_Qf_p_mat0_global = zeros(N_q_all,N_element);
old_Qf_p_lift2_mat_global = zeros(N_q_all,3*N_element);
Qf_p_global_a = old_Qf_p_global;
Qf_p_mat_global_a = old_Qf_p_mat_global;
Qf_p_mat0_global_a = old_Qf_p_mat0_global;
Qf_p_lift2_mat_global_a = old_Qf_p_lift2_mat_global;
Qf_p_global = old_Qf_p_global;
Qf_p_mat_global = old_Qf_p_mat_global;
Qf_p_mat0_global = old_Qf_p_mat0_global;
Qf_p_lift2_mat_global = old_Qf_p_lift2_mat_global;


%% [2] Energy and work-rate time histories
h_E_inertia = zeros(1,length( time_m));
h_E_em = h_E_inertia;
h_E_ek = h_E_inertia;
h_E_Ja = h_E_inertia;

h_W_dm = h_E_inertia;
h_W_dk = h_E_inertia;
h_W_em2 = h_E_inertia;
h_W_ek2 = h_E_inertia;
h_W_d_theta = h_E_inertia;

h_W_f_ext = zeros(1,N_fluid_time);
h_W_f_ext_Xdist = zeros(N_fluid_time,Nx);