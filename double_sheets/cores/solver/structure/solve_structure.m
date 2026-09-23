%% Extract variables
X_vec = h_X_vec(:,i_time);

q_vec = X_vec(1:N_q_all,1);
dt_q_vec = X_vec(2*N_q_all+1:3*N_q_all,1);

q_vec_1 = X_vec(N_q_all+1:2*N_q_all,1);
dt_q_vec_1 = X_vec(3*N_q_all+1:4*N_q_all,1);



%% Time integration (Euler predictor-corrector)

%%[*] Fluid force vector: [p] = p_lift + Mf1*dt^2_q + (Mf2_1*(dt_r - Vin - Vwake)^T*dt_ni + Mf2_2)

%%[0] Linear interpolation of the fluid force (Mf2_2)
Qf_p_global_t = @( time) (Qf_p_global - old_Qf_p_global)*(time - time_fluid)/d_t_wake + Qf_p_global_a;
Qf_p_global_t_1 = @( time) (Qf_p_global_1 - old_Qf_p_global_1)*(time - time_fluid)/d_t_wake + Qf_p_global_a_1;

Qf_p_global_tv_n =  Qf_p_global_t( time);
Qf_p_global_t_1v_n =  Qf_p_global_t_1( time);

h_Qf_p_global_t(:,i_time) = Qf_p_global_tv_n;
h_Qf_p_global_t_1(:,i_time) = Qf_p_global_t_1v_n;


%%[1] ( -(τx*dxΓ + τy*dyΓ)*dt_rc ) 
Qf_p_lift2_global_t = @( time) (Qf_p_lift2_mat_global - old_Qf_p_lift2_mat_global)*(time - time_fluid)/d_t_wake + Qf_p_lift2_mat_global_a;
Qf_p_lift2_global_t_1 = @( time) (Qf_p_lift2_mat_global_1 - old_Qf_p_lift2_mat_global_1)*(time - time_fluid)/d_t_wake + Qf_p_lift2_mat_global_a_1;

Qf_p_lift2_global_tv_n = Qf_p_lift2_global_t( time);
Qf_p_lift2_global_t_1v_n = Qf_p_lift2_global_t_1( time);


%%[2] Linear interpolation of the added-mass matrix (Mf1)
Qf_p_mat_global_t = @( time) (Qf_p_mat_global - old_Qf_p_mat_global)*(time - time_fluid)/d_t_wake + Qf_p_mat_global_a;
Qf_p_mat_global_t_1 = @( time) (Qf_p_mat_global_1 - old_Qf_p_mat_global_1)*(time - time_fluid)/d_t_wake + Qf_p_mat_global_a_1;

Qf_p_mat_global_tv_n = Qf_p_mat_global_t( time);
Qf_p_mat_global_t_1v_n = Qf_p_mat_global_t_1( time);


%%[3] Linear interpolation of the fluid force (Mf2_1)
Qf_p_mat0_global_t = @( time) (Qf_p_mat0_global - old_Qf_p_mat0_global)*(time - time_fluid)/d_t_wake + Qf_p_mat0_global_a;
Qf_p_mat0_global_t_1 = @( time) (Qf_p_mat0_global_1 - old_Qf_p_mat0_global_1)*(time - time_fluid)/d_t_wake + Qf_p_mat0_global_a_1;

Qf_p_mat0_global_tv_n = Qf_p_mat0_global_t( time);
Qf_p_mat0_global_t_1v_n = Qf_p_mat0_global_t_1( time);



%% Fluid force evaluation
%%[*] Unit normal vector evaluation
generate_dt_n_vec;


%%[*] Displacement velocity at the collocation points
generate_r_panel;




if ~exist( 'Gamma_wake', 'var')
    
   V_wake_plate = 0; 
   dt_Amat2_Gamma = 0;
   dt_Amat1 = 0;
   Gamma_all = 0;
end


Qf_p_mat0_global_t_n = Qf_p_mat0_global_tv_n*( sum( (dt_rc_vec_all - V_in_all - V_wake_plate - dt_Amat2_Gamma).*dt_n_vec_i_all, 2) - dt_Amat1*Gamma_all);
Qf_p_mat0_global_t_n_1 = Qf_p_mat0_global_t_1v_n*( sum( (dt_rc_vec_all - V_in_all - V_wake_plate - dt_Amat2_Gamma).*dt_n_vec_i_all, 2) - dt_Amat1*Gamma_all);
Qf_p_lift2_global_t_n = Qf_p_lift2_global_tv_n*reshape( dt_rc_vec.', [], 1); 
Qf_p_lift2_global_t_n_1 = Qf_p_lift2_global_t_1v_n*reshape( dt_rc_vec_1.', [], 1); 

%% Predictor (solve for the state at the new time)


%%[0] Assemble the viscoelastic vector Qe^(n)
flag_output = 1;                                                            %% Enable stiffness matrix evaluation
generate_stiff_matrices;                                                    %% Stiffness matrix of sheet 1
generate_stiff_matrices_1;                                                  %% Stiffness matrix of sheet 2

Qk_global_n = Qk_global;                                                    %% Stiffness matrix of sheet 1
Qe_global_n = Qe_global;
dq_Qe_global_n = dq_Qe_global;
Qd_global_n = Qd_global;

Qk_global_n_1 = Qk_global_1;                                                %% Stiffness matrix of sheet 2
Qe_global_n_1 = Qe_global_1;
dq_Qe_global_n_1 = dq_Qe_global_1;
Qd_global_n_1 = Qd_global_1;


%%[1] Predictor
%%[1-0] Mass matrix      : M
m_global_struct.M_global = M_global - Qf_p_mat_global_tv_n(:,1:end/2);
m_global_struct.M_global_1 = M_global_1 - Qf_p_mat_global_t_1v_n(:,end/2+1:end);
m_global_struct.M_global1 =  -Qf_p_mat_global_tv_n(:,end/2+1:end);
m_global_struct.M_global1_1 = -Qf_p_mat_global_t_1v_n(:,1:end/2);
%%[1-1] Body force       : Q_f
qf_global_struct.Qf_global = Qf_global + Qf_time_global*q_in_norm( time) + Qf_p_global_tv_n + Qf_p_mat0_global_t_n + Qf_p_lift2_global_t_n;
qf_global_struct.Qf_global_1 = Qf_global + Qf_time_global*q_in_norm_1( time) + Qf_p_global_t_1v_n + Qf_p_mat0_global_t_n_1 + Qf_p_lift2_global_t_n_1;
%%[1-2] dq_Qe(q(n))
dq_qe_global_struct.dq_Qe_global = dq_Qe_global_n;
dq_qe_global_struct.dq_Qe_global_1 = dq_Qe_global_n_1;
%%[1-3] Elastic force    : Q_e
qe_global_struct.Qe_global = Qe_global_n + Qk_global_n;
qe_global_struct.Qe_global_1 = Qe_global_n_1 + Qk_global_n_1;
%%[1-4] Damping force    : Q_d
qd_global_struct.Qd_global = Qd_global_n;
qd_global_struct.Qd_global_1 = Qd_global_n_1;



[ X_vec_p, out1] = new_X_func_FAST( X_vec, m_global_struct, qf_global_struct, dq_qe_global_struct, qe_global_struct, qd_global_struct, var_param, 0, []);

                                


%% Corrector (solve with the elastic force Qe^(n+1) at the new time)

q_vec = X_vec_p(1:N_q_all,1);
dt_q_vec = X_vec_p(2*N_q_all+1:3*N_q_all,1);

q_vec_1 = X_vec_p(N_q_all+1:2*N_q_all,1);
dt_q_vec_1 = X_vec_p(3*N_q_all+1:4*N_q_all,1);


%%[*] Unit normal vector evaluation
generate_dt_n_vec;


%%[*] Displacement velocity at the collocation points

dt_rc_vec = Sc_mat_col_global*dt_q_vec;
dt_rc_vec = reshape( dt_rc_vec, 3, []).';

dt_rc_vec_1 = Sc_mat_col_global*dt_q_vec_1;
dt_rc_vec_1 = reshape( dt_rc_vec_1, 3, []).';

dt_rc_vec_all = [   dt_rc_vec;
                    dt_rc_vec_1];


%%[*] Predict the fluid force at t+dt
%%[*] Fluid force vector: [p] = p_lift + Mf1*dt^2_q + (Mf2_1*(dt_r - Vin - Vwake)^T*dt_ni + Mf2_2)
Qf_p_global_tv_np1 =  Qf_p_global_t( time + d_t);
Qf_p_global_t_1v_np1 =  Qf_p_global_t_1( time + d_t);

%%[1] ( -(τx*dxΓ + τy*dyΓ)*dt_rc ) 
Qf_p_lift2_global_tv_np1 = Qf_p_lift2_global_t( time + d_t);
Qf_p_lift2_global_t_1v_np1 = Qf_p_lift2_global_t_1( time + d_t);


%%[2] Linear interpolation of the added-mass matrix (Mf1)
Qf_p_mat_global_tv_np1 = Qf_p_mat_global_t( time + d_t);
Qf_p_mat_global_t_1v_np1 = Qf_p_mat_global_t_1( time + d_t);


%%[3] Linear interpolation of the fluid force (Mf2_1)
Qf_p_mat0_global_tv_np1 = Qf_p_mat0_global_t( time + d_t);
Qf_p_mat0_global_t_1v_np1 = Qf_p_mat0_global_t_1( time + d_t);



Qf_p_mat0_global_t_np1 = Qf_p_mat0_global_tv_np1*( sum( (dt_rc_vec_all - V_in_all - V_wake_plate - dt_Amat2_Gamma).*dt_n_vec_i_all, 2) - dt_Amat1*Gamma_all);
Qf_p_mat0_global_t_np1_1 = Qf_p_mat0_global_t_1v_np1*( sum( (dt_rc_vec_all - V_in_all - V_wake_plate - dt_Amat2_Gamma).*dt_n_vec_i_all, 2) - dt_Amat1*Gamma_all);
Qf_p_lift2_global_t_np1 = Qf_p_lift2_global_tv_np1*reshape( dt_rc_vec.', [], 1); 
Qf_p_lift2_global_t_np1_1 = Qf_p_lift2_global_t_1v_np1*reshape( dt_rc_vec_1.', [], 1); 



%%[0] Assemble the viscoelastic vector Qe^(n+1)
flag_output = 0;                                                            %% Disable stiffness matrix evaluation
generate_stiff_matrices;                                                    %% Stiffness matrix of sheet 1
generate_stiff_matrices_1;                                                  %% Stiffness matrix of sheet 2

Qk_global_np1 = Qk_global;                                                  %% Stiffness matrix of sheet 1
Qk_global_np1_1 = Qk_global_1;                                              %% Stiffness matrix of sheet 2




%%[1-0] Mass matrix      : M
m_global_struct.M_global = M_global - ( Qf_p_mat_global_tv_n(:,1:end/2) + Qf_p_mat_global_tv_np1(:,1:end/2) )/2;
m_global_struct.M_global_1 = M_global_1 - ( Qf_p_mat_global_t_1v_n(:,end/2+1:end) + Qf_p_mat_global_t_1v_np1(:,end/2+1:end) )/2;
m_global_struct.M_global1 =  -( Qf_p_mat_global_tv_n(:,end/2+1:end) + Qf_p_mat_global_tv_np1(:,end/2+1:end) )/2;
m_global_struct.M_global1_1 = -( Qf_p_mat_global_t_1v_n(:,1:end/2) + Qf_p_mat_global_t_1v_np1(:,1:end/2) )/2;
%%[1-1] Body force       : Q_f
qf_global_struct.Qf_global = Qf_global + Qf_time_global*q_in_norm( time) ...
                                        + ( Qf_p_global_tv_n + Qf_p_global_tv_np1 )/2 ...
                                        + ( Qf_p_mat0_global_t_n + Qf_p_mat0_global_t_np1 )/2 ...
                                        + ( Qf_p_lift2_global_t_n + Qf_p_lift2_global_t_np1 )/2;
qf_global_struct.Qf_global_1 = Qf_global + Qf_time_global*q_in_norm_1( time) ...
                                            + ( Qf_p_global_t_1v_n + Qf_p_global_t_1v_np1 )/2 ...
                                            + ( Qf_p_mat0_global_t_n_1 + Qf_p_mat0_global_t_np1_1 )/2 ...
                                            + ( Qf_p_lift2_global_t_n_1 + Qf_p_lift2_global_t_np1_1 )/2;
%%[1-2] dq_Qe(q(n))
dq_qe_global_struct.dq_Qe_global = dq_Qe_global_n;
dq_qe_global_struct.dq_Qe_global_1 = dq_Qe_global_n_1;
%%[1-3] Elastic force    : Q_e
%%[*] Qe(q(n+1)) = Qe(q(n)) + Δt*dq_Qe(q(n))*dt_q(n+1)
qe_global_struct.Qe_global = Qe_global_n + (Qk_global_n + Qk_global_np1)/2;
qe_global_struct.Qe_global_1 = Qe_global_n_1 + (Qk_global_n_1 + Qk_global_np1_1)/2;
%%[1-4] Damping force    : Q_d
qd_global_struct.Qd_global = Qd_global_n;
qd_global_struct.Qd_global_1 = Qd_global_n_1;


new_X_vec = new_X_func_FAST( X_vec, m_global_struct, qf_global_struct, dq_qe_global_struct, qe_global_struct, qd_global_struct, var_param, 1, out1);

h_X_vec(:,i_time+1) = new_X_vec;                                            %% X(n+1) solved with (Qe^(n)+Qe^(n+1))/2


%% Stop the analysis once it diverges
if sum( isnan( X_vec))

    warndlg( 'Divergence!!') 
end