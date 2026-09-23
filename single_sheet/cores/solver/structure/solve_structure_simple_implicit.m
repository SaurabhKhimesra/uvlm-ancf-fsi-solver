%% Extract variables
X_vec = h_X_vec(:,i_time);
q_vec = X_vec(1:N_q_all,1);

%% Stiffness matrix assembly
generate_stiff_matrices;

%% Time integration (Euler predictor-corrector)

%%[0] Linear interpolation of the fluid force
Qf_p_global_t = (Qf_p_global - old_Qf_p_global)*(time - time_fluid)/d_t_wake + old_Qf_p_global;
h_Qf_p_global_t(:,i_time) = Qf_p_global_t;

%%[1] Linear interpolation of the added-mass matrix
Qf_p_mat_global_t = (Qf_p_mat_global - old_Qf_p_mat_global)*(time - time_fluid)/d_t_wake + old_Qf_p_mat_global;


new_X_vec = new_X_func( X_vec, M_global-Qf_p_mat_global_t, Qf_global+Qf_p_global_t, Qe_global, Qd_global, var_param);

h_X_vec(:,i_time+1) = new_X_vec;

%% Stop the analysis once it diverges
if sum( isnan( X_vec))

    warndlg( 'Divergence!!') 
    break;
end