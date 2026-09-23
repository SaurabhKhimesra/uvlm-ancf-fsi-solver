

q_vec = X_vec(1:N_q_all);
dt_q_vec = X_vec(N_q_all+1:end);







%% Panel node point evaluation [m]
                         
generate_r_panel;

%%[0] Time history of the panel node points [-]
%%% 1:element index, 2:coordinate component, 3:time
h_r_panel_vec(:,:,i_wake_time) = [ r_panel_vec_1;
                                   r_panel_vec_2;
                                   r_panel_vec_3;
                                   r_panel_vec_4];
                               

                            
%% Collocation point evaluation [m]
%%% 1:element index, 2:coordinate component, 3:time
rc_vec = Sc_mat_col_global*q_vec;
rc_vec = reshape( rc_vec, 3, []).';
h_rcol_vec(:,:,i_wake_time) = rc_vec;





%% Unit normal vector evaluation

%%[1-0] n_vec
generate_dt_n_vec;




%%% 1:element index, 2:coordinate component, 3:time
h_n_vec(:,:,i_wake_time) = n_vec_i;
h_dt_n_vec(:,:,i_wake_time) = dt_n_vec_i;


%% Generation of influence coefficient matrix
%%%
%%%    ^ Y
%%%    | (1)-----2--------(2)
%%%    | |               |
%%%    | 1       X       2
%%%    | |               |
%%%    | (4)-----4--------(3)
%%%    |--------------------------> X
%%%


%%[2-0] Generation of influence coefficient matrix
q1234_mat = generate_q1234_mat( rc_vec, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, var_param, 1);                                    


%%[2-1] Generation of influence coefficient matrix
n_vec_i_mat = repmat( n_vec_i, [ 1 N_element]);
A_mat = inner_mat( q1234_mat, n_vec_i_mat);
A_mat = A_mat(:,1:3:end);                                                   %% Copying only the 3 components inside inner_mat is unnecessary.


%% Initial wake parameters [-]
if ~exist( 'Gamma_wake', 'var')

    old_Gamma = zeros(N_element,1);    
    Gamma = zeros(N_element,1);
    Gamma_trail = zeros(Ny,1);
    Qf_p_global = zeros(N_q_all,1);
    
end


%% Displacement velocity at the collocation points [-]


%%% 1:element index, 2:coordinate component
dt_rc_vec = Sc_mat_col_global*dt_q_vec;
dt_rc_vec = reshape( dt_rc_vec, 3, []).';





%% Wake generation

generate_wake;


%% Wake-induced velocity [-]

[ V_wake_plate, q1234_wake_mat] = V_wake_func( rc_vec, r_wake_1, r_wake_2, r_wake_3, r_wake_4, Gamma_wake, var_param, 1);
    
%%[*] --------------------------------------------------------------------added (18:36 2017/03/30)
%%[*] Velocity induced on the plate by the vortex sheet with the freshly shed circulation removed [-]
n_vec_i_mat2 = repmat( n_vec_i, [ 1 size( r_wake_1, 1)]);
q1234_wake_mat_n = inner_mat( q1234_wake_mat, n_vec_i_mat2);
q1234_wake_mat_n = q1234_wake_mat_n(:,1:3:end);                        	%% Copying only the 3 components inside inner_mat is unnecessary.
Gamma_wake2 = zeros(size( r_wake_1, 1),1);
Gamma_wake2(Ny+1:end) = Gamma_wake(Ny+1:end);
V_wake_plate2_n = sum( q1234_wake_mat_n.*(ones(N_element,1)*Gamma_wake2.'), 2);




%% Normal inlet and body velocity [-]


V_normal = sum( (dt_rc_vec - V_in - V_wake_plate).*n_vec_i, 2);



%%  dt_u_wake = B*dt_Γ + Σ{Γwake*(dt_q_mat)^T*n}
% n_vec_i_wake_mat = repmat( n_vec_i, [ 1 size( q1234_wake_mat, 2)/3]);
% q1234_wake_n = inner_mat( q1234_wake_mat, n_vec_i_wake_mat);
% q1234_wake_n = q1234_wake_n(:,1:3:end);                                     %% Copying only the 3 components inside inner_mat is unnecessary.
% 
% B_mat = [ zeros(N_element,N_element-Ny) q1234_wake_n(:,1:Ny)];


q1234_wake_trail_mat = generate_q1234_mat( rc_vec, r_wake_1(1:Ny,:), r_wake_2(1:Ny,:), r_wake_3(1:Ny,:), r_wake_4(1:Ny,:), var_param, 1);   
n_vec_i_trail_wake_mat = repmat( n_vec_i, [ 1 size( q1234_wake_trail_mat, 2)/3]);
q1234_wake_trail_n = inner_mat( q1234_wake_trail_mat, n_vec_i_trail_wake_mat);
q1234_wake_trail_n = q1234_wake_trail_n(:,1:3:end);                       	%% Copying only the 3 components inside inner_mat is unnecessary.

B_mat = [ zeros(N_element,N_element-Ny) q1234_wake_trail_n];




%% Circulation evaluation

Gamma = A_mat\V_normal;

%%% Circulation of the trailing-edge panel (used for the wake circulation via the unsteady Kutta condition: Dt_Γ=V^T*∇Γ+∂t_Γ=0 -> vorticity convection)
Gamma_trail = old_Gamma(end-Ny+1:end);    

%%--------------------------------------------------------------------added (18:36 2017/03/30)
% V_normal2 = sum( (dt_rc_vec - V_in).*n_vec_i, 2) - V_wake_plate2_n;

% Gamma = (A_mat + B_mat)\V_normal2;
%%% Circulation of the trailing-edge panel (used for the wake circulation via the Kutta condition)
% Gamma_trail = Gamma(end-Ny+1:end);  


h_Gamma(i_wake_time) = { Gamma };


%% Update the wake-induced velocity
Gamma_wake(1:Ny) = Gamma_trail;
Gamma_wake_no_trail = Gamma_wake;
Gamma_wake_no_trail(1:Ny) = 0;

V_wake_plate_trail = V_wake_func( rc_vec, r_wake_1(1:Ny,:), r_wake_2(1:Ny,:), r_wake_3(1:Ny,:), r_wake_4(1:Ny,:), Gamma_trail, var_param, 1);
q_gamma_wake_no_trail = q1234_wake_mat.*( ones(N_element,1)*kron( Gamma_wake_no_trail.', ones(1,3)) );    
V_wake_plate_no_trail = [ sum( q_gamma_wake_no_trail(:,1:3:end), 2) sum( q_gamma_wake_no_trail(:,2:3:end), 2) sum( q_gamma_wake_no_trail(:,3:3:end), 2)];
V_wake_plate = V_wake_plate_trail + V_wake_plate_no_trail;





%% Surface velocity distribution (at the collocation points)

q_gamma = q1234_mat.*( ones(N_element,1)*kron( Gamma.', ones(1,3)) );    
V_gamma = [ sum( q_gamma(:,1:3:end), 2) sum( q_gamma(:,2:3:end), 2) sum( q_gamma(:,3:3:end), 2)];
V_surf = V_gamma + V_wake_plate + V_in - dt_rc_vec;
V_surf1 = V_gamma + V_wake_plate + V_in;


h_V_surf(:,:,i_wake_time) = V_surf;




%% Fluid force evaluation
%%%
%%% dp^* = (Vw^* + Vb^* - dt_r^*)*(τx*dx_Γ^* + τy*dy_Γ^*) + dt_Γ^*
%%%

calc_fluid_force;
