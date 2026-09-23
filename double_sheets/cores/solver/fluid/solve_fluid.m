

q_vec = X_vec(1:N_q_all);
dt_q_vec = X_vec(2*N_q_all+1:3*N_q_all);


q_vec_1 = X_vec(N_q_all+1:2*N_q_all,1);
dt_q_vec_1 = X_vec(3*N_q_all+1:4*N_q_all,1);


dt_q_vec_all = [ dt_q_vec;
                 dt_q_vec_1];


%% Panel node point evaluation [m]
                         
generate_r_panel;

%%[0] Time history of the panel node points [-]
%%% 1:element index, 2:coordinate component, 3:time            (upper surface)
h_r_panel_vec(:,:,i_wake_time) = [      r_panel_vec_1;
                                        r_panel_vec_2;
                                        r_panel_vec_3;
                                        r_panel_vec_4];
 %%% 1:element index, 2:coordinate component, 3:time            (lower surface)
h_r_panel_vec_1(:,:,i_wake_time) = [	r_panel_vec_1_1;
                                        r_panel_vec_2_1;
                                        r_panel_vec_3_1;
                                        r_panel_vec_4_1];
                               

                            
%% Collocation point evaluation [m]
%%% 1:element index, 2:coordinate component, 3:time    (upper surface)
h_rcol_vec(:,:,i_wake_time) = rc_vec;

%%% 1:element index, 2:coordinate component, 3:time    (lower surface)
h_rcol_vec_1(:,:,i_wake_time) = rc_vec_1;





%% Unit normal vector evaluation

%%[1-0] n_vec
generate_dt_n_vec;




%%% 1:element index, 2:coordinate component, 3:time    (upper surface)
h_n_vec(:,:,i_wake_time) = n_vec_i;
h_dt_n_vec(:,:,i_wake_time) = dt_n_vec_i;

%%% 1:element index, 2:coordinate component, 3:time    (lower surface)
h_n_vec_1(:,:,i_wake_time) = n_vec_i_1;
h_dt_n_vec_1(:,:,i_wake_time) = dt_n_vec_i_1;


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
q1234_mat = generate_q1234_mat( rc_vec_all, r_panel_vec_1_all, r_panel_vec_2_all, r_panel_vec_3_all, r_panel_vec_4_all, var_param, 1);                                    


%%[2-1] Generation of influence coefficient matrix
n_vec_i_mat = repmat( n_vec_i_all, [ 1 size( n_vec_i_all, 1)]);
A_mat = inner_mat( q1234_mat, n_vec_i_mat);
A_mat = A_mat(:,1:3:end);                                                   %% Copying only the 3 components inside inner_mat is unnecessary.


%% Initial wake parameters [-]
if ~exist( 'Gamma_wake_all', 'var')    
    
    %% Upper surface
    old_Gamma = zeros(N_element,1);    
    Gamma = zeros(N_element,1);   
    %% Lower surface
    old_Gamma_1 = old_Gamma;    
    Gamma_1 = Gamma;

    Gamma_all = [   Gamma;
                    Gamma_1];
    Gamma_trail = zeros(Ny,1);
    Gamma_trail_1 = Gamma_trail;
end




%% Wake generation

generate_wake;


%% Wake-induced velocity [-]

[ V_wake_plate, q1234_wake_mat] = V_wake_func( rc_vec_all, r_wake_1_all, r_wake_2_all, r_wake_3_all, r_wake_4_all, Gamma_wake_all, var_param, 1);
    
%%[*] --------------------------------------------------------------------added (18:36 2017/03/30)
%%[*] Velocity induced on the plate by the vortex sheet with the freshly shed circulation removed [-]
%%[*] --------------------------------------------------------------------
n_vec_i_mat2 = repmat( n_vec_i_all, [ 1 size( r_wake_1_all, 1)]);
q1234_wake_mat_n = inner_mat( q1234_wake_mat, n_vec_i_mat2);
q1234_wake_mat_n = q1234_wake_mat_n(:,1:3:end);                        	%% Copying only the 3 components inside inner_mat is unnecessary.
%%[*-0] (sheet 1)
Gamma_wake2 = zeros(size( r_wake_1, 1),1);
Gamma_wake2(Ny+1:end) = Gamma_wake(Ny+1:end);
%%[*-1] (sheet 2)
Gamma_wake2_1 = zeros(size( r_wake_1_1, 1),1);
Gamma_wake2_1(Ny+1:end) = Gamma_wake_1(Ny+1:end);
%%[*-2] Merge
Gamma_wake2_all = [ Gamma_wake2;
                    Gamma_wake2_1];
V_wake_plate2_n = sum( q1234_wake_mat_n.*(ones(2*N_element,1)*Gamma_wake2_all.'), 2);




%% Normal inlet and body velocity [-]


V_normal = sum( (dt_rc_vec_all - V_in_all - V_wake_plate).*n_vec_i_all, 2);






%% Circulation evaluation

Gamma_all = A_mat\V_normal;
Gamma = Gamma_all(1:end/2);
Gamma_1 = Gamma_all(end/2+1:end);

%%% Circulation of the trailing-edge panel (used for the wake circulation via the unsteady Kutta condition: Dt_Γ=V^T*∇Γ+∂t_Γ=0 -> vorticity convection)
Gamma_trail = old_Gamma(end-Ny+1:end);    
Gamma_trail_1 = old_Gamma_1(end-Ny+1:end);    





%% Update the wake-induced velocity
%%[*-0] (sheet 1)
Gamma_wake(1:Ny) = Gamma_trail;
Gamma_wake_no_trail = Gamma_wake;
Gamma_wake_no_trail(1:Ny) = 0;
%%[*-1] (sheet 2)
Gamma_wake_1(1:Ny) = Gamma_trail_1;
Gamma_wake_no_trail_1 = Gamma_wake_1;
Gamma_wake_no_trail_1(1:Ny) = 0;
%%[*-2] Merge
Gamma_wake_no_trail_all = [ Gamma_wake_no_trail;
                            Gamma_wake_no_trail_1];

r_wake_trail_1_all = [  r_wake_1(1:Ny,:);
                        r_wake_1_1(1:Ny,:)];
r_wake_trail_2_all = [  r_wake_2(1:Ny,:);
                        r_wake_2_1(1:Ny,:)];
r_wake_trail_3_all = [  r_wake_3(1:Ny,:);
                        r_wake_3_1(1:Ny,:)];
r_wake_trail_4_all = [  r_wake_4(1:Ny,:);
                        r_wake_4_1(1:Ny,:)];
                    
Gamma_trail_all = [ Gamma_trail;
                    Gamma_trail_1];

V_wake_plate_trail = V_wake_func( rc_vec_all, r_wake_trail_1_all, r_wake_trail_2_all, r_wake_trail_3_all, r_wake_trail_4_all, Gamma_trail_all, var_param, 1);
q_gamma_wake_no_trail = q1234_wake_mat.*( ones(2*N_element,1)*kron( Gamma_wake_no_trail_all.', ones(1,3)) );    
V_wake_plate_no_trail = [ sum( q_gamma_wake_no_trail(:,1:3:end), 2) sum( q_gamma_wake_no_trail(:,2:3:end), 2) sum( q_gamma_wake_no_trail(:,3:3:end), 2)];
V_wake_plate = V_wake_plate_trail + V_wake_plate_no_trail;

%%[*-3] Refresh the circulation of every wake panel for the fluid force that follows (brings the wake just behind the trailing edge up to date)
Gamma_wake_all = [  Gamma_wake;
                	Gamma_wake_1];
             


%% Surface velocity distribution (at the collocation points)

q_gamma = q1234_mat.*( ones(2*N_element,1)*kron( Gamma_all.', ones(1,3)) );    
V_gamma = [ sum( q_gamma(:,1:3:end), 2) sum( q_gamma(:,2:3:end), 2) sum( q_gamma(:,3:3:end), 2)];
V_surf = V_gamma + V_wake_plate + V_in_all - dt_rc_vec_all;
V_surf1 = V_gamma + V_wake_plate + V_in_all;


h_V_surf(:,:,i_wake_time) = V_surf;




%% Fluid force evaluation
%%%
%%% dp^* = (Vw^* + Vb^* - dt_r^*)*(τx*dx_Γ^* + τy*dy_Γ^*) + dt_Γ^*
%%%

calc_fluid_force;





%% Update the value held from the previous step

%%[*] The previous step's circulation is used so the Kutta condition is satisfied.
old_Gamma = Gamma;  
old_Gamma_1 = Gamma_1;  
