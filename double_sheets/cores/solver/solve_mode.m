%% solve mode

%% Equilibrium point  q0 = q(0), dt_q0 = 0
%%[*] Sheet 1
q_vec = h_X_vec(1:N_q_all,1);
dt_q_vec = 0*q_vec;

%%[*] Sheet 2
q_vec_1 = h_X_vec(N_q_all+1:2*N_q_all,1);
dt_q_vec_1 = 0*q_vec_1;


%% Assemble the viscoelastic vector Qe^(n)
flag_output = 1;                                                            %% Enable stiffness matrix evaluation
theta_a_vec = [ 1 1];                                                     	%% Enable evaluation of ∂q_ε^TD∂q_ε and ∂q_κ^TD∂q_κ

generate_stiff_matrices;                                                    %% Sheet 1
generate_stiff_matrices_1;                                                  %% Sheet 2

theta_a_vec = var_param.theta_a_vec;



%% Linearised stiffness matrix

%%[0-0] Membrane stiffness (sheet 1)
dq_e_Dp_dq_e_global = sparse(N_q_all,N_q_all);
    
for ii = 1:N_element

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};  

    int_dq_e_Dp_dq_e = Qd_eps_mat_i(:,:,ii);    
    dq_e_Dp_dq_e_global(i_vec,i_vec) = dq_e_Dp_dq_e_global(i_vec,i_vec) + squeeze( zeta_m*int_dq_e_Dp_dq_e);
    
end   
    
K_e_m_mat =  dq_e_Dp_dq_e_global;

%%[0-1] Bending stiffness
dq_k_Dp_dq_k_global = sparse(N_q_all,N_q_all);
    
for ii = 1:N_element

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};  

    int_dq_k_Dp_dq_k = Qd_k_mat_i(:,:,ii);    
    dq_k_Dp_dq_k_global(i_vec,i_vec) = dq_k_Dp_dq_k_global(i_vec,i_vec) + squeeze( eta_m*int_dq_k_Dp_dq_k);
    
end   
    
K_e_k_mat =  dq_k_Dp_dq_k_global;





%%[0-1] Membrane stiffness (sheet 2)
dq_e_Dp_dq_e_global_1 = sparse(N_q_all,N_q_all);
    
for ii = 1:N_element

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};  

    int_dq_e_Dp_dq_e = Qd_eps_mat_i_1(:,:,ii);    
    dq_e_Dp_dq_e_global_1(i_vec,i_vec) = dq_e_Dp_dq_e_global_1(i_vec,i_vec) + squeeze( zeta_m*int_dq_e_Dp_dq_e);
    
end   
    
K_e_m_mat_1 =  dq_e_Dp_dq_e_global_1;

%%[0-1] Bending stiffness
dq_k_Dp_dq_k_global_1 = sparse(N_q_all,N_q_all);
    
for ii = 1:N_element

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};  

    int_dq_k_Dp_dq_k = Qd_k_mat_i_1(:,:,ii);    
    dq_k_Dp_dq_k_global_1(i_vec,i_vec) = dq_k_Dp_dq_k_global_1(i_vec,i_vec) + squeeze( eta_m*int_dq_k_Dp_dq_k);
    
end   
    
K_e_k_mat_1 =  dq_k_Dp_dq_k_global_1;


%% Boundary conditions

%%[*] Fixed at zero (sheet 1)
%%[1-0] Displacement boundary condition
if ~isempty( node_r_0)
    i_r = repmat( ( N_qi*(node_r_0 - 1)+1 ).', [ 1 3]) + repmat( 0:2, [ length( node_r_0) 1]);                  %% x,y displacement component indices of the displacement-constrained nodes (z=0 [m])
    i_r = reshape(i_r.',1,[]);
else
    i_r = [];
end

%%[1-1] Slope boundary condition (x direction)
if ~isempty( node_dxr_0)
    i_dx_r = repmat( ( N_qi*(node_dxr_0 - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_0) 1]);           %% dx_r = [1 0 0]^T.
    i_dx_r = reshape(i_dx_r.',1,[]);
else
    i_dx_r = [];
end

%%[1-2] Slope boundary condition (y direction)
if ~isempty( node_dyr_0)
    i_dy_r = repmat( ( N_qi*(node_dyr_0 - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_0) 1]);           %% dy_r = [0 1 0]^T.
    i_dy_r = reshape(i_dy_r.',1,[]);
else
    i_dy_r = [];
end




%%[*] Fixed at zero (sheet 2)
%%[1-0] Displacement boundary condition
if ~isempty( node_r_0_1)
    i_r_1 = repmat( ( N_qi*(node_r_0_1 - 1)+1 ).', [ 1 3]) + repmat( 0:2, [ length( node_r_0_1) 1]);                  %% x,y displacement component indices of the displacement-constrained nodes (z=0 [m])
    i_r_1 = reshape(i_r_1.',1,[]);
else
    i_r_1 = [];
end

%%[1-1] Slope boundary condition (x direction)
if ~isempty( node_dxr_0_1)
    i_dx_r_1 = repmat( ( N_qi*(node_dxr_0_1 - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_0_1) 1]);           %% dx_r = [1 0 0]^T.
    i_dx_r_1 = reshape(i_dx_r_1.',1,[]);
else
    i_dx_r_1 = [];
end

%%[1-2] Slope boundary condition (y direction)
if ~isempty( node_dyr_0_1)
    i_dy_r_1 = repmat( ( N_qi*(node_dyr_0_1 - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_0_1) 1]);           %% dy_r = [0 1 0]^T.
    i_dy_r_1 = reshape(i_dy_r_1.',1,[]);
else
    i_dy_r_1 = [];
end




%%[*] Shared nodal values
%%[1-0] Displacement boundary condition
if ~isempty( node_r_marge)
    i_r_marge = repmat( ( N_qi*(node_r_marge - 1)+1 ).', [ 1 3]) + repmat( 0:2, [ length( node_r_marge) 1]);                  
    i_r_marge = reshape(i_r_marge.',1,[]);
else
    i_r_marge = [];
end

%%[1-1] Slope boundary condition (x direction)
if ~isempty( node_dxr_marge)
    i_dx_r_marge = repmat( ( N_qi*(node_dxr_marge - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_marge) 1]);           
    i_dx_r_marge = reshape(i_dx_r_marge.',1,[]);
else
    i_dx_r_marge = [];
end

%%[1-2] Slope boundary condition (y direction)
if ~isempty( node_dyr_marge)
    i_dy_r_marge = repmat( ( N_qi*(node_dyr_marge - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_marge) 1]);           
    i_dy_r_marge = reshape(i_dy_r_marge.',1,[]);
else
    i_dy_r_marge = [];
end

%%[*] Displacement constraint
%%[*-0] Sheet 1
i_vec = [ i_r i_dx_r i_dy_r];

not_i_vec = (1:N_q_all);
not_i_vec(i_vec) = [];

%%[*-1] Sheet 2
i_vec_1 = [ i_r_1 i_dx_r_1 i_dy_r_1];

not_i_vec_1 = (1:N_q_all);
not_i_vec_1(i_vec_1) = [];


%%[*] Shared nodes
i_vec_marge = [ i_r_marge i_dx_r_marge i_dy_r_marge];



i_vec_all = union( i_vec_1, i_vec_marge).';
if size( i_vec_all, 2) == 1
    i_vec_all = i_vec_all.';
end

not_i_vec_all = (1:N_q_all);
not_i_vec_all(i_vec_all) = [];

%% Apply the boundary conditions

%%[*] Mass matrix
M_mat_BC = M_global;
M_mat_BC_1 = M_mat_BC;
%%[*] Shared-node boundary condition
%%[*] Sheet 1
U_M_mat = 0*M_mat_BC;
U_M_mat(i_vec_marge,:) = M_mat_BC_1(i_vec_marge,:);
M_mat_BC(i_vec_marge,i_vec_marge) = M_mat_BC(i_vec_marge,i_vec_marge) + M_mat_BC_1(i_vec_marge,i_vec_marge);

%%[*] Sheet 2
L_M_mat = 0*M_mat_BC_1;
L_M_mat(:,i_vec_marge) = M_mat_BC_1(:,i_vec_marge);

%%[*] Fixed-nodal-value boundary condition
M_global_BC = [ M_mat_BC  	U_M_mat;
                L_M_mat    	M_mat_BC_1];
            
M_global_BC([ i_vec N_q_all+i_vec_all],:) = [];
M_global_BC(:,[ i_vec N_q_all+i_vec_all]) = [];



%%[*] Stiffness matrix
K_e_m_mat_BC = K_e_m_mat;
K_e_k_mat_BC = K_e_k_mat;

K_e_m_mat_BC_1 = K_e_m_mat_1;
K_e_k_mat_BC_1 = K_e_k_mat_1;

%%[*] Shared-node boundary condition
%%[*] Sheet 1
U_e_m_mat = 0*K_e_m_mat_BC;
U_e_m_mat(i_vec_marge,:) = K_e_m_mat_BC_1(i_vec_marge,:);
K_e_m_mat_BC(i_vec_marge,i_vec_marge) = K_e_m_mat_BC(i_vec_marge,i_vec_marge) + K_e_m_mat_BC_1(i_vec_marge,i_vec_marge);

U_e_k_mat = 0*K_e_k_mat_BC;
U_e_k_mat(i_vec_marge,:) = K_e_k_mat_BC_1(i_vec_marge,:);
K_e_k_mat_BC(i_vec_marge,i_vec_marge) = K_e_k_mat_BC(i_vec_marge,i_vec_marge) + K_e_k_mat_BC_1(i_vec_marge,i_vec_marge);

%%[*] Sheet 2
L_e_m_mat = 0*K_e_m_mat_BC_1;
L_e_m_mat(:,i_vec_marge) = K_e_m_mat_BC_1(:,i_vec_marge);

L_e_k_mat = 0*K_e_k_mat_BC_1;
L_e_k_mat(:,i_vec_marge) = K_e_k_mat_BC_1(:,i_vec_marge);


%%[*] Fixed-nodal-value boundary condition
%%[*] Sheet 1
K_e_m_global_BC = [ K_e_m_mat_BC    U_e_m_mat;
                    L_e_m_mat       K_e_m_mat_BC_1];
K_e_k_global_BC = [ K_e_k_mat_BC	U_e_k_mat;
                    L_e_k_mat       K_e_k_mat_BC_1];
                
K_e_m_global_BC([ i_vec N_q_all+i_vec_all],:) = [];               
K_e_m_global_BC(:,[ i_vec N_q_all+i_vec_all]) = [];
K_e_k_global_BC([ i_vec N_q_all+i_vec_all],:) = [];
K_e_k_global_BC(:,[ i_vec N_q_all+i_vec_all]) = [];


%% modal analysis

%%[2-0] Compute the mode_num smallest eigenvalues
[ Phi_dq_mat, omega_a2] = eigs( (mu_m*M_global_BC)\(K_e_m_global_BC + K_e_k_global_BC), mode_num, 'SM');   


%%[2-1] Nondimensional natural frequency: ω^* := L/U_in*ω [-]
omega_a = sqrt( diag( omega_a2));

%%[2-1] Nondimensional displacement eigenmode: corresponds to Δq(t)
Phi_dq_mat_BC = zeros(2*N_q_all,mode_num);
Phi_dq_mat_BC([ not_i_vec N_q_all+not_i_vec_all],:) = Phi_dq_mat;

Phi_dq_mat_BC(N_q_all+i_vec_marge,:) = Phi_dq_mat_BC(i_vec_marge,:);

%%[2-2] Nondimensional position eigenmode: corresponds to q(t) = q0 + Δq(t)
Phi_q_mat_BC = repmat( [ q_vec; q_vec_1], [ 1 mode_num]) + Phi_dq_mat_BC;



%% Modal damping ratio

%%[3-0] Modal damping ratio [-]
zeta_n = theta_a_vec(1)*omega_a/2;                      



disp( 'Nondimensional natural frequency: ')
disp( omega_a)
disp( 'Modal damping ratio: ')
disp( zeta_n)


