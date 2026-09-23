function [ out, out1] = new_X_func_FAST( X_vec, m_global_struct, qf_global_struct, dq_qe_global_struct, qe_global_struct, qd_global_struct, var_param, stage, out1)


%% [0] Extract variables
node_r_0 = var_param.node_r_0;
node_dxr_0 = var_param.node_dxr_0;
node_dyr_0 = var_param.node_dyr_0;
node_dxr_theta_c = var_param.node_dxr_theta_c;
node_dyr_theta_c = var_param.node_dyr_theta_c;
coordinates = var_param.coordinates;

d_t = var_param.d_t;
alpha_v = var_param.alpha_v;
theta_a= var_param.theta_a ;

N_qi = var_param.N_qi; 
N_q_all = var_param.N_q_all;


%%[*] parameters
M_global = m_global_struct.M_global;                                        %% Mass matrix      : M (mass + added mass from the sheet itself [ Madd_1 Madd_2]*d_t^2[ q_1^T q_2^T]^T -> Madd_1*d_t^2 q_1)
Qf_global = qf_global_struct.Qf_global;                                     %% Body force        : Q_f
dq_Qe_global = dq_qe_global_struct.dq_Qe_global;                            %% dq_Qe(q(n))
Qe_global = qe_global_struct.Qe_global;                                     %% Elastic force     : Q_e
Qd_global = qd_global_struct.Qd_global;                                     %% Damping force     : Q_d




%% [1] Boundary conditions

%%[*] Fixed at zero
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



%%[*] Constrained to a common value
Fmat = speye( N_q_all); 

%%[1-3] Slope boundary condition (x direction)
if ~isempty( node_dxr_theta_c)    
    
    %%[1-3-0] Node index of the reference value [-]
    i_dx_r_theta_c1 = repmat( ( N_qi*(node_dxr_theta_c(1) - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_theta_c(1)) 1]);	
    i_dx_r_theta_c1 = reshape( i_dx_r_theta_c1.', 1, []);    
    
    %%[1-3-1] Node index being modified [-]
    i_dx_r_theta_c = repmat( ( N_qi*(node_dxr_theta_c(2:end) - 1)+4 ).', [ 1 3]) + repmat( 0:2, [ length( node_dxr_theta_c(2:end)) 1]);	
    i_dx_r_theta_c = reshape( i_dx_r_theta_c.', 1, []);    


    Fmat(i_dx_r_theta_c,:) = 0;
    Fmat(i_dx_r_theta_c,i_dx_r_theta_c1) = repmat( diag( [ 1 0 1]), [ length( node_dxr_theta_c(2:end)) 1]);
    
    i_dx_r_theta_c1_2 = i_dx_r_theta_c1(2);
else
    i_dx_r_theta_c1 = [];
    i_dx_r_theta_c1_2 = [];
    i_dx_r_theta_c = [];
end

%%[1-4] Slope boundary condition (y direction)
if ~isempty( node_dyr_theta_c)    
    
    %%[1-4-0] Node index of the reference value [-]
    i_dy_r_theta_c1 = repmat( ( N_qi*(node_dyr_theta_c(1) - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_theta_c(1)) 1]);	
    i_dy_r_theta_c1 = reshape( i_dy_r_theta_c1.', 1, []);    
    
    %%[1-4-1] Node index being modified [-]
    i_dy_r_theta_c = repmat( ( N_qi*(node_dyr_theta_c(2:end) - 1)+7 ).', [ 1 3]) + repmat( 0:2, [ length( node_dyr_theta_c(2:end)) 1]);	
    i_dy_r_theta_c = reshape( i_dy_r_theta_c.', 1, []);    


    Fmat(i_dy_r_theta_c,:) = 0;
    Fmat(i_dy_r_theta_c,i_dy_r_theta_c1) = repmat( diag( [ 0 1 1]), [ length( node_dyr_theta_c(2:end)) 1]);
    
    i_dy_r_theta_c1_1 = i_dy_r_theta_c1(1);
else
    i_dy_r_theta_c1 = [];
    i_dy_r_theta_c1_1 = [];
    i_dy_r_theta_c = [];
end

i_vec = [ i_r i_dx_r i_dy_r i_dx_r_theta_c1_2 i_dx_r_theta_c i_dy_r_theta_c1_1 i_dy_r_theta_c];

%% [2] Evaluate the acceleration dtt_q

Q_global = (Qf_global - Qe_global);                                                                     %% Internal + external force terms

M_global = M_global*Fmat;
if ~isempty( i_dx_r_theta_c1)
    M_global_mat = permute( reshape( permute( full( M_global([i_dx_r_theta_c1 i_dx_r_theta_c],:)), [1 3 2]), 3, [], N_q_all), [1 3 2]);
    M_global(i_dx_r_theta_c1,:) = sum( M_global_mat, 3);
end
if ~isempty( i_dy_r_theta_c1)
    M_global_mat = permute( reshape( permute( full( M_global([i_dy_r_theta_c1 i_dy_r_theta_c],:)), [1 3 2]), 3, [], N_q_all), [1 3 2]);
    M_global(i_dy_r_theta_c1,:) = sum( M_global_mat, 3);
end
M_global(i_vec,:) = [];
M_global(:,i_vec) = [];
if ~isempty( i_dx_r_theta_c1)
    Q_global_mat = permute( reshape( permute( Q_global([i_dx_r_theta_c1 i_dx_r_theta_c]), [1 3 2]), 3, [], 1), [1 3 2]);
    Q_global(i_dx_r_theta_c1) = sum( Q_global_mat, 3);
end
if ~isempty( i_dy_r_theta_c1)
    Q_global_mat = permute( reshape( permute( Q_global([i_dy_r_theta_c1 i_dy_r_theta_c]), [1 3 2]), 3, [], 1), [1 3 2]);
    Q_global(i_dy_r_theta_c1) = sum( Q_global_mat, 3);
end
Q_global(i_vec) = []; 

Qd_global = Qd_global*Fmat;
if ~isempty( i_dx_r_theta_c1)
    Qd_global_mat = permute( reshape( permute( full( Qd_global([i_dx_r_theta_c1 i_dx_r_theta_c],:)), [1 3 2]), 3, [],N_q_all), [1 3 2]);
    Qd_global(i_dx_r_theta_c1,:) = sum( Qd_global_mat, 3);
end
if ~isempty( i_dy_r_theta_c1)
    Qd_global_mat = permute( reshape( permute( full( Qd_global([i_dy_r_theta_c1 i_dy_r_theta_c],:)), [1 3 2]), 3, [], N_q_all), [1 3 2]);
    Qd_global(i_dy_r_theta_c1,:) = sum( Qd_global_mat, 3);
end
Qd_global(i_vec,:) = [];
Qd_global(:,i_vec) = [];
dq_Qe_global = dq_Qe_global*Fmat;
if ~isempty( i_dx_r_theta_c1)
    dq_Qe_global_mat = permute( reshape( permute( full( dq_Qe_global([i_dx_r_theta_c1 i_dx_r_theta_c],:)), [1 3 2]), 3, [], N_q_all), [1 3 2]);
    dq_Qe_global(i_dx_r_theta_c1,:) = sum( dq_Qe_global_mat, 3);
end
if ~isempty( i_dy_r_theta_c1)
    dq_Qe_global_mat = permute( reshape( permute( full( dq_Qe_global([i_dy_r_theta_c1 i_dy_r_theta_c],:)), [1 3 2]), 3, [], N_q_all), [1 3 2]);
    dq_Qe_global(i_dy_r_theta_c1,:) = sum( dq_Qe_global_mat, 3);
end
dq_Qe_global(i_vec,:) = [];
dq_Qe_global(:,i_vec) = [];


eye_mat = speye( N_q_all);
eye_mat(i_vec,:) = [];
eye_mat(:,i_vec) = [];
zero_mat = sparse( N_q_all-length( i_vec),  N_q_all-length( i_vec));
if stage == 0                                                                                           %% Reuse the value computed on the first pass.
    C_damp = (theta_a == 0)*2 + ~(theta_a == 0)*1;                                                      %% Semi-implicit when damping is present
    D_matrix = [ eye_mat                                    zero_mat;
                 (Qd_global + C_damp*d_t/2*dq_Qe_global)  	M_global];                                  %% Reduces to the identity matrix at zero damping.
     out1.D_matrix = D_matrix;
else
    D_matrix = out1.D_matrix;
end
X2_matrix = [ zero_mat      eye_mat;
              zero_mat      zero_mat];
A_mat1 = D_matrix - alpha_v*d_t*X2_matrix;
A_mat2 = D_matrix + (1 - alpha_v)*d_t*X2_matrix;




%%[2-0] At the clamped nodes dt_q = 0, dtt_q = 0
not_i_vec = (1:N_q_all);
not_i_vec(i_vec) = [];


%% [3] State vector update

if stage == 0                                                                                           %% Reuse the value computed on the first pass.
    out1.A1_A2_Xn = A_mat1\( A_mat2*X_vec([ not_i_vec N_q_all+not_i_vec]) );
end
out_0 = out1.A1_A2_Xn + A_mat1\( d_t*[ zero_mat(:,1);
                                       Q_global]);
                    
out = X_vec;
out([ not_i_vec N_q_all+not_i_vec]) = out_0;
out = blkdiag( Fmat, Fmat)*out;

end