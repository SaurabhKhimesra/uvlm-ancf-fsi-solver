%% Stiffness matrix assembly

%%[*] Gauss-Legendre weights as a matrix: 1:1~N_q, 2:1~N_q, 3:ξ quadrature points, 4:η quadrature points, 5:component index used when differencing ∂qκ
ww_mat = repmat( permute( w_vec.'*w_vec, [ 3 4 1 2]), [ N_q 1 1 1]);
ww_mat2 = repmat( permute( w_vec.'*w_vec, [ 3 4 1 2]), [ N_q N_q 1 1]);
ww_mat3 = repmat( permute( w_vec.'*w_vec, [ 3 4 1 2]), [ 1 1 1 1]);
ww_mat4 = repmat( permute( w_vec, [ 3 4 1 2]), [ N_q N_q 1 1]);



%%[*] Stiffness copied to the Gauss-Legendre quadrature points only: 1:1~N_q, 2:1~N_q, 3:ξ quadrature points, 4:η quadrature points, 5:component index used when differencing ∂qκ
Dp_mat2 = repmat( Dp_mat, [ 1 1 length( p_vec) length( p_vec)]);


%% [0] Normal strain (Qeε = t∫_A (∂q_ε)^t*Dm*ε dA ∈ R^36)
if flag_output == 1     
    
    int_e_Dp_e_vec = zeros(1,N_element);

    Qe_eps_vec_i = zeros(N_q,N_element);
    dq_Qe_eps_mat_i = zeros(N_q,N_q,N_element);                           	%% ∫_A (dqdq_ε^T*D*ε + dq_ε^T*D*dq_ε) dA 
    Qd_eps_mat_i = zeros(N_q,N_q,N_element);                            	%% ∫_A dq_ε^T*D*dq_ε dA (dq_ε*dt_q = dt_ε)    

 
    for ii = 1:N_element

        dL = dL_vec(ii);                        %% Length of element ii [-]
        dW = dW_vec(ii);                        %% Width of element ii [-]


        %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
        %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
        i_vec = i_vec_v{ii};
        q_i_vec = q_vec(i_vec);

        dx_n_sc_struct = dx_n_Sc_struct(ii);       


        [ dq_eps_Dp_eps_vec, dq_eps_Dp_dq_eps_mat, dqdq_eps_Dp_eps_P_dq_eps_Dp_dq_eps_mat, eps_Dp_eps] = dq_eps_Dm_eps( q_i_vec, Dp_mat2, dx_n_sc_struct, p_vec);
        int_dq_e_Dp_e = dL*dW/4*ww_mat.*dq_eps_Dp_eps_vec;
        int_dq_e_Dp_e = sum( sum( int_dq_e_Dp_e, 4), 3).';        

        Qe_eps_vec_i(:,ii) = int_dq_e_Dp_e;                                 %% External force matrix of element ii

        int_e_Dp_e = dL*dW/4*ww_mat3.*eps_Dp_eps;
        int_e_Dp_e_vec(ii) = sum( sum( int_e_Dp_e, 4), 3).';
        

        int_dqdq_e_Dp_e_P_dq_e_Dp_dq_e = dL*dW/4*ww_mat2.*dqdq_eps_Dp_eps_P_dq_eps_Dp_dq_eps_mat;
        int_dqdq_e_Dp_e_P_dq_e_Dp_dq_e = sum( sum( int_dqdq_e_Dp_e_P_dq_e_Dp_dq_e, 4), 3).';

        dq_Qe_eps_mat_i(:,:,ii) = int_dqdq_e_Dp_e_P_dq_e_Dp_dq_e;       	%% External force matrix of element ii

        if theta_a ~= 0                                                     %% Evaluated only when structural damping is present (saves computation)
            int_dq_e_Dp_dq_e = dL*dW/4*ww_mat2.*dq_eps_Dp_dq_eps_mat;
            int_dq_e_Dp_dq_e = sum( sum( int_dq_e_Dp_dq_e, 4), 3).';

            Qd_eps_mat_i(:,:,ii) = int_dq_e_Dp_dq_e;                        %% External force matrix of element ii
        end    


    end
end

%% [1] Bending strain (Qeε = ∫_A (∂q_κ)^t*Dκ*κ dA ∈ R^36)
int_k_Dp_k_vec = zeros(1,N_element);

Qe_k_vec_i = zeros(N_q,N_element);
Qd_k_mat_i = zeros(N_q,N_q,N_element);                                      %% ∫_A dq_κ^T*D*dq_κ dA (dq_κ*dt_q = dt_κ)


%%[1-1] Perturbation used when differencing ∂qκ: 1:1~N_q, 2:1~N_q, 3:ξ quadrature points, 4:η quadrature points, 5:component index used when differencing ∂qκ
% eye_N_q = repmat( eye(N_q), [ 1 1 length( p_vec)  length( p_vec)]);
% d_q = 1e-10;
% dq_vec = d_q*eye_N_q;
% dq_vec = permute( dq_vec, [ 5 1 3 4 2]);
% dq_vec = reshape( dq_vec, 3, [], length( p_vec), length( p_vec), N_q);	 %% Match the shape-function row count (for the product). Dim 1: coordinate component (x,y,z), dim 2: direction (r,dx_r,dy_r), dims 3-4: Gauss quadrature points, dim 5: component index used when differencing ∂qκ



for ii = 1:N_element
    
    dL = dL_vec(ii);                %% Length of element ii [-]
    dW = dW_vec(ii);                %% Width of element ii [-]
    
    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    q_i_vec = q_vec(i_vec);    
    
    dx_n_sc_struct = dx_n_Sc_struct(ii);            

    [ dq_k_Dp_k_vec, dq_k_Dp_dq_k_mat, k_Dp_k] = dq_k_Dk_k_FAST( q_i_vec, Dp_mat2, dx_n_sc_struct, p_vec, theta_a);
    
    int_k_Dp_k = dL*dW/4*ww_mat3.*k_Dp_k;
    int_k_Dp_k = sum( sum( int_k_Dp_k, 4), 3).'; 
    int_k_Dp_k_vec(ii) = int_k_Dp_k;
    
    int_dq_k_Dp_k = dL*dW/4*ww_mat.*dq_k_Dp_k_vec;
    int_dq_k_Dp_k = sum( sum( int_dq_k_Dp_k, 4), 3).';  
    Qe_k_vec_i(:,ii) = int_dq_k_Dp_k;                                       %% External force matrix of element ii
    
    if theta_a ~= 0                                                         %% Evaluated only when structural damping is present (saves computation)
        int_dq_k_Dp_dq_k = dL*dW/4*ww_mat2.*dq_k_Dp_dq_k_mat;
        int_dq_k_Dp_dq_k = sum( sum( int_dq_k_Dp_dq_k, 4), 3).';
        
        Qd_k_mat_i(:,:,ii) = int_dq_k_Dp_dq_k;
    end 

end





%% [2] Rotational damper [-] (Qdθ = -∫_ΔW dq_(n/||n||)^T*dx_r/||dx_r||*∂q_θ0y dy ∈ R^36x36)
 

Qd_theta_mat_i = zeros(N_q,N_q,N_element);                                  %% ∫ dq(n/|n|)^T*(dx_r/|dx_r|)*dq_θ_0y dy 
dt_theta_0y = zeros(size( element_C_theta));

for ii = element_C_theta

    dW = dW_vec(ii);                %% Width of element ii [-]

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    q_i_vec = q_vec(i_vec); 
    dt_q_i_vec = dt_q_vec(i_vec);   
    

    dx_n_sc_struct = dx_n_Sc_struct(ii);     

    [ q_d_theta_mat, dq_theta_0y_v] = dq_n_dx_r_dq_theta_0y( q_i_vec, dx_n_sc_struct, p_vec);

    int_q_d_theta_mat = dW/2*ww_mat4.*q_d_theta_mat;                        %% Constant in x (x=0 [-]), integrated in y
    int_q_d_theta_mat = sum( sum( int_q_d_theta_mat, 4), 3);   
    
    int_dq_theta_0y_v = dW/2*ww_mat4(1,:,:,:).*dq_theta_0y_v;
    int_dq_theta_0y_v = sum( sum( int_dq_theta_0y_v, 4), 3);
    

    Qd_theta_mat_i(:,:,ii) = int_q_d_theta_mat;
    
    dt_theta_0y(ii) = int_dq_theta_0y_v*dt_q_i_vec;
end





%% [3] Second term of the rotational inertia matrix [-] (-∫_ΔW dq_(n/||n||)^T*dx_r/||dx_r||*dt∂q_θ0y dy ∈ R^36x36)
 

J2_mat_i = zeros(N_q,N_q,N_element);                                  %% ∫ dq(n/|n|)^T*(dx_r/|dx_r|)*dq_θ_0y dy 


for ii = element_C_theta

    dW = dW_vec(ii);                %% Width of element ii [-]

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    q_i_vec = q_vec(i_vec);  
    dt_q_i_vec = dt_q_vec(i_vec);    

    dx_n_sc_struct = dx_n_Sc_struct(ii);     

    J2_mat = dq_n_dx_r_dtdq_theta_0y( q_i_vec, dt_q_i_vec, dx_n_sc_struct, p_vec);

    int_J2_mat = dW/2*ww_mat4.*J2_mat;                        %% Constant in x (x=0 [-]), integrated in y
    int_J2_mat = sum( sum( int_J2_mat, 4), 3);           

    J2_mat_i(:,:,ii) = int_J2_mat;
end





%% [4] Global matrix assembly

%%[3-0] Elastic force
Qe_global = zeros(N_q_all,1);
Qk_global = zeros(N_q_all,1);
for ii = 1:N_element

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};    

    if flag_output == 1 
        Qe_global(i_vec,1) = Qe_global(i_vec,1) + zeta_m*Qe_eps_vec_i(:,ii);
    end
    Qk_global(i_vec,1) = Qk_global(i_vec,1) + eta_m*Qe_k_vec_i(:,ii);
end


%%[4-1] Stiffness matrix + structural damping
dq_Qe_global = sparse(N_q_all,N_q_all);
Qd_global = sparse(N_q_all,N_q_all);

dq_Qe_global_vec = zeros(1,length( dq_Qe_eps_mat_i(:)));
Qd_global_vec = zeros(1,length( Qd_eps_mat_i(:)));
if flag_output == 1 
    
    for ii = 1:N_element

        %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
        %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
        i_vec = i_vec_v{ii};        
        
        dq_Qe_global_vec_i = zeta_m*dq_Qe_eps_mat_i(:,:,ii);
        dq_Qe_global_vec((ii - 1)*N_q^2+(1:N_q^2)) = dq_Qe_global_vec_i(:).';
        
%         dq_Qe_global(i_vec,i_vec) = dq_Qe_global(i_vec,i_vec) + squeeze( zeta_m*dq_Qe_eps_mat_i(:,:,ii));

        if theta_a ~= 0                                                         %% Evaluated only when structural damping is present (saves computation)
            
%             Qd_global(i_vec,i_vec) = Qd_global(i_vec,i_vec) + squeeze( theta_a*zeta_m*Qd_eps_mat_i(:,:,ii) +  theta_a*eta_m*Qd_k_mat_i(:,:,ii));
            
            Qd_global_vec_i = theta_a*zeta_m*Qd_eps_mat_i(:,:,ii) +  theta_a*eta_m*Qd_k_mat_i(:,:,ii);
            Qd_global_vec((ii - 1)*N_q^2+(1:N_q^2)) = Qd_global_vec_i(:).';
        end     
    end

    dq_Qe_global = setsparse( dq_Qe_global, I_vec, J_vec, dq_Qe_global_vec, @plus);
    if theta_a ~= 0                                                         %% Evaluated only when structural damping is present (saves computation)
        Qd_global = setsparse( Qd_global, I_vec, J_vec, Qd_global_vec, @plus);
    end
    
end


%%[4-2] Rotational damping
Qd_theta_global = sparse(N_q_all,N_q_all);

    
for ii = element_C_theta

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};

    Qd_theta_global(i_vec,i_vec) = Qd_theta_global(i_vec,i_vec) + squeeze( C_theta_a*Qd_theta_mat_i(:,:,ii));
end


%%[4-3] Rotational inertia: J1*dt^2_q + J2*dt_q
J_global_1 = sparse(N_q_all,N_q_all);    
for ii = element_C_theta

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};  

    J_global_1(i_vec,i_vec) = J_global_1(i_vec,i_vec) + squeeze( J_a*Qd_theta_mat_i(:,:,ii));
end

J_global_2 = sparse(N_q_all,N_q_all);    
for ii = element_C_theta

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii}; 

    J_global_2(i_vec,i_vec) = J_global_2(i_vec,i_vec) + squeeze( J_a*J2_mat_i(:,:,ii));
end

    

