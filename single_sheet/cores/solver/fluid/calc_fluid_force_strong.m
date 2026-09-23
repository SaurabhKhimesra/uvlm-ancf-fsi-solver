%%[*] Fluid force vector: [p] = Mf1*dt^2_q + (Mf2_1*(dt_r - Vin - Vwake)^T*dt_ni + Mf2_2)

%% [2-0] Fluid force vector (Mf2_2)
%%[*] Nodal values used for the linear interpolation of the fluid force
dp_nvec_i = zeros(N_element,3,3);        
dp_nvec_i(:,:,2) = ( (dp_lift1 + Mf2_vec1)*ones(1,3) ).*n_vec_i;
dp_nvec_i(Ny+1:end,:,1) = dp_nvec_i(1:end-Ny,:,2);
dp_nvec_i(1:end-Ny,:,3) = dp_nvec_i(Ny+1:end,:,2);
        

Qf_p_vec_i = zeros(N_q,N_element);
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% Length of element ii [-]
    dW = dW_vec(ii);                        %% Width of element ii [-]
            
    int_StF = zeros(N_q,1);
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre quadrature
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] Linear interpolation of the fluid force
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 3 2 1]);
        dp_ni_interp = sum( p_interp_vec(:,ones(1,3),:).*dp_nvec_i(ii,:,:), 3);
        
        
        i_eta_a = 1;
        for eta_a = p_vec           
            
            
            %%[*] Pressure acts along the normal vector, so the sign is positive.
            StF = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_ni_interp.';
            int_StF = int_StF + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StF;
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_vec_i(:,ii) = int_StF;            %% External force matrix of element ii (per unit area)
end



%% [2-1] Fluid force matrix ( -(τx*dxΓ + τy*dyΓ)*dt_rc )
%%[*] Nodal values used for the linear interpolation of the fluid force
dp_lift2_nvec_i = zeros(3,3,N_element,3);      
n_vec_mat = permute( n_vec_i, [ 2 3 1]);
dp_lift2 = permute( dp_lift2, [3 2 1]);
niT_dp_lift2 = mntimes2( n_vec_mat, dp_lift2);

dp_lift2_nvec_i(:,:,:,2) = niT_dp_lift2;
dp_lift2_nvec_i(:,:,Ny+1:end,1) = dp_lift2_nvec_i(:,:,1:end-Ny,2);
dp_lift2_nvec_i(:,:,1:end-Ny,3) = dp_lift2_nvec_i(:,:,Ny+1:end,2);

Qf_p_lift2_mat_i = zeros(N_q,3*N_element,N_element);
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% Length of element ii [-]
    dW = dW_vec(ii);                        %% Width of element ii [-]
    
   
            
    int_StM = zeros(N_q,3*N_element);
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre quadrature
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] Linear interpolation of the fluid force
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 4 2 3 1]);
        dp_lift2_nvec_i_interp = sum( p_interp_vec(ones(1,3),ones(1,3),:,:).*dp_lift2_nvec_i(:,:,ii,:), 4);
        
        i_eta_a = 1;
        dp_lift2_nvec_i_interp_v = zeros(3,3*N_element);
        dp_lift2_nvec_i_interp_v(:,3*ii-2:3*ii) = dp_lift2_nvec_i_interp;
        for eta_a = p_vec
            
            %%[*] Pressure acts along the normal vector, so the sign is positive.
            StM = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_lift2_nvec_i_interp_v;
            int_StM = int_StM + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM;
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_lift2_mat_i(:,:,ii) = int_StM;            %% External force matrix of element ii (per unit area)
end


%% [2-2] Fluid force matrix (Mf2_1)
%%[*] Nodal values used for the linear interpolation of the fluid force
dp_add_nvec0_i = zeros(3,N_element,N_element,3);      
n_vec_mat = permute( n_vec_i, [ 2 3 1]);
Mf2_mat_i = permute( Mf2_mat, [3 2 1]);
niT_Mf2_mat = mntimes2( n_vec_mat, Mf2_mat_i);

dp_add_nvec0_i(:,:,:,2) = niT_Mf2_mat;
dp_add_nvec0_i(:,:,Ny+1:end,1) = dp_add_nvec0_i(:,:,1:end-Ny,2);
dp_add_nvec0_i(:,:,1:end-Ny,3) = dp_add_nvec0_i(:,:,Ny+1:end,2);

Qf_p_mat0_i = zeros(N_q,N_element,N_element);
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% Length of element ii [-]
    dW = dW_vec(ii);                        %% Width of element ii [-]
    
   
            
    int_StM = zeros(N_q,N_element);
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre quadrature
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] Linear interpolation of the fluid force
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 4 2 3 1]);
        dp_add_ni0_interp = sum( p_interp_vec(ones(1,3),ones(1,N_element),:,:).*dp_add_nvec0_i(:,:,ii,:), 4);
        
        i_eta_a = 1;
        for eta_a = p_vec
            
            %%[*] Pressure acts along the normal vector, so the sign is positive.
            StM = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_add_ni0_interp;
            int_StM = int_StM + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM;
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_mat0_i(:,:,ii) = int_StM;            %% External force matrix of element ii (per unit area)
end




%% [2-3] Added-mass matrix (Mf1)
%%[*] Nodal values used for the linear interpolation of the fluid force
dp_add_nvec_i = zeros(3,N_q_all,N_element,3);      
n_vec_mat = permute( n_vec_i, [ 2 3 1]);
Mf1_mat_i = permute( Mf1_mat, [3 2 1]);
niT_Mf1_mat = mntimes2( n_vec_mat, Mf1_mat_i);

dp_add_nvec_i(:,:,:,2) = niT_Mf1_mat;
dp_add_nvec_i(:,:,Ny+1:end,1) = dp_add_nvec_i(:,:,1:end-Ny,2);
dp_add_nvec_i(:,:,1:end-Ny,3) = dp_add_nvec_i(:,:,Ny+1:end,2);

Qf_p_mat_i = zeros(N_q,N_q,N_element);
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% Length of element ii [-]
    dW = dW_vec(ii);                        %% Width of element ii [-]
    
    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    
            
    int_StM = zeros(N_q);
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre quadrature
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] Linear interpolation of the fluid force
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 4 2 3 1]);
        dp_add_ni_interp = sum( p_interp_vec(ones(1,3),ones(1,N_q),:,:).*dp_add_nvec_i(:,i_vec,ii,:), 4);
        
        i_eta_a = 1;
        for eta_a = p_vec
            
            %%[*] Pressure acts along the normal vector, so the sign is positive.
            StM = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_add_ni_interp;
            int_StM = int_StM + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM;
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_mat_i(:,:,ii) = int_StM;            %% External force matrix of element ii (per unit area)
end

