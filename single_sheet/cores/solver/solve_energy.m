%% (Nondimensional) energy and work-rate evaluation
%%[*] Extract variables
q_vec = X_vec(1:N_q_all,1);
dt_q_vec = X_vec(N_q_all+1:end,1);





%% [0] Kinetic energy [-]
h_E_inertia(i_time) = 1/2*dt_q_vec.'*M_global*dt_q_vec;


%% [1] Membrane strain [-]
sum_e_Dp_e = 0;
sum_dt_e_Dp_dt_e = 0;
sum_dt_e_Dp_e = 0;
for ii = 1:N_element
    
    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    q_i_vec = q_vec(i_vec);
    dt_q_i_vec = dt_q_vec(i_vec);

    
    int_e_Dp_e = int_e_Dp_e_vec(ii);
    int_dq_e_Dp_e = Qe_eps_vec_i(:,ii);

     
 
    sum_e_Dp_e = sum_e_Dp_e + int_e_Dp_e;

    int_dt_e_Dp_e = dt_q_i_vec.'*int_dq_e_Dp_e; 
    sum_dt_e_Dp_e = sum_dt_e_Dp_e + int_dt_e_Dp_e;
    
    
    if theta_a ~= 0  
        int_dq_e_Dp_dq_e = Qd_eps_mat_i(:,:,ii);    
        
        int_dt_e_Dp_dt_e = dt_q_i_vec.'*int_dq_e_Dp_dq_e*dt_q_i_vec;     
        sum_dt_e_Dp_dt_e = sum_dt_e_Dp_dt_e + int_dt_e_Dp_dt_e;
    end
end

%%[*] Membrane strain energy
h_E_em(i_time) = 1/2*zeta_m*sum_e_Dp_e;
%%[*] Bending strain work rate (dt_E_em)
h_W_em2(i_time) = zeta_m*sum_dt_e_Dp_e; 

%%[*] Dissipation from the membrane strain rate
h_W_dm(i_time) = theta_a*zeta_m*sum_dt_e_Dp_dt_e;





%% [2] Bending strain [-]
sum_k_Dp_k = 0;
sum_dt_k_Dp_dt_k = 0;
sum_dt_k_Dp_k = 0;
for ii = 1:N_element
    
    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    q_i_vec = q_vec(i_vec);
    dt_q_i_vec = dt_q_vec(i_vec);

    
    int_k_Dp_k = int_k_Dp_k_vec(ii);
    int_dq_k_Dp_k = Qe_k_vec_i(:,ii);

     
 
    sum_k_Dp_k = sum_k_Dp_k + int_k_Dp_k;

    int_dt_k_Dp_k = dt_q_i_vec.'*int_dq_k_Dp_k; 
    sum_dt_k_Dp_k = sum_dt_k_Dp_k + int_dt_k_Dp_k;
    
    
    if theta_a ~= 0  
        int_dq_k_Dp_dq_k = Qd_k_mat_i(:,:,ii);    
        
        int_dt_k_Dp_dt_k = dt_q_i_vec.'*int_dq_k_Dp_dq_k*dt_q_i_vec;     
        sum_dt_k_Dp_dt_k = sum_dt_k_Dp_dt_k + int_dt_k_Dp_dt_k;
    end
end
%%[*] Bending strain energy
h_E_ek(i_time) = 1/2*eta_m*sum_k_Dp_k;    
%%[*] Bending strain work rate (dt_E_ek)
h_W_ek2(i_time) = eta_m*sum_dt_k_Dp_k;   

%%[*] Dissipation from the bending strain rate
h_W_dk(i_time) = theta_a*eta_m*sum_dt_k_Dp_dt_k; 



%% [3] Rotational damper [-]
sum_dt_q_Qd_theta_dt_q = 0;
for ii = 1:N_element
    
    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    q_i_vec = q_vec(i_vec);
    dt_q_i_vec = dt_q_vec(i_vec);
    
    
    
    if C_theta_a ~= 0
        
        int_Qd_theta = Qd_theta_mat_i(:,:,ii);
        int_dt_q_Qd_theta_dt_q = dt_q_i_vec.'*int_Qd_theta*dt_q_i_vec;
        sum_dt_q_Qd_theta_dt_q = sum_dt_q_Qd_theta_dt_q + int_dt_q_Qd_theta_dt_q;
    end
end

%%[*] Dissipation through the rotational damper
h_W_d_theta(i_time) = C_theta_a*sum_dt_q_Qd_theta_dt_q; 




%% [3] Rotational inertia [-]


%%[*] Rotational inertia energy
h_E_Ja(i_time) = 1/2*J_a*sum( dt_theta_0y.^2); 




%% [5] Work rate of the external force [-]

if mod( i_time, dt_wake_per_dt) == 1        

    
    time_wake_m(i_wake_time) = time;          	%% Fluid parameter update time [-]
    
    %%[*] Nodal values used for the linear interpolation of the fluid force
    dp_vec = h_dp_vec(:,i_wake_time-1);
    n_vec_i = h_n_vec(:,:,i_wake_time-1);

    dp_nvec_i = zeros(N_element,3,3);        
    dp_nvec_i(:,:,2) = ( dp_vec*ones(1,3) ).*n_vec_i;
    dp_nvec_i(Ny+1:end,:,1) = dp_nvec_i(1:end-Ny,:,2);
    dp_nvec_i(1:end-Ny,:,3) = dp_nvec_i(Ny+1:end,:,2);

    sum_dt_r_f_ext = 0;
    dt_r_f_ext_Xdist = zeros(N_element,1);
    for ii = 1:N_element

        %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
        %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
        i_vec = i_vec_v{ii};
        dt_q_i_vec = dt_q_vec(i_vec);


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
        dt_r_f_ext = dt_q_i_vec.'*int_StF;            %% External force matrix of element ii (per unit area)
        sum_dt_r_f_ext = sum_dt_r_f_ext + dt_r_f_ext;
        dt_r_f_ext_Xdist(ii) = dt_r_f_ext;
    end

    h_W_f_ext(i_wake_time-1) = sum_dt_r_f_ext;
    
    dt_r_f_ext_Xdist = mean( reshape( dt_r_f_ext_Xdist, Ny, []).', 2);	%% Average the inflow energy distribution over the Y direction.
    h_W_f_ext_Xdist(i_wake_time-1,:) = dt_r_f_ext_Xdist;

end