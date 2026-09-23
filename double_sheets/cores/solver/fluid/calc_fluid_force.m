%% Fluid force evaluation
%%%
%%% dp/ρf = (Vin + Vw + Vb - dt_r)*(τx*dx_Γ + τy*dy_Γ) + dt_Γ
%%%


%% [0] Steady fluid force
%%%
%%% [p]_lift := (Vin + Vw + Vb - dt_r)*(τx*dx_Γ + τy*dy_Γ)
%%% Γ := (A+B)^-1*[ (dt_rc - Vin - Σ_(Ny+1)^Nwake Γ_wake*Σq_wake^T)*n_i ]

%%[*] Upper surface velocity
V_surf_0 = V_surf(1:end/2,:);
V_surf1_0 = V_surf1(1:end/2,:);
%%[*] Lower surface velocity
V_surf_1 = V_surf(end/2+1:end,:);
V_surf1_1 = V_surf1(end/2+1:end,:);


%%[0-0] Static fluid force          (upper surface)
%%[0-0-0] Unit tangent vector [-]
r21_vec = r_panel_vec_2 - r_panel_vec_1;
r34_vec = r_panel_vec_3 - r_panel_vec_4;
r14_vec = r_panel_vec_1 - r_panel_vec_4;
r23_vec = r_panel_vec_2 - r_panel_vec_3;

tau_x = (r21_vec + r34_vec)/2;
tau_y = (r14_vec + r23_vec)/2;


d_x_vec = norm_mat( tau_x);
d_y_vec = norm_mat( tau_y);

tau_x = tau_x./d_x_vec;
tau_y = tau_y./d_y_vec;

%%[0-0-1] Circulation gradients: dx_Γ, dy_Γ
d_x_mat = reshape( d_x_vec(:,1), Ny, []).';
d_y_mat = reshape( d_y_vec(:,1), Ny, []).';

%%[*] For the leading term of the difference, dx_Γ1 = ( (Γ_i - Γ_{i-1})/dx )|_{i=1} = Γ_1/dx.
%%% M. Ghommem, Modeling and Analysis for Optimization of Unsteady Aeroelastic
%%% Systems, Doctoral dissertation of Virginia Polytechnic Institute and
%%% State University, p. 138, 2011. 
Gamma_mat = reshape( Gamma, Ny, []).';                                      
dx_Gamma = [ Gamma_mat(1,:);
             diff( Gamma_mat, [], 1)]./d_x_mat;

Gamma_mat2 = [ zeros(Nx,1)  Gamma_mat   zeros(Nx,1)];                       %% Central difference in y (keeps the pressure distribution symmetric)
dy_Gamma = (Gamma_mat2(:,3:end) - Gamma_mat2(:,1:end-2))./(2*d_y_mat);      %% Central difference in y
dy_Gamma(:,1) = Gamma_mat(:,1)./d_y_mat(:,1);                               %% One-sided differences at both ends
dy_Gamma(:,end) = -Gamma_mat(:,end)./d_y_mat(:,end);                        %% One-sided differences at both ends

%%[0-0-2] Nondimensional fluid force evaluation [-]
tau_x_dx_Gamma = tau_x.*( reshape( dx_Gamma.', [], 1)*ones(1,3) );
tau_y_dy_Gamma = tau_y.*( reshape( dy_Gamma.', [], 1)*ones(1,3) );


dp_add = (Gamma - old_Gamma)/d_t_wake;                                      %% Added-mass contribution
dp_lift = sum( V_surf_0.*(tau_x_dx_Gamma + tau_y_dy_Gamma), 2);               %% Steady fluid force contribution
dp_lift1 = sum( V_surf1_0.*(tau_x_dx_Gamma + tau_y_dy_Gamma), 2);             %% Steady fluid force contribution
dp_lift2 = -(tau_x_dx_Gamma + tau_y_dy_Gamma);                              %% -(τx*dxΓ + τy*dyΓ)*dt_rc 


dp_vec = dp_lift + dp_add; 
dp_mat = reshape( dp_vec, Ny, []).';




%%[0-1] Static fluid force          (lower surface)
%%[0-1-0] Unit tangent vector [-]
r21_vec_1 = r_panel_vec_2_1 - r_panel_vec_1_1;
r34_vec_1 = r_panel_vec_3_1 - r_panel_vec_4_1;
r14_vec_1 = r_panel_vec_1_1 - r_panel_vec_4_1;
r23_vec_1 = r_panel_vec_2_1 - r_panel_vec_3_1;

tau_x_1 = (r21_vec_1 + r34_vec_1)/2;
tau_y_1 = (r14_vec_1 + r23_vec_1)/2;


d_x_vec_1 = norm_mat( tau_x_1);
d_y_vec_1 = norm_mat( tau_y_1);

tau_x_1 = tau_x_1./d_x_vec_1;
tau_y_1 = tau_y_1./d_y_vec_1;

%%[0-1-1] Circulation gradients: dx_Γ, dy_Γ
d_x_mat_1 = reshape( d_x_vec_1(:,1), Ny, []).';
d_y_mat_1 = reshape( d_y_vec_1(:,1), Ny, []).';

%%[*] For the leading term of the difference, dx_Γ1 = ( (Γ_i - Γ_{i-1})/dx )|_{i=1} = Γ_1/dx.
%%% M. Ghommem, Modeling and Analysis for Optimization of Unsteady Aeroelastic
%%% Systems, Doctoral dissertation of Virginia Polytechnic Institute and
%%% State University, p. 138, 2011. 
Gamma_mat_1 = reshape( Gamma_1, Ny, []).';                                      
dx_Gamma_1 = [  Gamma_mat_1(1,:);
                diff( Gamma_mat_1, [], 1)]./d_x_mat_1;

Gamma_mat2_1 = [ zeros(Nx,1)  Gamma_mat_1   zeros(Nx,1)];                           %% Central difference in y (keeps the pressure distribution symmetric)
dy_Gamma_1 = (Gamma_mat2_1(:,3:end) - Gamma_mat2_1(:,1:end-2))./(2*d_y_mat_1);      %% Central difference in y
dy_Gamma_1(:,1) = Gamma_mat_1(:,1)./d_y_mat_1(:,1);                                 %% One-sided differences at both ends
dy_Gamma_1(:,end) = -Gamma_mat_1(:,end)./d_y_mat_1(:,end);                          %% One-sided differences at both ends

%%[0-1-2] Nondimensional fluid force evaluation [-]
tau_x_dx_Gamma_1 = tau_x_1.*( reshape( dx_Gamma_1.', [], 1)*ones(1,3) );
tau_y_dy_Gamma_1 = tau_y_1.*( reshape( dy_Gamma_1.', [], 1)*ones(1,3) );


dp_add_1 = (Gamma_1 - old_Gamma_1)/d_t_wake;                                      %% Added-mass contribution
dp_lift_1 = sum( V_surf_1.*(tau_x_dx_Gamma_1 + tau_y_dy_Gamma_1), 2);               %% Steady fluid force contribution
dp_lift1_1 = sum( V_surf1_1.*(tau_x_dx_Gamma_1 + tau_y_dy_Gamma_1), 2);             %% Steady fluid force contribution
dp_lift2_1 = -(tau_x_dx_Gamma_1 + tau_y_dy_Gamma_1);                              %% -(τx*dxΓ + τy*dyΓ)*dt_rc 


dp_vec_1 = dp_lift_1 + dp_add_1; 
dp_mat_1 = reshape( dp_vec_1, Ny, []).';


%%[*] Merge
dp_lift_all = [ dp_lift;
                dp_lift_1];
dp_lift1_all = [    dp_lift1;
                    dp_lift1_1];
dp_vec_all = [ dp_vec;
               dp_vec_1];                
                


h_dp_add(:,i_wake_time) = dp_add;                                           %% Added-mass contribution
h_dp_lift(:,i_wake_time) = dp_lift_all;                                         %% Steady fluid force contribution
h_dp_vec(:,i_wake_time) = dp_vec_all;



%% Added-mass matrix evaluation




%%[1-1] Added-mass matrix evaluation (Mf2) [-]
%%% dtΓ := (A+B)^-1*[ n_i^T*Sc_i ]dt^2_q 
%%%         + (A+B)^-1*( [ (-Σ_(i=1)^Nwake Γ_wake*Σdt_q_wake^T)*n_i] 
%%%                         + [(dt_rc - V_wake - Vin)^T*n_i]
%%%                         - dt_A*Γ)
%%%


if ~exist( 'old_dt_q_vec_all', 'var')

   old_dt_q_vec_all = dt_q_vec_all;
end

%%[1-1-0] Σ{Γwake*(dt_q_mat)^T*n}
dt_q1234_wake_mat = dt_generate_q1234_mat( rc_vec_all, r_wake_1_all, r_wake_2_all, r_wake_3_all, r_wake_4_all, ...
                                           dt_rc_vec_all, dt_r_wake_1_all, dt_r_wake_2_all, dt_r_wake_3_all, dt_r_wake_4_all);
 
Gamma_wake_dt_q1234 = dt_q1234_wake_mat.*( ones(2*N_element,1)*kron( Gamma_wake_all.', ones(1,3)) );    
Gamma_wake_dt_q1234 = [ sum( Gamma_wake_dt_q1234(:,1:3:end), 2) sum( Gamma_wake_dt_q1234(:,2:3:end), 2) sum( Gamma_wake_dt_q1234(:,3:3:end), 2)];

Gamma_wake_dt_q1234_n = sum( Gamma_wake_dt_q1234.*n_vec_i_all, 2);

                                 





%%[1-1-2] dt_A = (dt_q_mat)^T*n + q_mat^T*dt_n
dt_q1234_mat = dt_generate_q1234_mat( rc_vec_all, r_panel_vec_1_all, r_panel_vec_2_all, r_panel_vec_3_all, r_panel_vec_4_all, ...
                                      dt_rc_vec_all, dt_r_panel_vec_1_all, dt_r_panel_vec_2_all, dt_r_panel_vec_3_all, dt_r_panel_vec_4_all);
                                  
dt_q_mat_ni = inner_mat( dt_q1234_mat, n_vec_i_mat);
dt_q_mat_ni = dt_q_mat_ni(:,1:3:end);                                      	%% Copying only the 3 components inside inner_mat is unnecessary.

dt_n_vec_i_mat = repmat( dt_n_vec_i_all, [ 1 2*N_element]);
q_mat_dt_ni = inner_mat( q1234_mat, dt_n_vec_i_mat);
q_mat_dt_ni = q_mat_dt_ni(:,1:3:end);                                      	%% Copying only the 3 components inside inner_mat is unnecessary.

dt_Amat = dt_q_mat_ni + q_mat_dt_ni;
h_dt_Amat(:,:,i_wake_time) = dt_Amat;
h_Amat(:,:,i_wake_time) = A_mat;


dt_Amat1 = dt_q_mat_ni;


q_mat_Gamma_vec = q1234_mat.*( ones(2*N_element,1)*kron( Gamma_all.', ones(1,3)) );
q_mat_Gamma_vec = [ sum( q_mat_Gamma_vec(:,1:3:end), 2)  sum( q_mat_Gamma_vec(:,2:3:end), 2) sum( q_mat_Gamma_vec(:,3:3:end), 2)];

dt_Amat2_Gamma = q_mat_Gamma_vec;


%%[1-1-3] M_f2
%%[*] Fluid force vector: [p] = Mf1*dt^2_q + (Mf2_1*(dt_r - Vin - Vwake)^T*dt_ni + Mf2_2)
Mf2_vec1 = A_mat\(  -Gamma_wake_dt_q1234_n  );                                  %% dt_A = [Σ{ (dt_q_wake)^T*ni + q_wake^T*dt_ni }] is folded into the structural model.


                            
Mf2_mat = inv( A_mat);                            


%%[1-2] Added-mass matrix evaluation (Mf1) [-]

nvec_Sc_global = zeros(N_element,N_q_all);
nvec_Sc_global_1 = nvec_Sc_global;
for ii = 1:N_element

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);
    
    j_vec = 3*ii-2:3*ii;
    
    nvec_Sc_global(ii,i_vec) = nvec_Sc_global(ii,i_vec) + n_vec_i(ii,:)*Sc_mat_col_global(j_vec,i_vec);
    nvec_Sc_global_1(ii,i_vec) = nvec_Sc_global_1(ii,i_vec) + n_vec_i_1(ii,:)*Sc_mat_col_global(j_vec,i_vec);
end

nvec_Sc_global_all = blkdiag( nvec_Sc_global, nvec_Sc_global_1);

Mf1_mat = A_mat\nvec_Sc_global_all;




%% added mass effect

% h_dp_add_estimate(:,i_wake_time) = Mf1_mat*(dt_q_vec - old_dt_q_vec)/d_t_wake ...
%                                     + Mf2_vec;                                %% Added-mass contribution
h_dp_add_estimate(:,i_wake_time) = Mf1_mat*(dt_q_vec_all - old_dt_q_vec_all)/d_t_wake ...
                                    ...
                                    + Mf2_mat*( sum( (dt_rc_vec_all - V_in_all - V_wake_plate - dt_Amat2_Gamma).*dt_n_vec_i_all, 2) - dt_Amat1*Gamma_all)...
                                    + Mf2_vec1;                                 %% Added-mass contribution

                                
old_dt_q_vec_all = dt_q_vec_all;
old_Amat = A_mat;







%% Fluid force matrix assembly
%%[*] Fluid force vector: [p] = Mf1*dt^2_q + (Mf2_1*(dt_r - Vin - Vwake)^T*dt_ni + Mf2_2)



%%[2-0] Fluid force vector (Mf2_2)
%%[2-0-0] Nodal values for the linear interpolation of the fluid force (sheet 1)
dp_nvec_i = zeros(N_element,3,3);        
dp_nvec_i(:,:,2) = ( (dp_lift1 + Mf2_vec1(1:end/2))*ones(1,3) ).*n_vec_i;
dp_nvec_i(Ny+1:end,:,1) = dp_nvec_i(1:end-Ny,:,2);
dp_nvec_i(1:end-Ny,:,3) = dp_nvec_i(Ny+1:end,:,2);
%%[2-0-1] Nodal values for the linear interpolation of the fluid force (sheet 2)
dp_nvec_i_1 = zeros(N_element,3,3);        
dp_nvec_i_1(:,:,2) = ( (dp_lift1_1 + Mf2_vec1(end/2+1:end))*ones(1,3) ).*n_vec_i_1;
dp_nvec_i_1(Ny+1:end,:,1) = dp_nvec_i_1(1:end-Ny,:,2);
dp_nvec_i_1(1:end-Ny,:,3) = dp_nvec_i_1(Ny+1:end,:,2);
        

Qf_p_vec_i = zeros(N_q,N_element);
Qf_p_vec_i_1 = Qf_p_vec_i;
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% Length of element ii [-]
    dW = dW_vec(ii);                        %% Width of element ii [-]
            
    int_StF = zeros(N_q,1);
    int_StF_1 = int_StF;
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre quadrature
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] Linear interpolation of the fluid force
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 3 2 1]);
        %%[*-0] Linear interpolation of the fluid force (sheet 1)
        dp_ni_interp = sum( p_interp_vec(:,ones(1,3),:).*dp_nvec_i(ii,:,:), 3);
        %%[*-1] Linear interpolation of the fluid force (sheet 2)
        dp_ni_interp_1 = sum( p_interp_vec(:,ones(1,3),:).*dp_nvec_i_1(ii,:,:), 3);
        
        i_eta_a = 1;
        for eta_a = p_vec           
            
            
            %%[*] Pressure acts along the normal vector, so the sign is positive.
            %%[*-0] (sheet 1)
            StF = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_ni_interp.';
            int_StF = int_StF + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StF;
            %%[*-1] (sheet 2)
            StF_1 = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_ni_interp_1.';
            int_StF_1 = int_StF_1 + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StF_1;
            
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_vec_i(:,ii) = int_StF;            %% External force matrix of element ii (per unit area)
    Qf_p_vec_i_1(:,ii) = int_StF_1;            %% External force matrix of element ii (per unit area)
end



%%[2-1] Fluid force matrix ( -(τx*dxΓ + τy*dyΓ)*dt_rc )
%%[2-1-0] Nodal values for the linear interpolation of the fluid force (sheet 1)
dp_lift2_nvec_i = zeros(3,3,N_element,3);      
n_vec_mat = permute( n_vec_i, [ 2 3 1]);
dp_lift2 = permute( dp_lift2, [3 2 1]);
niT_dp_lift2 = mntimes2( n_vec_mat, dp_lift2);

dp_lift2_nvec_i(:,:,:,2) = niT_dp_lift2;
dp_lift2_nvec_i(:,:,Ny+1:end,1) = dp_lift2_nvec_i(:,:,1:end-Ny,2);
dp_lift2_nvec_i(:,:,1:end-Ny,3) = dp_lift2_nvec_i(:,:,Ny+1:end,2);

%%[2-1-1] Nodal values for the linear interpolation of the fluid force (sheet 2)
dp_lift2_nvec_i_1 = zeros(3,3,N_element,3);      
n_vec_mat_1 = permute( n_vec_i_1, [ 2 3 1]);
dp_lift2_1 = permute( dp_lift2_1, [3 2 1]);
niT_dp_lift2_1 = mntimes2( n_vec_mat_1, dp_lift2_1);

dp_lift2_nvec_i_1(:,:,:,2) = niT_dp_lift2_1;
dp_lift2_nvec_i_1(:,:,Ny+1:end,1) = dp_lift2_nvec_i_1(:,:,1:end-Ny,2);
dp_lift2_nvec_i_1(:,:,1:end-Ny,3) = dp_lift2_nvec_i_1(:,:,Ny+1:end,2);



Qf_p_lift2_mat_i = zeros(N_q,3*N_element,N_element);
Qf_p_lift2_mat_i_1 = Qf_p_lift2_mat_i;
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% Length of element ii [-]
    dW = dW_vec(ii);                        %% Width of element ii [-]
    
   
            
    int_StM = zeros(N_q,3*N_element);
    int_StM_1 = int_StM;
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre quadrature
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] Linear interpolation of the fluid force
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 4 2 3 1]);
        %%[*-0] Linear interpolation of the fluid force (sheet 1)
        dp_lift2_nvec_i_interp = sum( p_interp_vec(ones(1,3),ones(1,3),:,:).*dp_lift2_nvec_i(:,:,ii,:), 4);
        %%[*-1] Linear interpolation of the fluid force (sheet 2)
        dp_lift2_nvec_i_interp_1 = sum( p_interp_vec(ones(1,3),ones(1,3),:,:).*dp_lift2_nvec_i_1(:,:,ii,:), 4);
        
        i_eta_a = 1;
        %%[*-0] (sheet 1)
        dp_lift2_nvec_i_interp_v = zeros(3,3*N_element);
        dp_lift2_nvec_i_interp_v(:,3*ii-2:3*ii) = dp_lift2_nvec_i_interp;
        %%[*-1] (sheet 2)
        dp_lift2_nvec_i_interp_v_1 = 0*dp_lift2_nvec_i_interp_v;
        dp_lift2_nvec_i_interp_v_1(:,3*ii-2:3*ii) = dp_lift2_nvec_i_interp_1;
        for eta_a = p_vec
            
            %%[*] Pressure acts along the normal vector, so the sign is positive.
            %%[*-0] (sheet 1)
            StM = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_lift2_nvec_i_interp_v;
            int_StM = int_StM + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM;
            %%[*-1] (sheet 2)
            StM_1 = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_lift2_nvec_i_interp_v_1;
            int_StM_1 = int_StM_1 + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM_1;
            
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_lift2_mat_i(:,:,ii) = int_StM;                 %% External force matrix of element ii (per unit area)
    Qf_p_lift2_mat_i_1(:,:,ii) = int_StM_1;             %% External force matrix of element ii (per unit area)
end





%%[2-2] Fluid force matrix (Mf2_1)
%%[*] Nodal values used for the linear interpolation of the fluid force
dp_add_nvec0_i = zeros(3,2*N_element,2*N_element,3);      
n_vec_mat = permute( n_vec_i_all, [ 2 3 1]);
Mf2_mat_i = permute( Mf2_mat, [3 2 1]);
niT_Mf2_mat = mntimes2( n_vec_mat, Mf2_mat_i);

dp_add_nvec0_i(:,:,:,2) = niT_Mf2_mat;
dp_add_nvec0_i(:,:,[ Ny+1:N_element N_element+(Ny+1:N_element)],1) = dp_add_nvec0_i(:,:,[ 1:N_element-Ny N_element+(1:N_element-Ny)],2);
dp_add_nvec0_i(:,:,[ 1:N_element-Ny N_element+(1:N_element-Ny)],3) = dp_add_nvec0_i(:,:,[ Ny+1:N_element N_element+(Ny+1:N_element)],2);

Qf_p_mat0_i = zeros(N_q,2*N_element,N_element);
Qf_p_mat0_i_1 = Qf_p_mat0_i;
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% Length of element ii [-]
    dW = dW_vec(ii);                        %% Width of element ii [-]
    
   
            
    int_StM = zeros(N_q,2*N_element);
    int_StM_1 = int_StM;
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre quadrature
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] Linear interpolation of the fluid force
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 4 2 3 1]);
        %%[*-0] Linear interpolation of the fluid force (sheet 1)
        dp_add_ni0_interp = sum( p_interp_vec(ones(1,3),ones(1,2*N_element),:,:).*dp_add_nvec0_i(:,:,ii,:), 4);
        %%[*-1] Linear interpolation of the fluid force (sheet 2)
        dp_add_ni0_interp_1 = sum( p_interp_vec(ones(1,3),ones(1,2*N_element),:,:).*dp_add_nvec0_i(:,:,ii+N_element,:), 4);
        
        i_eta_a = 1;
        for eta_a = p_vec
            
            %%[*] Pressure acts along the normal vector, so the sign is positive.
            %%[*-0] (sheet 1)
            StM = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_add_ni0_interp;
            int_StM = int_StM + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM;
            %%[*-1] (sheet 2)
            StM_1 = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_add_ni0_interp_1;
            int_StM_1 = int_StM_1 + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM_1;
            
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_mat0_i(:,:,ii) = int_StM;                  %% External force matrix of element ii (per unit area)
    Qf_p_mat0_i_1(:,:,ii) = int_StM_1;             	%% External force matrix of element ii (per unit area)
end




%%[2-3] Added-mass matrix (Mf1)
%%[*] Nodal values used for the linear interpolation of the fluid force
dp_add_nvec_i = zeros(3,2*N_q_all,2*N_element,3);      
n_vec_mat = permute( n_vec_i_all, [ 2 3 1]);
Mf1_mat_i = permute( Mf1_mat, [3 2 1]);
niT_Mf1_mat = mntimes2( n_vec_mat, Mf1_mat_i);

dp_add_nvec_i(:,:,:,2) = niT_Mf1_mat;
dp_add_nvec_i(:,:,[ Ny+1:N_element N_element+(Ny+1:N_element)],1) = dp_add_nvec_i(:,:,[ 1:N_element-Ny N_element+(1:N_element-Ny)],2);
dp_add_nvec_i(:,:,[ 1:N_element-Ny N_element+(1:N_element-Ny)],3) = dp_add_nvec_i(:,:,[ Ny+1:N_element N_element+(Ny+1:N_element)],2);

Qf_p_mat_i = zeros(N_q,2*N_q,N_element);
Qf_p_mat_i_1 = Qf_p_mat_i;
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);

    dL = dL_vec(ii);                        %% Length of element ii [-]
    dW = dW_vec(ii);                        %% Width of element ii [-]
    
    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    
            
    int_StM = zeros(N_q,2*N_q);
    int_StM_1 = int_StM;
    i_xi_a = 1;
    for xi_a = p_vec                        %% Gauss-Legendre quadrature
        
        x_i = dL*(xi_a + 1)/2;
        %%[*] Linear interpolation of the fluid force
        p_interp_vec = p_interp( x_i, ii, dL_vec, Nx, Ny);
        p_interp_vec = permute( p_interp_vec, [ 4 2 3 1]);
        %%[*-0] Linear interpolation of the fluid force (sheet 1)
        dp_add_ni_interp = sum( p_interp_vec(ones(1,3),ones(1,2*N_q),:,:).*dp_add_nvec_i(:,[ i_vec i_vec+N_q_all],ii,:), 4);
        %%[*-1] Linear interpolation of the fluid force (sheet 2)
        dp_add_ni_interp_1 = sum( p_interp_vec(ones(1,3),ones(1,2*N_q),:,:).*dp_add_nvec_i(:,[ i_vec i_vec+N_q_all],ii+N_element,:), 4);
        
        i_eta_a = 1;
        for eta_a = p_vec
            
            %%[*] Pressure acts along the normal vector, so the sign is positive.
            %%[*-0] (sheet 1)
            StM = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_add_ni_interp;
            int_StM = int_StM + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM;
            %%[*-1] (sheet 2)
            StM_1 = Sc_mat_v(:,:,i_xi_a,i_eta_a,ii).'*dp_add_ni_interp_1;
            int_StM_1 = int_StM_1 + dL*dW/4*w_vec(i_xi_a)*w_vec(i_eta_a)*StM_1;
            
            i_eta_a = i_eta_a+1;
        end
        i_xi_a = i_xi_a+1;
    end
    Qf_p_mat_i(:,:,ii) = int_StM;                   %% External force matrix of element ii (per unit area)
    Qf_p_mat_i_1(:,:,ii) = int_StM_1;             	%% External force matrix of element ii (per unit area)
end






%% [3] Global matrix assembly

Qf_p_global = zeros(N_q_all,1);
Qf_p_global_1 = Qf_p_global;
for ii = 1:N_element

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    
    Qf_p_global(i_vec,1) = Qf_p_global(i_vec,1) + Qf_p_vec_i(:,ii);
    Qf_p_global_1(i_vec,1) = Qf_p_global_1(i_vec,1) + Qf_p_vec_i_1(:,ii);
end

Qf_p_mat_global = sparse(N_q_all,2*N_q_all);
Qf_p_mat_global_1 = Qf_p_mat_global;
Qf_p_mat0_global = zeros(N_q_all,2*N_element);
Qf_p_mat0_global_1 = Qf_p_mat0_global;
Qf_p_lift2_mat_global = zeros(N_q_all,3*N_element);
Qf_p_lift2_mat_global_1 = Qf_p_lift2_mat_global;
for ii = 1:N_element

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = i_vec_v{ii};
    
    Qf_p_mat0_global(i_vec,:) = Qf_p_mat0_global(i_vec,:) + squeeze( Qf_p_mat0_i(:,:,ii));
    Qf_p_mat0_global_1(i_vec,:) = Qf_p_mat0_global_1(i_vec,:) + squeeze( Qf_p_mat0_i_1(:,:,ii));
    
    Qf_p_lift2_mat_global(i_vec,:) = Qf_p_lift2_mat_global(i_vec,:) + squeeze( Qf_p_lift2_mat_i(:,:,ii));
    Qf_p_lift2_mat_global_1(i_vec,:) = Qf_p_lift2_mat_global_1(i_vec,:) + squeeze( Qf_p_lift2_mat_i_1(:,:,ii));
    
    Qf_p_mat_global(i_vec,[ i_vec i_vec+N_q_all]) = Qf_p_mat_global(i_vec,[ i_vec i_vec+N_q_all]) + squeeze( Qf_p_mat_i(:,:,ii));
    Qf_p_mat_global_1(i_vec,[ i_vec i_vec+N_q_all]) = Qf_p_mat_global_1(i_vec,[ i_vec i_vec+N_q_all]) + squeeze( Qf_p_mat_i_1(:,:,ii));
end
                                                         
