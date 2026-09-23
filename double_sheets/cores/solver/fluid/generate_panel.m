%% Vortex lattice method grid generation


%% [0] Panel node point coordinates [-]( r_i = Sc(x_i,y_i)*q_i, i∈{1,...,N_element} )

Sc_mat_v_panel_1 = zeros(3,N_q,N_element);
Sc_mat_v_panel_4 = Sc_mat_v_panel_1;
for ii = 1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);
    
    dL = dL_vec(ii);                        %% Length of element ii [-]
    dW = dW_vec(ii);                        %% Width of element ii [-]
            
    %%[*] Node points in element coordinates [-]
    x_i = dL*[ 1/4  1/4];   %% Points 1 and 4
    y_i = dW*[ 1    0];     %% Points 1 and 4

    Sc_mat_v_panel_1(:,:,ii) = Sc_mat( x_i(1), y_i(1), dL, dW);         
    Sc_mat_v_panel_4(:,:,ii) = Sc_mat( x_i(2), y_i(2), dL, dW);         
end

%% [1] Plate tip node coordinates [m]( r_i = Sc(x_i,y_i)*q_i, i∈{1,...,N_element} )

Sc_mat_v_panel_end_2 = zeros(3,N_q,N_element);
Sc_mat_v_panel_end_3 = Sc_mat_v_panel_end_2;
%%[*] Element numbering increases in Y, so the elements at the plate tip are N_element-Ny+1 ... N_element.
for ii = N_element-Ny+1:N_element
    
    disp( [ 'Node number:', int2str( ii), '/', int2str( N_element)]);
    
    dL = dL_vec(ii);                        %% Length of element ii [-]
    dW = dW_vec(ii);                        %% Width of element ii [-]
            
    %%[*] Node points in element coordinates [-]
    x_i = dL*[ 1  1];   %% Points 2 and 3
    y_i = dW*[ 1  0];   %% Points 2 and 3

    Sc_mat_v_panel_end_2(:,:,ii) = Sc_mat( x_i(1), y_i(1), dL, dW);         
    Sc_mat_v_panel_end_3(:,:,ii) = Sc_mat( x_i(2), y_i(2), dL, dW);         
end


%% [2] Global matrix assembly

Sc_mat_panel_global_1 = sparse(3*N_element,N_q_all);
Sc_mat_panel_global_4 = Sc_mat_panel_global_1;
for ii = 1:N_element

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);
    
    j_vec = 3*ii-2:3*ii;
    
    Sc_mat_panel_global_1(j_vec,i_vec) = Sc_mat_panel_global_1(j_vec,i_vec) + squeeze( Sc_mat_v_panel_1(:,:,ii));
    Sc_mat_panel_global_4(j_vec,i_vec) = Sc_mat_panel_global_4(j_vec,i_vec) + squeeze( Sc_mat_v_panel_4(:,:,ii));
end



%%[*] Element numbering increases in Y, so the elements at the plate tip are N_element-Ny+1 ... N_element.
Sc_mat_panel_end_global_2 = sparse(3*N_element,N_q_all);
Sc_mat_panel_end_global_3 = Sc_mat_panel_end_global_2;
for ii = N_element-Ny+1:N_element

    %% 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %% 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec = repmat( ( N_qi*(nodes(ii,:) - 1)+1 ).', [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( nodes(ii,:)) 1]);
    i_vec = reshape(i_vec.',1,[]);
    
    j_vec = 3*ii-2:3*ii;
    
    Sc_mat_panel_end_global_2(j_vec,i_vec) = Sc_mat_panel_end_global_2(j_vec,i_vec) + squeeze( Sc_mat_v_panel_end_2(:,:,ii));
    Sc_mat_panel_end_global_3(j_vec,i_vec) = Sc_mat_panel_end_global_3(j_vec,i_vec) + squeeze( Sc_mat_v_panel_end_3(:,:,ii));
end


%%[2-0] Point 2 [m] (panel node points at the plate tip are interpolated from the tip nodal coordinates)
Sc_mat_panel_global_2 = [ Sc_mat_panel_global_1(3*Ny+1:end,:);
                          4/3*( Sc_mat_panel_end_global_2(end-3*Ny+1:end,:) - Sc_mat_panel_global_1(end-3*Ny+1:end,:) ) + Sc_mat_panel_global_1(end-3*Ny+1:end,:)];

%%[2-1] Point 3 [m] (panel node points at the plate tip are interpolated from the tip nodal coordinates)
Sc_mat_panel_global_3 = [ Sc_mat_panel_global_4(3*Ny+1:end,:);
                          4/3*( Sc_mat_panel_end_global_3(end-3*Ny+1:end,:) - Sc_mat_panel_global_4(end-3*Ny+1:end,:) ) + Sc_mat_panel_global_4(end-3*Ny+1:end,:)];

                      
%% [3] Collocation point evaluation [-]
% Sc_mat_col_global = ( Sc_mat_panel_global_1 + Sc_mat_panel_global_2 + Sc_mat_panel_global_3 + Sc_mat_panel_global_4 )/4;

%%[*] Handles non-uniform element lengths
Sc_mat_col_global = sparse(3*N_element,N_q_all);
for ii = 1:N_element    
     
    j_vec = 3*ii-2:3*ii;                                %% [ rx ry rz]^T
    
    if ii <= N_element-Ny
        dL_i = dL_vec(ii);                              %% Length of element ii [-]
        dL_ip1 = dL_vec(ii+Ny);                         %% Length of the element behind element ii [-]
        Sc_mat_col_global(j_vec,:) = Sc_mat_col_global(j_vec,:) ...
                                            + diag( [ (dL_i + dL_ip1)/(3*dL_i + dL_ip1)     1/2     (dL_i + dL_ip1)/(3*dL_i + dL_ip1)])*( Sc_mat_panel_global_1(j_vec,:) + Sc_mat_panel_global_4(j_vec,:) )/2 ...
                                           	+ diag( [ 2*dL_i/(3*dL_i + dL_ip1)              1/2     2*dL_i/(3*dL_i + dL_ip1)])*( Sc_mat_panel_global_2(j_vec,:) + Sc_mat_panel_global_3(j_vec,:) )/2;
    else
        Sc_mat_col_global(j_vec,:) = Sc_mat_col_global(j_vec,:) ...
                                            + ( Sc_mat_panel_global_1(j_vec,:) + Sc_mat_panel_global_2(j_vec,:) + Sc_mat_panel_global_3(j_vec,:) + Sc_mat_panel_global_4(j_vec,:) )/4;
    end
end


%% [4] Normal vector evaluation [-]

Sc_mat_31 = Sc_mat_panel_global_3 - Sc_mat_panel_global_1;
Sc_mat_24 = Sc_mat_panel_global_2 - Sc_mat_panel_global_4;



