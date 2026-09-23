function [ X_center_disp, Z_center_disp] = r_center_disp( data)

N_element = data.N_element;
Nx = data.Nx;
Ny = data.Ny;
N_qi = data.N_qi;
N_q = data.N_q;
nodes = data.nodes;
h_X_vec = data.h_X_vec;


%% [1] Extract the nodal displacement data (mid-span)

if mode( Ny, 2)
    %%[*] Odd number of spanwise elements (average the position at each node)
    
    %% r_X
    %%[*] 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %%[*] 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    idx_center_element = ceil( Ny/2):Ny:N_element-ceil( Ny/2);
    node_num = reshape( nodes(idx_center_element,:).', [], 1);

    i_vec = repmat( ( N_qi*(node_num - 1)+1 ), [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( node_num) 1]);
    i_vec = reshape( i_vec.',1,[]);
    %%%
    %%% ^  Y
    %%% | q_(4) ---- q_(3)
    %%% |  |           |
    %%% | q_(1) ---- q_(2)
    %%% O -------------------> X 
    %%%
    i_vec_1elem_X = i_vec(1:N_q:end);                       %% X coordinate of node 1 within one element
    i_vec_2elem_X = i_vec(N_qi+1:N_q:end);              	%% X coordinate of node 2 within one element
    i_vec_3elem_X = i_vec(2*N_qi+1:N_q:end);              	%% X coordinate of node 3 within one element
    i_vec_4elem_X = i_vec(3*N_qi+1:N_q:end);              	%% X coordinate of node 4 within one element
    
    X_center_disp1 = h_X_vec([ i_vec_1elem_X i_vec_2elem_X(end)],:);
    X_center_disp2 = h_X_vec([ i_vec_4elem_X i_vec_3elem_X(end)],:);

    X_center_disp =(X_center_disp1 + X_center_disp2)/2;
    
    %% r_Z
    %%[*] 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %%[*] 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec_1elem_Z = i_vec(3:N_q:end);                       %% Z coordinate of node 1 within one element
    i_vec_2elem_Z = i_vec(N_qi+3:N_q:end);              	%% Z coordinate of node 2 within one element
    i_vec_3elem_Z = i_vec(2*N_qi+3:N_q:end);              	%% Z coordinate of node 3 within one element
    i_vec_4elem_Z = i_vec(3*N_qi+3:N_q:end);              	%% Z coordinate of node 4 within one element
    
    Z_center_disp1 = h_X_vec([ i_vec_1elem_Z i_vec_2elem_Z(end)],:);
    Z_center_disp2 = h_X_vec([ i_vec_4elem_Z i_vec_3elem_Z(end)],:);
    
    Z_center_disp = (Z_center_disp1 + Z_center_disp2)/2;
else
    %%[*] Even number of spanwise elements
    
    %% r_X
    %%[*] 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %%[*] 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    idx_center_element = Ny/2:Ny:N_element-Ny/2;
    node_num = reshape( nodes(idx_center_element,:).', [], 1);

    i_vec = repmat( ( N_qi*(node_num - 1)+1 ), [ 1 N_qi]) + repmat( 0:N_qi-1, [ length( node_num) 1]);
    i_vec = reshape( i_vec.',1,[]);
    %%%
    %%% ^  Y
    %%% | q_(4) ---- q_(3)
    %%% |  |           |
    %%% | q_(1) ---- q_(2)
    %%% O -------------------> X 
    %%%
    i_vec_3elem_X = i_vec(2*N_qi+1:N_q:end);              	%% X coordinate of node 3 within one element
    i_vec_4elem_X = i_vec(3*N_qi+1:N_q:end);              	%% X coordinate of node 4 within one element
    
    X_center_disp = h_X_vec([ i_vec_4elem_X i_vec_3elem_X(end)],:);


    %% r_Z
    %%[*] 9 components per node ( q_i = [ rx_i ry_i rz_i : dx_rx_i dx_ry_i dx_rz_i : dy_rx_i dy_ry_i dy_rz_i]^T ∈ R^9 )
    %%[*] 36 components per element ( q := [ q_i1^T q_i2^T q_i3^T q_i4^T]^T ∈ R^36 )
    i_vec_3elem_Z = i_vec(2*N_qi+3:N_q:end);              	%% Z coordinate of node 3 within one element
    i_vec_4elem_Z = i_vec(3*N_qi+3:N_q:end);              	%% Z coordinate of node 4 within one element
    
    Z_center_disp = h_X_vec([ i_vec_4elem_Z i_vec_3elem_Z(end)],:);
end


