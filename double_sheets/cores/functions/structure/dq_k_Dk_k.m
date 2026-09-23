function [ out, out1] = dq_k_Dk_k( q_i_vec, Dk_mat, dx_n_Sc_struct, dq_vec, p_vec, theta_a)


N_q = length( q_i_vec);                                             %% Number of nodal components per element
lgth_p = length( p_vec);


%% Slope evaluation

%%[*] Zero components excluded (Sc*q = [ S1*I S2*I ... S12*I]*q = S1*q1_r + S2*q1_dx_r + ... + S12*q4_dy_r, I∈R^3*3, qi_r,qi_dx_r,qi_dy_r∈R^3)
q_i_mat = q_i_vec(:,ones(1,N_q),ones(1,lgth_p),ones(1,lgth_p));
q_i_mat = permute( q_i_mat, [ 1 5 3 4 2]);                       	%% Move the component index used when differencing ∂qκ into dimension 5.
q_i_mat = reshape( q_i_mat, 3, [], lgth_p, lgth_p, N_q);            %% Match the shape-function row count (for the product). Dim 1: coordinate component (x,y,z), dim 2: direction (r,dx_r,dy_r), dims 3-4: Gauss quadrature points, dim 5: component index used when differencing ∂qκ



q_i_vec = q_i_vec(:,1,ones(1,lgth_p),ones(1,lgth_p));
q_i_vec = reshape( q_i_vec, 3, [], lgth_p, lgth_p);                 %% Match the shape-function row count (for the product). Dim 1: coordinate component (x,y,z), dim 2: direction (r,dx_r,dy_r), dims 3-4: Gauss quadrature points


n_vec_norm_n1 = n_vec_norm_n_f( q_i_mat + dq_vec, dx_n_Sc_struct);
n_vec_norm_n = n_vec_norm_n_f( q_i_vec, dx_n_Sc_struct);

k_v1_mat = k_v( q_i_mat + dq_vec, n_vec_norm_n1, dx_n_Sc_struct);
k_v_vec = k_v( q_i_vec, n_vec_norm_n, dx_n_Sc_struct);

dq_k_v = ( k_v1_mat - k_v_vec(:,ones(1,N_q),:,:) )/dq_vec(1,1,1,1,1);


dq_k_v_Dk_mat = mntimes2( permute( dq_k_v, [ 2 1 3 4]), Dk_mat);
out = mntimes2( dq_k_v_Dk_mat, k_v_vec);  

if theta_a ~= 0                                                     %% Evaluated only when structural damping is present (saves computation)
    out1 = mntimes2( dq_k_v_Dk_mat, dq_k_v);  
else
    out1 = 0;
end

end


%% Normal vector (n = dx_r x dy_r)
function out = n_vec_f( q_vec, dx_n_Sc_struct)

dx_Sc_mat_v = dx_n_Sc_struct.dx_Sc_mat_v;
dy_Sc_mat_v = dx_n_Sc_struct.dy_Sc_mat_v;
dx_Sc_mat_v2 = dx_n_Sc_struct.dx_Sc_mat_v2;
dy_Sc_mat_v2 = dx_n_Sc_struct.dy_Sc_mat_v2;


if size( q_vec, 5) == 1
    
    dx_r = mntimes2_fast( dx_Sc_mat_v, q_vec);
    dy_r = mntimes2_fast( dy_Sc_mat_v, q_vec);
else

    dx_r = mntimes2_fast( dx_Sc_mat_v2, q_vec);
    dy_r = mntimes2_fast( dy_Sc_mat_v2, q_vec);
end
    
out =  cross_fast( dx_r, dy_r);

end

%% n/||n||^1
function out = n_vec_norm_n_f( q_vec, dx_n_Sc_struct)


n_vec = n_vec_f( q_vec, dx_n_Sc_struct);

%%[*] Is n/||n||^3 wrong?
%%% (Hui Wan et al., Study of Strain Energy in Deformed Insect Wings, Dynamic Behavior of Materials, 
%%%  Proceedings of the 2011 AnnualConference on Experimental and Applied
%%%  Mechanics, Vol. 1, 6 pages.)
norm_n_3 = norm2( n_vec).^1;

out = n_vec./norm_n_3(ones(3,1),:,:,:);

end


%% norm
function out = norm2( a)

out = sqrt( sum( a.^2, 1));

end

%% Curvatures κxx, κyy, κxy
function out = k_v( q_vec, n_vec_norm_n, dx_n_Sc_struct)


dx2_Sc_mat_v = dx_n_Sc_struct.dx2_Sc_mat_v;
dy2_Sc_mat_v = dx_n_Sc_struct.dy2_Sc_mat_v;
dxy_Sc_mat_v = dx_n_Sc_struct.dxy_Sc_mat_v;

dx2_Sc_mat_v2 = dx_n_Sc_struct.dx2_Sc_mat_v2;
dy2_Sc_mat_v2 = dx_n_Sc_struct.dy2_Sc_mat_v2;
dxy_Sc_mat_v2 = dx_n_Sc_struct.dxy_Sc_mat_v2;
              


if size( q_vec, 5) == 1
    out = [	sum( mntimes2_fast( dx2_Sc_mat_v, q_vec).*n_vec_norm_n, 1);
            sum( mntimes2_fast( dy2_Sc_mat_v, q_vec).*n_vec_norm_n, 1);
            2*sum( mntimes2_fast( dxy_Sc_mat_v, q_vec).*n_vec_norm_n, 1)];
else
    out = [	sum( mntimes2_fast( dx2_Sc_mat_v2, q_vec).*n_vec_norm_n, 1);
            sum( mntimes2_fast( dy2_Sc_mat_v2, q_vec).*n_vec_norm_n, 1);
            2*sum( mntimes2_fast( dxy_Sc_mat_v2, q_vec).*n_vec_norm_n, 1)];
end

end


%% Cross product
function out = cross_fast( a, b)


out = a([2 3 1],:,:,:).*b([3 1 2],:,:,:) - a([3 1 2],:,:,:).*b([2 3 1],:,:,:);

end






