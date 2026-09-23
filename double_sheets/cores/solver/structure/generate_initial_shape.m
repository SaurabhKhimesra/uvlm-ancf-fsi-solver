%% Initial shape generation





%% Convert to element coordinate points and slope vectors, (Nx+1)*(Ny+1) of them


%%[0] Element coordinate points [ r_x r_y r_z]^T

coordinates_tmp = coordinates;
%%[0-0] Sheet 1
coordinates = [ coordinates_tmp(:,1).';                                      	%% r_x0
                coordinates_tmp(:,2).';                                        	%% r_y0
                0*coordinates_tmp(:,1).' + Height/2].';                        	%% r_z0
%%[0-1] Sheet 2
coordinates_1 = [	coordinates_tmp(:,1).';                                    	%% r_x0
                    coordinates_tmp(:,2).';                                    	%% r_y0
                    0*coordinates_tmp(:,1).' - Height/2].';                    	%% r_z0



%%[1-1] Slope vector [ dx_r_x dx_r_y dx_r_z]^T
%%[1-0] Sheet 1
dx_r0_vec = [  	ones(1,N_node);                                              	%% dX = 1
                zeros(1,N_node);                                                %% dY = 0
                zeros(1,N_node)].';                                             %% dZ = 0
dx_r0_vec = dx_r0_vec./norm_mat( dx_r0_vec);                                    %% dx_r (normalised)


%%[1-1] Sheet 2
dx_r0_vec_1 = [	ones(1,N_node);                                              	%% dX = 1
                zeros(1,N_node);                                                %% dY = 0
                zeros(1,N_node)].';                                             %% dZ = 0
dx_r0_vec_1 = dx_r0_vec_1./norm_mat( dx_r0_vec_1);                             	%% dx_r (normalised)










