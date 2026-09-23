%% Panel node point evaluation [m]


%%[0-0] Point 1 [-]
r_panel_vec_1v = Sc_mat_panel_global_1*q_vec;
r_panel_vec_1 = reshape( r_panel_vec_1v, 3, []).';

%%[0-1] Point 2 [-] (panel node points at the plate tip are interpolated from the tip nodal coordinates)
r_panel_vec_2v = Sc_mat_panel_global_2*q_vec;
r_panel_vec_2 = reshape( r_panel_vec_2v, 3, []).';

%%[0-2] Point 3 [-] (panel node points at the plate tip are interpolated from the tip nodal coordinates)
r_panel_vec_3v = Sc_mat_panel_global_3*q_vec;
r_panel_vec_3 = reshape( r_panel_vec_3v, 3, []).';

%%[0-3] Point 4 [-]
r_panel_vec_4v = Sc_mat_panel_global_4*q_vec;
r_panel_vec_4 = reshape( r_panel_vec_4v, 3, []).';



%%[0-4] Point 1 [-]
dt_r_panel_vec_1v = Sc_mat_panel_global_1*dt_q_vec;
dt_r_panel_vec_1 = reshape( dt_r_panel_vec_1v, 3, []).';

%%[0-5] Point 2 [-] (panel node points at the plate tip are interpolated from the tip nodal coordinates)
dt_r_panel_vec_2v = Sc_mat_panel_global_2*dt_q_vec;
dt_r_panel_vec_2 = reshape( dt_r_panel_vec_2v, 3, []).';

%%[0-6] Point 3 [-] (panel node points at the plate tip are interpolated from the tip nodal coordinates)
dt_r_panel_vec_3v = Sc_mat_panel_global_3*dt_q_vec;
dt_r_panel_vec_3 = reshape( dt_r_panel_vec_3v, 3, []).';

%%[0-7] Point 4 [-]
dt_r_panel_vec_4v = Sc_mat_panel_global_4*dt_q_vec;
dt_r_panel_vec_4 = reshape( dt_r_panel_vec_4v, 3, []).';