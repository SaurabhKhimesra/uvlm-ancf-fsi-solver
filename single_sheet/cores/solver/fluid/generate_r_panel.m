%% パネルノード点計算 [m]


%%[0-0] 点① [-]
r_panel_vec_1v = Sc_mat_panel_global_1*q_vec;
r_panel_vec_1 = reshape( r_panel_vec_1v, 3, []).';

%%[0-1] 点② [-]　(平板末端のパネルノード点は，平板末端ノード座標より補間)
r_panel_vec_2v = Sc_mat_panel_global_2*q_vec;
r_panel_vec_2 = reshape( r_panel_vec_2v, 3, []).';

%%[0-2] 点③ [-] (平板末端のパネルノード点は，平板末端ノード座標より補間)
r_panel_vec_3v = Sc_mat_panel_global_3*q_vec;
r_panel_vec_3 = reshape( r_panel_vec_3v, 3, []).';

%%[0-3] 点④ [-]
r_panel_vec_4v = Sc_mat_panel_global_4*q_vec;
r_panel_vec_4 = reshape( r_panel_vec_4v, 3, []).';



%%[0-4] 点① [-]
dt_r_panel_vec_1v = Sc_mat_panel_global_1*dt_q_vec;
dt_r_panel_vec_1 = reshape( dt_r_panel_vec_1v, 3, []).';

%%[0-5] 点② [-]　(平板末端のパネルノード点は，平板末端ノード座標より補間)
dt_r_panel_vec_2v = Sc_mat_panel_global_2*dt_q_vec;
dt_r_panel_vec_2 = reshape( dt_r_panel_vec_2v, 3, []).';

%%[0-6] 点③ [-] (平板末端のパネルノード点は，平板末端ノード座標より補間)
dt_r_panel_vec_3v = Sc_mat_panel_global_3*dt_q_vec;
dt_r_panel_vec_3 = reshape( dt_r_panel_vec_3v, 3, []).';

%%[0-7] 点④ [-]
dt_r_panel_vec_4v = Sc_mat_panel_global_4*dt_q_vec;
dt_r_panel_vec_4 = reshape( dt_r_panel_vec_4v, 3, []).';