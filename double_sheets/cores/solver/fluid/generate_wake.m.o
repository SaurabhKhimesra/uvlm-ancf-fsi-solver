%% Wakeƒpƒlƒ‹ƒm[ƒh“_ã‚Ì—¬‘¬‚ðŽZo+Wake¶¬


    
if i_wake_time == 1    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ŒvŽZŠJŽn’¼Œã  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%[2-0] Œã‰ƒpƒlƒ‹ƒm[ƒhÀ•WŽæ“¾ [-]
    i_trail = N_element-Ny+1:N_element;
    N_trail = length( i_trail);                                             %% Œã‰ƒpƒlƒ‹” [-]
    r_panel_vec_2_end = r_panel_vec_2(i_trail,:);
    r_panel_vec_3_end = r_panel_vec_3(i_trail,:);
    r_panel_vec_31_end = r_panel_vec_3_end(1,:);                             %% Y=0 [-]‚É‚¨‚¯‚éWakeƒm[ƒh


    %%[2-1] •½”Âã‚Ì‰Qƒpƒlƒ‹‚É‚æ‚éŒã‰ƒpƒlƒ‹ƒm[ƒh‚Ö‚Ì—U“±‘¬“xŒvŽZ [-]
    r_trail_vec = [ r_panel_vec_2_end;
                    r_panel_vec_31_end];       
    V_gamma_wake = V_wake_func( r_trail_vec, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    

    V_wake_2 = V_gamma_wake(1:end-1,:) + V_in(i_trail,:) - dt_r_panel_vec_2(i_trail,:);
    V_wake_31 = V_gamma_wake(end,:) + V_in(i_trail(1),:) - dt_r_panel_vec_3(i_trail(1),:);


    %%[2-2-0] Wakeƒpƒlƒ‹ƒm[ƒhÀ•W¶¬ [-]    
    r_wake_1 = r_panel_vec_2_end;
    r_wake_4 = r_panel_vec_3_end;
    r_wake_2 = r_panel_vec_2_end + V_wake_2*d_t_wake;
    r_wake_31 = r_panel_vec_31_end + V_wake_31*d_t_wake;
    r_wake_3 = [ r_wake_31;
                 r_wake_2(1:end-1,:)];

    h_r_wake(i_wake_time) = { [ r_wake_1;
                                r_wake_2;
                                r_wake_3;
                                r_wake_4] };
                            
    %%[2-2-1] Wakeƒm[ƒhˆÚ—¬‘¬“x [-]
    dt_r_wake_1 = dt_r_panel_vec_2(i_trail,:);
    dt_r_wake_4 = dt_r_panel_vec_3(i_trail,:);
    dt_r_wake_2 = V_wake_2;                            
    dt_r_wake_31 = V_wake_31;
    dt_r_wake_3 = [ dt_r_wake_31;                                           %% Y=0 [m]‚É‚¨‚¯‚éWakeƒm[ƒh‡B‚ÌˆÚ—¬‘¬“x
                    dt_r_wake_2(1:end-1,:)];                                %% ‚»‚Ì‘¼‚ÌWakeƒm[ƒh‡B‚ÌˆÚ—¬‘¬“x

    
    %%[2-3] zŠÂ‚Ì•Û‘¶
    Gamma_wake = Gamma_trail;
    
    h_Gamma_wake(i_wake_time) = { Gamma_wake };
    
    h_r_end(:,:,i_wake_time) = [ r_panel_vec_31_end;
                                 r_panel_vec_2_end];
    h_V_wake_end(:,:,i_wake_time) = [ V_wake_31(1,:);
                                      V_wake_2(1:Ny,:)];

elseif i_wake_time > 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2‰ñ–Ú‚ÌŒvŽZˆÈ~  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%[*] è‡’lˆÈ~‚Å‚Í•ÏŒ`‚³‚¹‚È‚¢ (Œã‰•”‚Ìwake‚ÌØ‚èŽæ‚è‚É‚æ‚é•Ï“®‚Ì–hŽ~)
    idx_R_wake_x_threshold_no_change = find( (r_wake_1(:,1) + r_wake_4(:,1))/2 > R_wake_x_threshold_no_change);    	%% R_wake_x_threshold_no_change: •ÏŒ`‚ð‹–—e‚·‚éWake––’[ˆÊ’u‚Ìè‡’l [m]
    if ~isempty( idx_R_wake_x_threshold_no_change)

        idx_R_wake_x_threshold_no_change = floor( idx_R_wake_x_threshold_no_change(1)/Ny)*Ny + 1;
        idx_R_wake_x_threshold_no_change_31 = floor( idx_R_wake_x_threshold_no_change(1)/Ny) + 1;
    end
    
    
    

    N_wake = size( r_wake_2, 1);                                            %% Wakeƒpƒlƒ‹” [-]

    %%[2-4] Wakeƒm[ƒhÀ•W
    r_panel_vec_2_end = r_panel_vec_2(i_trail,:);
    r_panel_vec_3_end = r_panel_vec_3(i_trail,:);
    r_panel_vec_31_end = r_panel_vec_3_end(1,:);
  
    %%[*] •½”ÂŒã‰‚Ìƒm[ƒh“_‚àŠÜ‚ß‚éD [-]
    old_r_wake_2 = [ r_panel_vec_2_end;             
                     r_wake_2];
    old_r_wake_31 = [ r_panel_vec_31_end;
                      r_wake_3(1:Ny:end,:)];                                %% Y=0 [m]‚É‚¨‚¯‚éWakeƒm[ƒh
    r_wake_23 = [ old_r_wake_2;
                  old_r_wake_31];
              
              
    %%[2-5] •½”ÂŒã‰‚ÌˆÚ“®‘¬“x [-]    
    ii_trail = 1:N_trail;
    N_wake_trail = N_wake+N_trail; 
    dt_r_panel_vec_2_wake = zeros(N_wake_trail,3);
    dt_r_panel_vec_31_wake = zeros(N_wake_trail/Ny,3);
    dt_r_panel_vec_2_wake(ii_trail,:) = dt_r_panel_vec_2(i_trail,:);
    dt_r_panel_vec_31_wake(1,:) = dt_r_panel_vec_3(i_trail(1),:);
    
    %%[2-6] Wakeƒpƒlƒ‹ƒm[ƒhÀ•W¶¬ [-]   
    %%[*] —z“IEuler–@‚Å‚Í•sˆÀ’è
    %%% H. ABEDI, Development of Vortex Filament Method for Aerodynamic
    %%% Loads on Rotor Blades, Doctor thesis of Chalmers University of Technology, p. 33, 2013.  
    
        
    
    %%[*] yk1z
    %%% RK4 method    
    %%[*] •½”Âã‚Ì‰Qƒpƒlƒ‹‚É‚æ‚éWakeƒpƒlƒ‹ƒm[ƒh‚Ö‚Ì—U“±‘¬“xŒvŽZ [-]           
    V_gamma_wake_k1 = V_wake_func( r_wake_23, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] Wakeƒpƒlƒ‹‚É‚æ‚éWakeƒpƒlƒ‹ƒm[ƒh‚Ö‚Ì—U“±‘¬“xŒvŽZ [-]     
    V_gamma_w_wake_k1 = V_wake_func( r_wake_23, r_wake_1, r_wake_2, r_wake_3, r_wake_4, Gamma_wake, var_param, 0);
         
    %%[*] •½”ÂŒã‰ + Wakeƒpƒlƒ‹ƒm[ƒh‚Å‚Ì—¬‘¬ [-]
    V_in_wake = ones(N_wake_trail,1)*U_in*[ 1 0 0]; 
    V_wake_2_k1 = V_gamma_wake_k1(1:N_wake_trail,:) + V_gamma_w_wake_k1(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k1 = V_gamma_wake_k1(N_wake_trail+1:end,:) + V_gamma_w_wake_k1(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;

        
        
    %%[*] yk2z
    r_wake_2_k2 = old_r_wake_2 + V_wake_2_k1*d_t_wake/2;
    r_wake_31_k2 = old_r_wake_31 + V_wake_31_k1*d_t_wake/2;   
    r_wake_3_k2 = zeros(N_wake_trail,3);
    r_wake_3_k2(1:Ny:end,:) = r_wake_31_k2;                                   %% Y=0 [-]‚É‚¨‚¯‚éWakeƒm[ƒh‡B
    idx_r2 = 1:N_wake_trail;
    idx_r2(1:Ny:end) = [];
    r_wake_3_k2(idx_r2,:) = r_wake_2_k2(idx_r2-1,:);                          %% ‚»‚Ì‘¼‚ÌWakeƒm[ƒh‡B
    r_wake_1_k2 = [ r_panel_vec_2_end;
                   r_wake_2_k2(1:end-Ny,:)];
    r_wake_4_k2 = [ r_panel_vec_3_end;
                   r_wake_3_k2(1:end-Ny,:)];  
               
    %%[*] zŠÂ‚Ì•Û‘¶
    %%% Kutta‚ÌðŒ‚ð–ž‚½‚·‚æ‚¤‚É•úo‰Q‚ðŒˆ‚ß‚éD
    Gamma_wake = [ Gamma_trail;
                   Gamma_wake];
             

    %%[*] •½”Âã‚Ì‰Qƒpƒlƒ‹‚É‚æ‚éWakeƒpƒlƒ‹ƒm[ƒh‚Ö‚Ì—U“±‘¬“xŒvŽZ [-]                       
    r_wake_23_k2 = [ r_wake_2_k2;
                    r_wake_31_k2];
    V_gamma_wake_k2 = V_wake_func( r_wake_23_k2, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] Wakeƒpƒlƒ‹‚É‚æ‚éWakeƒpƒlƒ‹ƒm[ƒh‚Ö‚Ì—U“±‘¬“xŒvŽZ [-]     
    V_gamma_w_wake_k2 = V_wake_func( r_wake_23_k2, r_wake_1_k2, r_wake_2_k2, r_wake_3_k2, r_wake_4_k2, Gamma_wake, var_param, 0);
                 
    %%[*] •½”ÂŒã‰ + Wakeƒpƒlƒ‹ƒm[ƒh‚Å‚Ì—¬‘¬ [-]
    V_wake_2_k2 = V_gamma_wake_k2(1:N_wake_trail,:) + V_gamma_w_wake_k2(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k2 = V_gamma_wake_k2(N_wake_trail+1:end,:) + V_gamma_w_wake_k2(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;
         
    
    
    
    %%[*] yk3z    
    %%[*] •½”ÂŒã‰ + Wakeƒpƒlƒ‹ƒm[ƒh [-]
    r_wake_2_k3 = old_r_wake_2 + V_wake_2_k2*d_t_wake/2;
    r_wake_31_k3 = old_r_wake_31 + V_wake_31_k2*d_t_wake/2;
    
    r_wake_3_k3 = zeros(N_wake_trail,3);
    r_wake_3_k3(1:Ny:end,:) = r_wake_31_k3;                                       %% Y=0 [-]‚É‚¨‚¯‚éWakeƒm[ƒh‡B
    r_wake_3_k3(idx_r2,:) = r_wake_2_k3(idx_r2-1,:);                              %% ‚»‚Ì‘¼‚ÌWakeƒm[ƒh‡B
    r_wake_1_k3 = [ r_panel_vec_2_end;
                    r_wake_2_k3(1:end-Ny,:)];
    r_wake_4_k3 = [ r_panel_vec_3_end;
                    r_wake_3_k3(1:end-Ny,:)];        
                

    r_wake_23_k3 = [ r_wake_2_k3;
                    r_wake_31_k3];
    V_gamma_wake_k3 = V_wake_func( r_wake_23_k3, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] Wakeƒpƒlƒ‹‚É‚æ‚éWakeƒpƒlƒ‹ƒm[ƒh‚Ö‚Ì—U“±‘¬“xŒvŽZ [-]     
    V_gamma_w_wake_k3 = V_wake_func( r_wake_23_k3, r_wake_1_k3, r_wake_2_k3, r_wake_3_k3, r_wake_4_k3, Gamma_wake, var_param, 0);
                
    %%[*] •½”ÂŒã‰ + Wakeƒpƒlƒ‹ƒm[ƒh‚Å‚Ì—¬‘¬ [-]
    V_wake_2_k3 = V_gamma_wake_k3(1:N_wake_trail,:) + V_gamma_w_wake_k3(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k3 = V_gamma_wake_k3(N_wake_trail+1:end,:) + V_gamma_w_wake_k3(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;
         
    
    
    
    %%[*] yk4z   
    %%[*] •½”ÂŒã‰ + Wakeƒpƒlƒ‹ƒm[ƒh [-]
    r_wake_2_k4 = old_r_wake_2 + V_wake_2_k3*d_t_wake;
    r_wake_31_k4 = old_r_wake_31 + V_wake_31_k3*d_t_wake;
    
    r_wake_3_k4 = zeros(N_wake_trail,3);
    r_wake_3_k4(1:Ny:end,:) = r_wake_31_k4;                                       %% Y=0 [-]‚É‚¨‚¯‚éWakeƒm[ƒh‡B
    r_wake_3_k4(idx_r2,:) = r_wake_2_k4(idx_r2-1,:);                              %% ‚»‚Ì‘¼‚ÌWakeƒm[ƒh‡B
    r_wake_1_k4 = [ r_panel_vec_2_end;
                    r_wake_2_k4(1:end-Ny,:)];
    r_wake_4_k4 = [ r_panel_vec_3_end;
                    r_wake_3_k4(1:end-Ny,:)];        
                

    r_wake_23_k4 = [ r_wake_2_k4;
                    r_wake_31_k4];
    V_gamma_wake_k4 = V_wake_func( r_wake_23_k4, r_panel_vec_1, r_panel_vec_2, r_panel_vec_3, r_panel_vec_4, Gamma, var_param, 0);
    
    %%[*] Wakeƒpƒlƒ‹‚É‚æ‚éWakeƒpƒlƒ‹ƒm[ƒh‚Ö‚Ì—U“±‘¬“xŒvŽZ [-]     
    V_gamma_w_wake_k4 = V_wake_func( r_wake_23_k4, r_wake_1_k4, r_wake_2_k4, r_wake_3_k4, r_wake_4_k4, Gamma_wake, var_param, 0);
                 
    %%[*] •½”ÂŒã‰ + Wakeƒpƒlƒ‹ƒm[ƒh‚Å‚Ì—¬‘¬ [-]
    V_wake_2_k4 = V_gamma_wake_k4(1:N_wake_trail,:) + V_gamma_w_wake_k4(1:N_wake_trail,:) + V_in_wake - dt_r_panel_vec_2_wake;
    V_wake_31_k4 = V_gamma_wake_k4(N_wake_trail+1:end,:) + V_gamma_w_wake_k4(N_wake_trail+1:end,:) + V_in_wake(1:N_wake_trail/Ny,:) - dt_r_panel_vec_31_wake;
    
    
     
    
    %%[*] ŽžŠÔ”­“W (RK4)
    V_wake_2 = (V_wake_2_k1 + 2*V_wake_2_k2 + 2*V_wake_2_k3 + V_wake_2_k4)/6;
    V_wake_31 = (V_wake_31_k1 + 2*V_wake_31_k2 + 2*V_wake_31_k3 + V_wake_31_k4)/6;
    
    
    %%[*] è‡’lˆÈ~‚Å‚Í•ÏŒ`‚³‚¹‚È‚¢ (Œã‰•”‚Ìwake‚ÌØ‚èŽæ‚è‚É‚æ‚é•Ï“®‚Ì–hŽ~)
    if ~isempty( idx_R_wake_x_threshold_no_change)
        
        idx_2 = idx_R_wake_x_threshold_no_change:N_wake_trail;
        idx_31 = idx_R_wake_x_threshold_no_change_31:size( old_r_wake_31, 1);
        
        V_wake_2(idx_2,:) = V_in_wake(idx_2,:);
        V_wake_31(idx_31,:) = V_in_wake(idx_31,:);
    end   
    
    
    

    r_wake_2 = old_r_wake_2 + V_wake_2*d_t_wake;
    r_wake_31 = old_r_wake_31 + V_wake_31*d_t_wake;
    
    r_wake_3 = zeros(N_wake_trail,3);
    r_wake_3(1:Ny:end,:) = r_wake_31;                                       %% Y=0 [m]‚É‚¨‚¯‚éWakeƒm[ƒh‡B
    r_wake_3(idx_r2,:) = r_wake_2(idx_r2-1,:);                              %% ‚»‚Ì‘¼‚ÌWakeƒm[ƒh‡B
    r_wake_1 = [ r_panel_vec_2_end;
             	 r_wake_2(1:end-Ny,:)];
    r_wake_4 = [ r_panel_vec_3_end;
           	     r_wake_3(1:end-Ny,:)];        
                
                
                
                

    h_r_wake(i_wake_time) = { [ r_wake_1;
                                r_wake_2;
                                r_wake_3;
                                r_wake_4] }; 
                            
    %%[*] Wakeƒm[ƒhˆÚ—¬‘¬“x [-]
    dt_r_wake_2 = V_wake_2;                            
    dt_r_wake_31 = V_wake_31;
    
    dt_r_wake_3 = zeros(N_wake_trail,3);
    dt_r_wake_3(1:Ny:end,:) = dt_r_wake_31;                                	%% Y=0 [m]‚É‚¨‚¯‚éWakeƒm[ƒh‡B‚ÌˆÚ—¬‘¬“x
    dt_r_wake_3(idx_r2,:) = dt_r_wake_2(idx_r2-1,:);                      	%% ‚»‚Ì‘¼‚ÌWakeƒm[ƒh‡B‚ÌˆÚ—¬‘¬“x
    dt_r_wake_1 = [ dt_r_panel_vec_2(i_trail,:);
                    dt_r_wake_2(1:end-Ny,:)];
    dt_r_wake_4 = [ dt_r_panel_vec_3(i_trail,:);
                    dt_r_wake_3(1:end-Ny,:)];                           
                            
                            
    h_Gamma_wake(i_wake_time) = { Gamma_wake };

    h_r_end(:,:,i_wake_time) = [ r_panel_vec_31_end;
                                 r_panel_vec_2_end];
    h_V_wake_end(:,:,i_wake_time) = [ dt_r_wake_31(1,:);
                                      dt_r_wake_2(1:Ny,:)];

end




%% Wake––’[íœ(ŒvŽZƒRƒXƒgíŒ¸)





idx_R_wake_x_threshold = find( (r_wake_1(:,1) + r_wake_4(:,1))/2 > R_wake_x_threshold);     %% R_wake_x_threshold: Wake––’[ˆÊ’u‚Ìè‡’l [m]
if ~isempty( idx_R_wake_x_threshold)
    idx_R_wake_x_threshold = floor( idx_R_wake_x_threshold(1)/Ny)*Ny + 1;

    r_wake_1(idx_R_wake_x_threshold:end,:) = [];
    r_wake_2(idx_R_wake_x_threshold:end,:) = [];
    r_wake_3(idx_R_wake_x_threshold:end,:) = [];
    r_wake_4(idx_R_wake_x_threshold:end,:) = [];
    
    dt_r_wake_1(idx_R_wake_x_threshold:end,:) = [];
    dt_r_wake_2(idx_R_wake_x_threshold:end,:) = [];
    dt_r_wake_3(idx_R_wake_x_threshold:end,:) = [];
    dt_r_wake_4(idx_R_wake_x_threshold:end,:) = [];

    Gamma_wake(idx_R_wake_x_threshold:end) = [];
end




