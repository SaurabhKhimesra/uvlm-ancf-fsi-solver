%% parameters
param_ver = 12.0;                                        %% parameter file version

%% Analysis settings
End_Time = 30;                                          %% Nondimensional analysis time [-]
% End_Time = 10;                                          %% Nondimensional analysis time [-]
d_t = 1.5e-3;                                           %% Nondimensional time step [-]
core_num = 8;                                           %% Number of CPUs [-]
speed_check = 0;                                        %% Speed check: 1:ON, 0:OFF [-]
alpha_v = 0.5;                                          %% 1:implicit, 0:explicit [-]
coupling_flag = 1;                                      %% 1: strong coupling (added mass computed explicitly), 0: weak coupling (staggered)


Ma = 1.0;                                              %% Mass ratio [-]
Ua = 25;                                                %% Nondimensional flow velocity [-]
% theta_a = 1e-2;                                       %% Structural damping parameter [-]: σ = D(ε + θ^* dtε), θ^* := (θ [1/s])*(U_in/L [s])
theta_v = 7.8e-4/25;                                   %% Dimensional structural damping coefficient [s] x factor [1/s] (θ^*×ω_n(7)^*/2 = 0.001, ω_n(7)^* = 2.574 [-] @ U^* = 25 ⇒ θ^* = 0.00078 @ U^* = 25)
theta_a = 0*(Ua/Ma^2)*theta_v;                           %% Structural damping parameter [-]: σ = D(ε + θ^* dtε), θ^* := (θ [1/s])*(U_in/L [s]) = (U^*/M*^2)*√(ρ_f^4 E/12ρ_m^5 H^2)θ ∝ (U^*/M*^2)*θ
C_theta_a = 0*1e-2;                                     	%% Leading-edge rotational damper coefficient [-]: Cθ^* := Cθ/(ρf*L^3*W*Uin) [-]
J_a = 0;                                                %% Nondimensional moment of inertia [-]: J^* := J/(ρf*L^4*W) [-]

dt_rz_end = 0*3e-1;                                     %% Nondimensional initial velocity at the tip (X=1 [-]) [-] (distributed linearly in X) ---------------- an initial velocity gives an impulsive response.
q_in_norm = @( time)( 0.5*sin( pi*time/0.2).*( time < 0.2 ) );              %% Disturbance body force [-] ------------------------------------------- gives a smooth response.
q_in_vec = [ 0 0 1].';                                                      %% Direction the disturbance body force acts in [-]

mode_num = 5;

%% Plot settings

i_snapshot = 50;                           %% Snapshot plot interval
Snapshot_tmin = 25;                         %% Snapshot plot start time [-]
Snapshot_tmax = End_Time;                   %% Snapshot plot end time [-]
panel_node_plot = 0;                        %% Plot the plate panel nodes, collocation points and fluid force vectors: 1:ON, 0:OFF [-]
pressure_interp_plot = 0;                   %% Interpolated fluid force vector [Pa]
% movie_format = 'mpeg';                      %% Movie format [-]
% movie_format = 'avi';
movie_format = 'wmv';                      



%% Plate parameters

mu_m = 1/Ma;                                    %% Nondimensional plate density [-]
nu = 0.3;                                       %% Poisson's ratio [-]

Length = 1.0;                                   %% Nondimensional length [-]
Width = 1.00*Length;                           	%% Nondimensional width (aspect ratio) [-]
% Width = 50.0*Length;                           	%% Nondimensional width (aspect ratio) [-]

thick = 1e-3;%(1.21/1.39e+3)*1/Ma;                    %% Nondimensional thickness [-] (a thicker plate allows a larger structural time step.)
Aa = Width*thick;                               %% Cross-sectional area [-]
Ia = Width*thick^3/12;                          %% Second moment of area [-]
eta_m = mu_m/Ua^2;                              %% Nondimensional bending stiffness [-]
zeta_m = Aa/Ia*eta_m*Length^2;                	%% Nondimensional extensional stiffness [-]
Nx = 15;                                        %% Number of elements in x [-]
Ny = 10;                                         %% Number of elements in y [-]
% Nx = 14;                                        %% Number of elements in x [-]
% Ny = 40;                                         %% Number of elements in y [-]

n_LW = 1.0;                                     %% Non-uniform mesh parameter [-]
N_gauss = 5;                                    %% Gauss-Legendre quadrature order [-]



k_gravity = 0.136^4/(1.21^3*4.989e-4);
% F_in = -9.807*k_gravity*(Ma/Ua)^2*[ 0 1 0].';                           %% Nondimensional body force [-]
F_in = 0*[ 0 1 0].';                           %% Nondimensional body force [-]


x_vec = ( (0:Nx)/Nx ).^n_LW*Length;             %% Step range of the nodal x coordinate (element coordinate system) [-]
y_vec = (0:Ny)/Ny*Width;                        %% Step range of the nodal y coordinate (element coordinate system) [-]
N_element = Nx*Ny;                             	%% Total number of elements [-]

Dp_mat = 1/(1 - nu^2)*[	1   nu  0;
                        nu  1   0;
                        0   0   (1 - nu)/2];



%% Fluid force only
flag_fluid_bench = 0;                                                       %% 1: fluid force only (the structure is treated as a rigid plate), 0: full FSI

k_omega = 1/2;
Theta_pitch = pi/6.0;                                                        %% [rad]

omega_pitch = 2*k_omega;

theta_pitch_time = @( time)( Theta_pitch*sin( omega_pitch*time) );                          %% Pitch angle [rad]
dt_theta_pitch_time = @( time)( Theta_pitch*omega_pitch.*cos( omega_pitch*time));          	%% Pitch angular velocity [rad/-]
dtt_theta_pitch_time = @( time)( -Theta_pitch*omega_pitch.^2.*sin( omega_pitch*time));     	%% Pitch angular acceleration [rad/-]

L_pitch_center =  Length/2;                                                 %% Pitch rotation centre [-]


R_pitch = @( theta)[ ...
	cos( theta)     0   sin( theta) ;
    0               1   0           ;
    -sin( theta)    0   cos( theta) ];                  %% 3D rotation matrix

dt_R_pitch = @( theta, omega) ...
    omega*[ ...
	-sin( theta)   	0   cos( theta) ;
    0               0   0           ;
    -cos( theta)    0   -sin( theta) ];                  %% 3D rotation matrix


H_1_2nd = @( k)besselh( 1, 2, k); 
H_0_2nd = @( k)besselh( 0, 2, k); 

C_theodorsen = @( k)( H_1_2nd( k)/( H_1_2nd( k) + 1i*H_0_2nd( k) ) ); 


%% Fluid parameters

U_in = 1.0;                                     %% Inflow velocity in X [-]

V_in = ones(N_element,1)*U_in*[ 1 0 0];         %% Inflow velocity vector [-]
%%[*] Remove the singularity caused by vortex points lying on top of each other
%%% W. Hoydonck et al, Validity of Viscous Core Correction Models for
%%% Self-Induced Velocity Calculations, the Journal of the American Helicopter Society, pp. 1-8, 2011.
r_eps.fine = 1e-6;                             %% Vortex core size based on the largest wake cell length [-]: (when evaluating circulation at the collocation points on the plate)
r_eps.rough = 10e-2;                          	%% Vortex core size based on the largest wake cell length [-]: (during wake time integration)
Ncore = 2;                                      %% Dimension of vortex core model
%%[*] Remove the singularity in the cross product in the denominator
%%% R. Leuthold, Multiple-Wake Vortex Lattice Method for Membrane-Wing
%%% Kites, Master of Science Thesis, p. 90, 2015.
%%% (Use a value below the vortex core radius but larger than 10^-10. Anything smaller can produce an asymmetric flow field.)
eps_v = 1e-9;         


dL_vec = diff( x_vec);
dt_wake = dL_vec(end)/U_in;                            	%% Wake shedding time step [-] (dt_wake = dL/Uin)
dt_wake_per_dt = ceil( dt_wake/d_t);                    %% Fluid solver time step / structural solver time step [-]

%% Plate boundary conditions

%%[0] Assign zero
node_r_0 = [ 1:Ny+1];                                     %% Node indices with a displacement constraint [-]
node_dxr_0 = [ 1:Ny+1];                                  	%% Node indices with an x-direction slope constraint [-]
node_dyr_0 = [ 1:Ny+1];                                   %% Node indices with a y-direction slope constraint [-]

%%[1] Constrained to a common value (flexible sheet clamped to a freely rotating rigid bar)
node_dxr_theta_c = [ ];                             %% Node indices whose x-direction slope is tied together [-]
node_dyr_theta_c = [ ];                             %% Node indices whose y-direction slope is tied together [-]

%%[2] Support region of the leading-edge rotational damper
element_C_theta = 1:Ny;                             	%% Element indices supported by the leading-edge rotational damper [-]

%% Wake analysis parameters
R_wake_x_threshold = 5.5*Length;                                    %% Threshold on the wake tip position [m]
R_wake_x_threshold_no_change = R_wake_x_threshold - 1.5*Length;     %% Threshold on the wake tip position that may still deform [m]


%% global variables
global var_param

var_param.Length = Length;
var_param.Nx = Nx;
var_param.d_t = d_t;
var_param.alpha_v = alpha_v;
var_param.theta_a = theta_a;
var_param.C_theta_a = C_theta_a;
var_param.J_a = J_a;

var_param.node_r_0 = node_r_0;
var_param.node_dxr_0 = node_dxr_0;
var_param.node_dyr_0 = node_dyr_0;
var_param.node_dxr_theta_c = node_dxr_theta_c;
var_param.node_dyr_theta_c = node_dyr_theta_c;

var_param.r_eps = r_eps;
var_param.Ncore = Ncore;
var_param.eps_v = eps_v;


