function New_C=fun_update3D(C,u_vel_t,v_vel_t,w_vel_t,Dt,Dspace,Source)
% fun_rule returns a New value field based on Cellular automata local rule
% Input C is the value field to be calculated
% Input parameter is an object contains all neede parameters.

%% Parameter initilization
[IX,IY,IZ]=size(u_vel_t);
New_C=zeros(IX,IY,IZ);

ii=2:IX-1;
jj=2:IY-1;
kk=2:IZ-1;


U_e_tp=-u_vel_t;
U_n_tp=-v_vel_t;
U_w_tp=u_vel_t;
U_s_tp=v_vel_t;
U_u_tp=w_vel_t;
U_d_tp=-w_vel_t;


U_ne_tp=-u_vel_t*cos(4/pi)-v_vel_t*cos(4/pi);
U_nw_tp=u_vel_t*cos(4/pi)-v_vel_t*cos(4/pi);
U_sw_tp=u_vel_t*cos(4/pi)+v_vel_t*cos(4/pi);
U_se_tp=-u_vel_t*cos(4/pi)+v_vel_t*cos(4/pi);

U_e=1/2*(abs(U_e_tp)+U_e_tp);
U_n=1/2*(abs(U_n_tp)+U_n_tp);
U_w=1/2*(abs(U_w_tp)+U_w_tp);
U_s=1/2*(abs(U_s_tp)+U_s_tp);
U_u=1/2*(abs(U_u_tp)+U_u_tp);
U_d=1/2*(abs(U_d_tp)+U_d_tp);

U_ne=1/2*(abs(U_ne_tp)+U_ne_tp);
U_nw=1/2*(abs(U_nw_tp)+U_nw_tp);
U_sw=1/2*(abs(U_sw_tp)+U_sw_tp);
U_se=1/2*(abs(U_se_tp)+U_se_tp);

clear U_e_tp U_n_tp U_w_tp U_s_tp U_ne_tp U_nw_tp U_se_tp U_sw_tp


wa=0.8;
wb=0.5;

diffusion_flag=1;   % 0 -> no dispersion

% dispersion coefficient

Kx=1e-3;  % horizental dispersion coefficient on x axis
Ky=1e-3;
Kxy=1e-3;
Kyx=1e-3;
Kz=1e-3;

% deposition and reaction
% delta=zeros(IX,IY);   % deposition
% lambda=zeros(IX,IY);  % reaction

%% Local configuration
%
%   (NW)       (N)       (NE)
%   i-1,j+1   i,j+1     i+1,j+1
%
%    (W)       (O)       (E)
%   i-1,j      i,j      i+1,j
%
%   (SW)       (S)       (SE)
%   i-1,j-1   i,j-1     i+1,j-1

%% Function body
if diffusion_flag == 0
% Local update in main body
New_C(ii,jj,kk)=C(ii,jj,kk)+...
    wa*Dt/Dspace*(...
    U_n(ii,jj,kk).*(C(ii,jj+1,kk)-C(ii,jj,kk))+...
    U_s(ii,jj,kk).*(C(ii,jj-1,kk)-C(ii,jj,kk))+...
    U_e(ii,jj,kk).*(C(ii+1,jj,kk)-C(ii,jj,kk))+...
    U_w(ii,jj,kk).*(C(ii-1,jj,kk)-C(ii,jj,kk))+...
    U_u(ii,jj,kk).*(C(ii,jj,kk+1)-C(ii,jj,kk))+...
    U_d(ii,jj,kk).*(C(ii,jj,kk-1)-C(ii,jj,kk)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_ne(ii,jj,kk).*(C(ii+1,jj+1,kk)-C(ii,jj,kk))+...
    U_nw(ii,jj,kk).*(C(ii-1,jj+1,kk)-C(ii,jj,kk))+...
    U_sw(ii,jj,kk).*(C(ii-1,jj-1,kk)-C(ii,jj,kk))+...
    U_se(ii,jj,kk).*(C(ii+1,jj-1,kk)-C(ii,jj,kk))  )+...
    Source(ii,jj,kk)*Dt;

%%% Boundary Conditions
% South Boundary  (Y=1)
New_C(ii,1,kk)=C(ii,1,kk)+...
    wa*Dt/Dspace*(...
    U_n(ii,1,kk).*(C(ii,2,kk)-C(ii,1,kk))+...
    U_e(ii,1,kk).*(C(ii+1,1,kk)-C(ii,1,kk))+...
    U_w(ii,1,kk).*(C(ii-1,1,kk)-C(ii,1,kk))+...
    U_u(ii,1,kk).*(C(ii,1,kk+1)-C(ii,1,kk))+...
    U_d(ii,1,kk).*(C(ii,1,kk-1)-C(ii,1,kk)) )+...
    wb*Dt/Dspace/sqrt(2)*(...
    U_ne(ii,1,kk).*(C(ii+1,2,kk)-C(ii,1,kk))+...
    U_nw(ii,1,kk).*(C(ii-1,2,kk)-C(ii,1,kk))  )+...
    Source(ii,1,kk)*Dt;

% North Boundary  (Y=IY)
New_C(ii,IY,kk)=C(ii,IY,kk)+...
    wa*Dt/Dspace*(...
    U_s(ii,IY,kk).*(C(ii,IY-1,kk)-C(ii,IY,kk))+...
    U_e(ii,IY,kk).*(C(ii+1,IY,kk)-C(ii,IY,kk))+...
    U_w(ii,IY,kk).*(C(ii-1,IY,kk)-C(ii,IY,kk))+...
    U_u(ii,IY,kk).*(C(ii,IY,kk+1)-C(ii,IY,kk))+...
    U_d(ii,IY,kk).*(C(ii,IY,kk-1)-C(ii,IY,kk)) )+...
    wb*Dt/Dspace/sqrt(2)*(...
    U_se(ii,IY,kk).*(C(ii+1,IY-1,kk)-C(ii,IY,kk))+...
    U_sw(ii,IY,kk).*(C(ii-1,IY-1,kk)-C(ii,IY,kk))  )+...
    Source(ii,IY,kk)*Dt;

% East Boundary  (X=IX)
New_C(IX,jj,kk)=C(IX,jj,kk)+...*
    wa*Dt/Dspace.*(...
    U_s(IX,jj,kk).*(C(IX,jj-1,kk)-C(IX,jj,kk))+...
    U_w(IX,jj,kk).*(C(IX-1,jj,kk)-C(IX,jj,kk))+...
    U_n(IX,jj,kk).*(C(IX,jj+1,kk)-C(IX,jj,kk))+...
    U_u(IX,jj,kk).*(C(IX,jj,kk+1)-C(IX,jj,kk))+...
    U_d(IX,jj,kk).*(C(IX,jj,kk-1)-C(IX,jj,kk)) )+...
    wb*Dt/Dspace/sqrt(2)*(...
    U_nw(IX,jj,kk).*(C(IX-1,jj+1,kk)-C(IX,jj,kk))+...
    U_sw(IX,jj,kk).*(C(IX-1,jj-1,kk-1)-C(IX,jj,kk))  )+...
    Source(IX,jj,kk)*Dt;

% West Boundary  (X=1)
New_C(1,jj,kk)=C(1,jj,kk)+...
    wa*Dt/Dspace.*(...
    U_s(1,jj,kk).*(C(1,jj-1,kk)-C(1,jj,kk))+...
    U_e(1,jj,kk).*(C(2,jj,kk)-C(1,jj,kk))+...
    U_n(1,jj,kk).*(C(1,jj+1,kk)-C(1,jj,kk))+...
    U_u(1,jj,kk).*(C(1,jj,kk+1)-C(1,jj,kk))+...
    U_d(1,jj,kk).*(C(1,jj,kk-1)-C(1,jj,kk)) )+...
    wb*Dt/Dspace/sqrt(2)*(...
    U_ne(1,jj,kk).*(C(2,jj+1,kk)-C(1,jj,kk))+...
    U_se(1,jj,kk).*(C(2,jj-1,kk)-C(1,jj,kk))  )+...
    Source(1,jj,kk)*Dt;

% Upper Boundary (Z=IZ)
New_C(ii,jj,IZ)=C(ii,jj,IZ)+...
    wa*Dt/Dspace*(...
    U_n(ii,jj,IZ).*(C(ii,jj+1,IZ)-C(ii,jj,IZ))+...
    U_s(ii,jj,IZ).*(C(ii,jj-1,IZ)-C(ii,jj,IZ))+...
    U_e(ii,jj,IZ).*(C(ii+1,jj,IZ)-C(ii,jj,IZ))+...
    U_w(ii,jj,IZ).*(C(ii-1,jj,IZ)-C(ii,jj,IZ))+...
    U_d(ii,jj,IZ).*(C(ii,jj,IZ-1)-C(ii,jj,IZ)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_ne(ii,jj,IZ).*(C(ii+1,jj+1,IZ)-C(ii,jj,IZ))+...
    U_nw(ii,jj,IZ).*(C(ii-1,jj+1,IZ)-C(ii,jj,IZ))+...
    U_sw(ii,jj,IZ).*(C(ii-1,jj-1,IZ)-C(ii,jj,IZ))+...
    U_se(ii,jj,IZ).*(C(ii+1,jj-1,IZ)-C(ii,jj,IZ))  )+...
    Source(ii,jj,IZ)*Dt;

% Lower Boundary  (Z=1)
New_C(ii,jj,1)=C(ii,jj,1)+...
    wa*Dt/Dspace*(...
    U_n(ii,jj,1).*(C(ii,jj+1,1)-C(ii,jj,1))+...
    U_s(ii,jj,1).*(C(ii,jj-1,1)-C(ii,jj,1))+...
    U_e(ii,jj,1).*(C(ii+1,jj,1)-C(ii,jj,1))+...
    U_w(ii,jj,1).*(C(ii-1,jj,1)-C(ii,jj,1))+...
    U_u(ii,jj,1).*(C(ii,jj,1+1)-C(ii,jj,1))  )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_ne(ii,jj,1).*(C(ii+1,jj+1,1)-C(ii,jj,1))+...
    U_nw(ii,jj,1).*(C(ii-1,jj+1,1)-C(ii,jj,1))+...
    U_sw(ii,jj,1).*(C(ii-1,jj-1,1)-C(ii,jj,1))+...
    U_se(ii,jj,1).*(C(ii+1,jj-1,1)-C(ii,jj,1))  )+...
    Source(ii,jj,1)*Dt;


% SW Lower Corner  (1,1,1)
New_C(1,1,1)=C(1,1,1)+...
    wa*Dt/Dspace.*(...
    U_n(1,1,1).*(C(1,2,1)-C(1,1,1))+...
    U_e(1,1,1).*(C(2,1,1)-C(1,1,1))+...
    U_u(1,1,1).*(C(1,1,2)-C(1,1,1)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_ne(1,1,1).*(C(2,2,1)-C(1,1,1))  )+...
    Source(1,1,1)*Dt;

% NE Lower Corner (IX,IY,1)
New_C(IX,IY,1)=C(IX,IY,1)+...
    wa*Dt/Dspace.*(...
    U_s(IX,IY,1).*(C(IX,IY-1,1)-C(IX,IY,1))+...
    U_w(IX,IY,1).*(C(IX-1,IY,1)-C(IX,IY,1))+...
    U_u(IX,IY,1).*(C(IX,IY,2)-C(IX,IY,1)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_sw(IX,IY,1).*(C(IX-1,IY-1,1)-C(IX,IY,1))  )+...
    Source(IX,IY,1)*Dt;

% NW Lower Corner  (1,IY,1)
New_C(1,IY,1)=C(1,IY,1)+...
    wa*Dt/Dspace.*(...
    U_s(1,IY,1).*(C(1,IY-1,1)-C(1,IY,1))+...
    U_e(1,IY,1).*(C(2,IY,1)-C(1,IY,1))+...
    U_u(1,IY,1).*(C(1,IY,2)-C(1,IY,1))  )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_se(1,IY,1).*(C(2,IY-1,1)-C(1,IY,1))  )+...
    Source(1,IY,1)*Dt;

% SE Lower Corner  (IX,1,1)
New_C(IX,1,1)=C(IX,1,1)+...
    wa*Dt/Dspace.*(...
    U_n(IX,1,1).*(C(IX,2,1)-C(IX,1,1))+...
    U_w(IX,1,1).*(C(IX-1,1,1)-C(IX,1,1))+...
    U_u(IX,1,1).*(C(IX,1,2)-C(IX,1,1))  )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_nw(IX,1,1).*(C(IX-1,2,1)-C(IX,1,1))  )+...
    Source(IX,1,1)*Dt;

% SW Upper Corner  (1,1,IZ)
New_C(1,1,IZ)=C(1,1,IZ)+...
    wa*Dt/Dspace.*(...
    U_n(1,1,IZ).*(C(1,2,IZ)-C(1,1,IZ))+...
    U_e(1,1,IZ).*(C(2,1,IZ)-C(1,1,IZ))+...
    U_d(1,1,IZ).*(C(1,1,IZ-1)-C(1,1,IZ)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_ne(1,1,IZ).*(C(2,2,IZ)-C(1,1,IZ))  )+...
    Source(1,1,IZ)*Dt;

% NE Upper Corner (IX,IY,IZ)
New_C(IX,IY,IZ)=C(IX,IY,IZ)+...
    wa*Dt/Dspace.*(...
    U_s(IX,IY,IZ).*(C(IX,IY-1,IZ)-C(IX,IY,IZ))+...
    U_w(IX,IY,IZ).*(C(IX-1,IY,IZ)-C(IX,IY,IZ))+...
    U_d(IX,IY,IZ).*(C(IX,IY,IZ-1)-C(IX,IY,IZ)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_sw(IX,IY,IZ).*(C(IX-1,IY-1,IZ)-C(IX,IY,IZ))  )+...
    Source(IX,IY,IZ)*Dt;

% NW Upper Corner  (1,IY,IZ)
New_C(1,IY,IZ)=C(1,IY,IZ)+...
    wa*Dt/Dspace.*(...
    U_s(1,IY,IZ).*(C(1,IY-1,IZ)-C(1,IY,IZ))+...
    U_e(1,IY,IZ).*(C(2,IY,IZ)-C(1,IY,IZ))+...
    U_d(1,IY,IZ).*(C(1,IY,IZ-1)-C(1,IY,IZ))  )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_se(1,IY,IZ).*(C(2,IY-1,IZ)-C(1,IY,IZ))  )+...
    Source(1,IY,IZ)*Dt;

% SE Upper Corner  (IX,1,IZ)
New_C(IX,1,IZ)=C(IX,1,IZ)+...
    wa*Dt/Dspace.*(...
    U_n(IX,1,IZ).*(C(IX,2,IZ)-C(IX,1,IZ))+...
    U_w(IX,1,IZ).*(C(IX-1,1,IZ)-C(IX,1,IZ))+...
    U_d(IX,1,IZ).*(C(IX,1,IZ-1)-C(IX,1,IZ))  )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_nw(IX,1,IZ).*(C(IX-1,2,IZ)-C(IX,1,IZ))  )+...
    Source(IX,1,IZ)*Dt;

% ============================================================
%% Diffusion considered
% =========================================================
elseif diffusion_flag == 1
    
% Local update in main body
New_C(ii,jj,kk)=C(ii,jj,kk)+...
    Kx*Dt/Dspace^2*(C(ii+1,jj,kk)+C(ii-1,jj,kk)-2*C(ii,jj,kk))+...
    Ky*Dt/Dspace^2*(C(ii,jj+1,kk)+C(ii,jj-1,kk)-2*C(ii,jj,kk))+...
    Kxy*Dt/Dspace^2/2*(C(ii+1,jj-1,kk)+C(ii-1,jj+1,kk)-2*C(ii,jj,kk))+...
    Kyx*Dt/Dspace^2/2*(C(ii-1,jj-1,kk)+C(ii+1,jj+1,kk)-2*C(ii,jj,kk))+...
    Kz*Dt/Dspace^2*(C(ii,jj,kk+1)+C(ii,jj,kk-1)-2*C(ii,jj,kk))+...
    wa*Dt/Dspace*(...
    U_n(ii,jj,kk).*(C(ii,jj+1,kk)-C(ii,jj,kk))+...
    U_s(ii,jj,kk).*(C(ii,jj-1,kk)-C(ii,jj,kk))+...
    U_e(ii,jj,kk).*(C(ii+1,jj,kk)-C(ii,jj,kk))+...
    U_w(ii,jj,kk).*(C(ii-1,jj,kk)-C(ii,jj,kk))+...
    U_u(ii,jj,kk).*(C(ii,jj,kk+1)-C(ii,jj,kk))+...
    U_d(ii,jj,kk).*(C(ii,jj,kk-1)-C(ii,jj,kk)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_ne(ii,jj,kk).*(C(ii+1,jj+1,kk)-C(ii,jj,kk))+...
    U_nw(ii,jj,kk).*(C(ii-1,jj+1,kk)-C(ii,jj,kk))+...
    U_sw(ii,jj,kk).*(C(ii-1,jj-1,kk)-C(ii,jj,kk))+...
    U_se(ii,jj,kk).*(C(ii+1,jj-1,kk)-C(ii,jj,kk))  )+...
    Source(ii,jj,kk)*Dt;

%%%%%%%%% Face Boundary Conditions  %%%%%%%%%%%%%%%%%%%%%%%%%%
% South Boundary  (Y=1)
New_C(ii,1,kk)=C(ii,1,kk)+...
    Kx*Dt/Dspace^2*(C(ii+1,1,kk)+C(ii-1,1,kk)-2*C(ii,1,kk))+...
    Kxy*Dt/Dspace^2/2*(C(ii-1,1+1,kk)-2*C(ii,1,kk))+...
    Kyx*Dt/Dspace^2/2*(C(ii+1,1+1,kk)-2*C(ii,1,kk))+...
    Kz*Dt/Dspace^2*(C(ii,1,kk+1)+C(ii,1,kk-1)-2*C(ii,1,kk))+...
    wa*Dt/Dspace*(...
    U_n(ii,1,kk).*(C(ii,2,kk)-C(ii,1,kk))+...
    U_e(ii,1,kk).*(C(ii+1,1,kk)-C(ii,1,kk))+...
    U_w(ii,1,kk).*(C(ii-1,1,kk)-C(ii,1,kk))+...
    U_u(ii,1,kk).*(C(ii,1,kk+1)-C(ii,1,kk))+...
    U_d(ii,1,kk).*(C(ii,1,kk-1)-C(ii,1,kk)) )+...
    wb*Dt/Dspace/sqrt(2)*(...
    U_ne(ii,1,kk).*(C(ii+1,2,kk)-C(ii,1,kk))+...
    U_nw(ii,1,kk).*(C(ii-1,2,kk)-C(ii,1,kk))  )+...
    Source(ii,1,kk)*Dt;

% North Boundary  (Y=IY)
New_C(ii,IY,kk)=C(ii,IY,kk)+...
        Kx*Dt/Dspace^2*(C(ii+1,IY,kk)+C(ii-1,IY,kk)-2*C(ii,IY,kk))+...
    Ky*Dt/Dspace^2*(C(ii,IY-1,kk)-2*C(ii,IY,kk))+...
    Kxy*Dt/Dspace^2/2*(C(ii+1,IY-1,kk)-2*C(ii,IY,kk))+...
    Kyx*Dt/Dspace^2/2*(C(ii-1,IY-1,kk)-2*C(ii,IY,kk))+...
    Kz*Dt/Dspace^2*(C(ii,IY,kk+1)+C(ii,IY,kk-1)-2*C(ii,IY,kk))+...
    wa*Dt/Dspace*(...
    U_s(ii,IY,kk).*(C(ii,IY-1,kk)-C(ii,IY,kk))+...
    U_e(ii,IY,kk).*(C(ii+1,IY,kk)-C(ii,IY,kk))+...
    U_w(ii,IY,kk).*(C(ii-1,IY,kk)-C(ii,IY,kk))+...
    U_u(ii,IY,kk).*(C(ii,IY,kk+1)-C(ii,IY,kk))+...
    U_d(ii,IY,kk).*(C(ii,IY,kk-1)-C(ii,IY,kk)) )+...
    wb*Dt/Dspace/sqrt(2)*(...
    U_se(ii,IY,kk).*(C(ii+1,IY-1,kk)-C(ii,IY,kk))+...
    U_sw(ii,IY,kk).*(C(ii-1,IY-1,kk)-C(ii,IY,kk))  )+...
    Source(ii,IY,kk)*Dt;

% East Boundary  (X=IX)
New_C(IX,jj,kk)=C(IX,jj,kk)+...
        Kx*Dt/Dspace^2*(C(IX-1,jj,kk)-2*C(IX,jj,kk))+...
    Ky*Dt/Dspace^2*(C(IX,jj+1,kk)+C(IX,jj-1,kk)-2*C(IX,jj,kk))+...
    Kxy*Dt/Dspace^2/2*(+C(IX-1,jj+1,kk)-2*C(IX,jj,kk))+...
    Kyx*Dt/Dspace^2/2*(C(IX-1,jj-1,kk)-2*C(IX,jj,kk))+...
    Kz*Dt/Dspace^2*(C(IX,jj,kk+1)+C(IX,jj,kk-1)-2*C(IX,jj,kk))+...
    wa*Dt/Dspace.*(...
    U_s(IX,jj,kk).*(C(IX,jj-1,kk)-C(IX,jj,kk))+...
    U_w(IX,jj,kk).*(C(IX-1,jj,kk)-C(IX,jj,kk))+...
    U_n(IX,jj,kk).*(C(IX,jj+1,kk)-C(IX,jj,kk))+...
    U_u(IX,jj,kk).*(C(IX,jj,kk+1)-C(IX,jj,kk))+...
    U_d(IX,jj,kk).*(C(IX,jj,kk-1)-C(IX,jj,kk)) )+...
    wb*Dt/Dspace/sqrt(2)*(...
    U_nw(IX,jj,kk).*(C(IX-1,jj+1,kk)-C(IX,jj,kk))+...
    U_sw(IX,jj,kk).*(C(IX-1,jj-1,kk-1)-C(IX,jj,kk))  )+...
    Source(IX,jj,kk)*Dt;

% West Boundary  (X=1)
New_C(1,jj,kk)=C(1,jj,kk)+...
        Kx*Dt/Dspace^2*(C(1+1,jj,kk)-2*C(1,jj,kk))+...
    Ky*Dt/Dspace^2*(C(1,jj+1,kk)+C(1,jj-1,kk)-2*C(1,jj,kk))+...
    Kxy*Dt/Dspace^2/2*(C(1+1,jj-1,kk)-2*C(1,jj,kk))+...
    Kyx*Dt/Dspace^2/2*(C(1+1,jj+1,kk)-2*C(1,jj,kk))+...
    Kz*Dt/Dspace^2*(C(1,jj,kk+1)+C(1,jj,kk-1)-2*C(1,jj,kk))+...
    wa*Dt/Dspace.*(...
    U_s(1,jj,kk).*(C(1,jj-1,kk)-C(1,jj,kk))+...
    U_e(1,jj,kk).*(C(2,jj,kk)-C(1,jj,kk))+...
    U_n(1,jj,kk).*(C(1,jj+1,kk)-C(1,jj,kk))+...
    U_u(1,jj,kk).*(C(1,jj,kk+1)-C(1,jj,kk))+...
    U_d(1,jj,kk).*(C(1,jj,kk-1)-C(1,jj,kk)) )+...
    wb*Dt/Dspace/sqrt(2)*(...
    U_ne(1,jj,kk).*(C(2,jj+1,kk)-C(1,jj,kk))+...
    U_se(1,jj,kk).*(C(2,jj-1,kk)-C(1,jj,kk))  )+...
    Source(1,jj,kk)*Dt;

% Upper Boundary (Z=IZ)
New_C(ii,jj,IZ)=C(ii,jj,IZ)+...
        Kx*Dt/Dspace^2*(C(ii+1,jj,IZ)+C(ii-1,jj,IZ)-2*C(ii,jj,IZ))+...
    Ky*Dt/Dspace^2*(C(ii,jj+1,IZ)+C(ii,jj-1,IZ)-2*C(ii,jj,IZ))+...
    Kxy*Dt/Dspace^2/2*(C(ii+1,jj-1,IZ)+C(ii-1,jj+1,IZ)-2*C(ii,jj,IZ))+...
    Kyx*Dt/Dspace^2/2*(C(ii-1,jj-1,IZ)+C(ii+1,jj+1,IZ)-2*C(ii,jj,IZ))+...
    Kz*Dt/Dspace^2*(C(ii,jj,IZ-1)-2*C(ii,jj,IZ))+...
    wa*Dt/Dspace*(...
    U_n(ii,jj,IZ).*(C(ii,jj+1,IZ)-C(ii,jj,IZ))+...
    U_s(ii,jj,IZ).*(C(ii,jj-1,IZ)-C(ii,jj,IZ))+...
    U_e(ii,jj,IZ).*(C(ii+1,jj,IZ)-C(ii,jj,IZ))+...
    U_w(ii,jj,IZ).*(C(ii-1,jj,IZ)-C(ii,jj,IZ))+...
    U_d(ii,jj,IZ).*(C(ii,jj,IZ-1)-C(ii,jj,IZ)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_ne(ii,jj,IZ).*(C(ii+1,jj+1,IZ)-C(ii,jj,IZ))+...
    U_nw(ii,jj,IZ).*(C(ii-1,jj+1,IZ)-C(ii,jj,IZ))+...
    U_sw(ii,jj,IZ).*(C(ii-1,jj-1,IZ)-C(ii,jj,IZ))+...
    U_se(ii,jj,IZ).*(C(ii+1,jj-1,IZ)-C(ii,jj,IZ))  )+...
    Source(ii,jj,IZ)*Dt;

% Lower Boundary  (Z=1)
New_C(ii,jj,1)=C(ii,jj,1)+...
        Kx*Dt/Dspace^2*(C(ii+1,jj,1)+C(ii-1,jj,1)-2*C(ii,jj,1))+...
    Ky*Dt/Dspace^2*(C(ii,jj+1,1)+C(ii,jj-1,1)-2*C(ii,jj,1))+...
    Kxy*Dt/Dspace^2/2*(C(ii+1,jj-1,1)+C(ii-1,jj+1,1)-2*C(ii,jj,1))+...
    Kyx*Dt/Dspace^2/2*(C(ii-1,jj-1,1)+C(ii+1,jj+1,1)-2*C(ii,jj,1))+...
    Kz*Dt/Dspace^2*(C(ii,jj,1+1)-2*C(ii,jj,1))+...
    wa*Dt/Dspace*(...
    U_n(ii,jj,1).*(C(ii,jj+1,1)-C(ii,jj,1))+...
    U_s(ii,jj,1).*(C(ii,jj-1,1)-C(ii,jj,1))+...
    U_e(ii,jj,1).*(C(ii+1,jj,1)-C(ii,jj,1))+...
    U_w(ii,jj,1).*(C(ii-1,jj,1)-C(ii,jj,1))+...
    U_u(ii,jj,1).*(C(ii,jj,1+1)-C(ii,jj,1))  )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_ne(ii,jj,1).*(C(ii+1,jj+1,1)-C(ii,jj,1))+...
    U_nw(ii,jj,1).*(C(ii-1,jj+1,1)-C(ii,jj,1))+...
    U_sw(ii,jj,1).*(C(ii-1,jj-1,1)-C(ii,jj,1))+...
    U_se(ii,jj,1).*(C(ii+1,jj-1,1)-C(ii,jj,1))  )+...
    Source(ii,jj,1)*Dt;


% SW Lower Corner  (1,1,1)
New_C(1,1,1)=C(1,1,1)+...
    Kx*Dt/Dspace^2*(C(1+1,1,1)-2*C(1,1,1))+...
    Ky*Dt/Dspace^2*(C(1,1+1,1)-2*C(1,1,1))+...
    Kxy*Dt/Dspace^2/2*(-2*C(1,1,1))+...
    Kyx*Dt/Dspace^2/2*(+C(1+1,1+1,1)-2*C(1,1,1))+...
    Kz*Dt/Dspace^2*(C(1,1,1+1)-2*C(1,1,1))+...
    wa*Dt/Dspace.*(...
    U_n(1,1,1).*(C(1,2,1)-C(1,1,1))+...
    U_e(1,1,1).*(C(2,1,1)-C(1,1,1))+...
    U_u(1,1,1).*(C(1,1,2)-C(1,1,1)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_ne(1,1,1).*(C(2,2,1)-C(1,1,1))  )+...
    Source(1,1,1)*Dt;

% NE Lower Corner (IX,IY,1)
New_C(IX,IY,1)=C(IX,IY,1)+...
    Kx*Dt/Dspace^2*(C(IX-1,IY,1)-2*C(IX,IY,1))+...
    Ky*Dt/Dspace^2*(C(IX,IY-1,1)-2*C(IX,IY,1))+...
    Kxy*Dt/Dspace^2/2*(-2*C(IX,IY,1))+...
    Kyx*Dt/Dspace^2/2*(C(IX-1,IY-1,1)-2*C(IX,IY,1))+...
    Kz*Dt/Dspace^2*(C(IX,IY,1+1)-2*C(IX,IY,1))+...
    wa*Dt/Dspace.*(...
    U_s(IX,IY,1).*(C(IX,IY-1,1)-C(IX,IY,1))+...
    U_w(IX,IY,1).*(C(IX-1,IY,1)-C(IX,IY,1))+...
    U_u(IX,IY,1).*(C(IX,IY,2)-C(IX,IY,1)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_sw(IX,IY,1).*(C(IX-1,IY-1,1)-C(IX,IY,1))  )+...
    Source(IX,IY,1)*Dt;

% NW Lower Corner  (1,IY,1)
New_C(1,IY,1)=C(1,IY,1)+...
        Kx*Dt/Dspace^2*(C(1+1,IY,1)-2*C(1,IY,1))+...
    Ky*Dt/Dspace^2*(C(1,IY-1,1)-2*C(1,IY,1))+...
    Kxy*Dt/Dspace^2/2*(C(1+1,IY-1,1)-2*C(1,IY,1))+...
    Kyx*Dt/Dspace^2/2*(-2*C(1,IY,1))+...
    Kz*Dt/Dspace^2*(C(1,IY,1+1)-2*C(1,IY,1))+...
    wa*Dt/Dspace.*(...
    U_s(1,IY,1).*(C(1,IY-1,1)-C(1,IY,1))+...
    U_e(1,IY,1).*(C(2,IY,1)-C(1,IY,1))+...
    U_u(1,IY,1).*(C(1,IY,2)-C(1,IY,1))  )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_se(1,IY,1).*(C(2,IY-1,1)-C(1,IY,1))  )+...
    Source(1,IY,1)*Dt;

% SE Lower Corner  (IX,1,1)
New_C(IX,1,1)=C(IX,1,1)+...
        Kx*Dt/Dspace^2*(C(IX-1,1,1)-2*C(IX,1,1))+...
    Ky*Dt/Dspace^2*(C(IX,1+1,1)-2*C(IX,1,1))+...
    Kxy*Dt/Dspace^2/2*(C(IX-1,1+1,1)-2*C(IX,1,1))+...
    Kyx*Dt/Dspace^2/2*(-2*C(IX,1,1))+...
    Kz*Dt/Dspace^2*(C(IX,1,1+1)-2*C(IX,1,1))+...
    wa*Dt/Dspace.*(...
    U_n(IX,1,1).*(C(IX,2,1)-C(IX,1,1))+...
    U_w(IX,1,1).*(C(IX-1,1,1)-C(IX,1,1))+...
    U_u(IX,1,1).*(C(IX,1,2)-C(IX,1,1))  )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_nw(IX,1,1).*(C(IX-1,2,1)-C(IX,1,1))  )+...
    Source(IX,1,1)*Dt;

% SW Upper Corner  (1,1,IZ)
New_C(1,1,IZ)=C(1,1,IZ)+...
        Kx*Dt/Dspace^2*(C(1+1,1,IZ)-2*C(1,1,IZ))+...
    Ky*Dt/Dspace^2*(C(1,1+1,IZ)-2*C(1,1,IZ))+...
    Kxy*Dt/Dspace^2/2*(-2*C(1,1,IZ))+...
    Kyx*Dt/Dspace^2/2*(C(1+1,1+1,IZ)-2*C(1,1,IZ))+...
    Kz*Dt/Dspace^2*(C(1,1,IZ-1)-2*C(1,1,IZ))+...
    wa*Dt/Dspace.*(...
    U_n(1,1,IZ).*(C(1,2,IZ)-C(1,1,IZ))+...
    U_e(1,1,IZ).*(C(2,1,IZ)-C(1,1,IZ))+...
    U_d(1,1,IZ).*(C(1,1,IZ-1)-C(1,1,IZ)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_ne(1,1,IZ).*(C(2,2,IZ)-C(1,1,IZ))  )+...
    Source(1,1,IZ)*Dt;

% NE Upper Corner (IX,IY,IZ)
New_C(IX,IY,IZ)=C(IX,IY,IZ)+...
        Kx*Dt/Dspace^2*(C(IX-1,IY,IZ)-2*C(IX,IY,IZ))+...
    Ky*Dt/Dspace^2*(C(IX,IY-1,IZ)-2*C(IX,IY,IZ))+...
    Kxy*Dt/Dspace^2/2*(-2*C(IX,IY,IZ))+...
    Kyx*Dt/Dspace^2/2*(C(IX-1,IY-1,IZ)-2*C(IX,IY,IZ))+...
    Kz*Dt/Dspace^2*(C(IX,IY,IZ-1)-2*C(IX,IY,IZ))+...
    wa*Dt/Dspace.*(...
    U_s(IX,IY,IZ).*(C(IX,IY-1,IZ)-C(IX,IY,IZ))+...
    U_w(IX,IY,IZ).*(C(IX-1,IY,IZ)-C(IX,IY,IZ))+...
    U_d(IX,IY,IZ).*(C(IX,IY,IZ-1)-C(IX,IY,IZ)) )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_sw(IX,IY,IZ).*(C(IX-1,IY-1,IZ)-C(IX,IY,IZ))  )+...
    Source(IX,IY,IZ)*Dt;

% NW Upper Corner  (1,IY,IZ)
New_C(1,IY,IZ)=C(1,IY,IZ)+...
        Kx*Dt/Dspace^2*(C(1+1,IY,IZ)-2*C(1,IY,IZ))+...
    Ky*Dt/Dspace^2*(C(1,IY-1,IZ)-2*C(1,IY,IZ))+...
    Kxy*Dt/Dspace^2/2*(C(1+1,IY-1,IZ)-2*C(1,IY,IZ))+...
    Kyx*Dt/Dspace^2/2*(-2*C(1,IY,IZ))+...
    Kz*Dt/Dspace^2*(C(1,IY,IZ-1)-2*C(1,IY,IZ))+...
    wa*Dt/Dspace.*(...
    U_s(1,IY,IZ).*(C(1,IY-1,IZ)-C(1,IY,IZ))+...
    U_e(1,IY,IZ).*(C(2,IY,IZ)-C(1,IY,IZ))+...
    U_d(1,IY,IZ).*(C(1,IY,IZ-1)-C(1,IY,IZ))  )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_se(1,IY,IZ).*(C(2,IY-1,IZ)-C(1,IY,IZ))  )+...
    Source(1,IY,IZ)*Dt;

% SE Upper Corner  (IX,1,IZ)
New_C(IX,1,IZ)=C(IX,1,IZ)+...
        Kx*Dt/Dspace^2*(C(IX-1,1,IZ)-2*C(IX,1,IZ))+...
    Ky*Dt/Dspace^2*(C(IX,1+1,IZ)-2*C(IX,1,IZ))+...
    Kxy*Dt/Dspace^2/2*(C(IX-1,1+1,IZ)-2*C(IX,1,IZ))+...
    Kyx*Dt/Dspace^2/2*(-2*C(IX,1,IZ))+...
    Kz*Dt/Dspace^2*(C(IX,1,IZ-1)-2*C(IX,1,IZ))+...
    wa*Dt/Dspace.*(...
    U_n(IX,1,IZ).*(C(IX,2,IZ)-C(IX,1,IZ))+...
    U_w(IX,1,IZ).*(C(IX-1,1,IZ)-C(IX,1,IZ))+...
    U_d(IX,1,IZ).*(C(IX,1,IZ-1)-C(IX,1,IZ))  )+...
    wb*Dt/Dspace/sqrt(2)*( ...
    U_nw(IX,1,IZ).*(C(IX-1,2,IZ)-C(IX,1,IZ))  )+...
    Source(IX,1,IZ)*Dt;


    
    
end


