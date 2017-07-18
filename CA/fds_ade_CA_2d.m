%% Load fdscoverted data from hard disk
clear
clc
addpath('d:/fdscov')
load('Spectra_001.mat')
u_vel=U_VELOCITY;
v_vel=V_VELOCITY;
w_vel=W_VELOCITY;
vel=VELOCITY;
clear U_VELOCITY V_VELOCITY W_VELOCITY
[IX,IY,IZ,IT]=size(vel);
IT=IT-1;
time=time(1:IT);

%% Illustration of wind velocity vector to determine the actual direction
[Y0,X0]=meshgrid(-60:120,-60:275);
for i=1:IT
    u_vel_qver=reshape(u_vel(:,:,1,i),IX,IY);
    v_vel_qver=reshape(v_vel(:,:,1,i),IX,IY);
    quiver(X0,Y0,u_vel_qver,v_vel_qver)
    title(['Time = ',num2str(time(i)),' (s)'])
    pause(1)
end


%% Some statistics
% time space analysis
Dtimearray=time(2:IT)-time(1:IT-1);
plot(Dtimearray)
Dt_test=mean(Dtimearray);
pause(3)
% mean wind velocity at time t
for i=IT:-1:1
    mean_vel_array(i)=mean(mean(vel(:,:,1,i)));
end
plot(mean_vel_array)
vel_mean_time=mean(mean_vel_array);

%% Determine Dt and Dspace given specific Dispersion coefficient
% Rule: mean_vel*Dt/Dspace <= 1/2
Dspace_test=1.0;
% 1. Given Dspace
Dt_upperbound=Dspace_test/2/vel_mean_time;
disp(['Maximum time interval is: <',num2str(Dt_upperbound),' given space interval = ',num2str(Dspace_test)])
% 2. Given Dtime
Dspace_lowbound=2*vel_mean_time*Dt_test;
disp(['Minimum space interval is: >',num2str(Dspace_lowbound),' given time interval = ',num2str(Dt_test)])

%% Velocity field space interpolation
Dt=0.36;
Dspace=
[X1,Y1]=meshgrid(-60:275,-60:120);
Dspace1=0.25;
[X1,Y1]=meshgrid(-60:Dspace1:275,-60:Dspace1:120);
vel0=vel(:,:,1,100);

% for it=IT:-1:1
vel1=interp2(X0,Y0,vel0',X1,Y1,'cubic');
% end
imagesc(vel1)


%% Velocity field time interpolation
tic;
Dt=0.1;
Dspace=1;
i=1;j=1;
time_q=0:Dt:360;
ITq=length(time_q);
for i=IX:-1:1
    for j=IY:-1:1
        u_vel_temp=reshape(u_vel(i,j,1,1:IT),1,IT);
        v_vel_temp=reshape(v_vel(i,j,1,1:IT),1,IT);
        
        u_vel_new=spline(time,u_vel_temp,time_q);
        v_vel_new=spline(time,v_vel_temp,time_q);
        vel_new=sqrt(u_vel_new.^2+v_vel_new.^2);
        
        u_vel_q(i,j,1,:)=reshape(u_vel_new,1,1,1,ITq);
        v_vel_q(i,j,1,:)=reshape(v_vel_new,1,1,1,ITq);
        vel_q(i,j,1,:)=reshape(vel_new,1,1,1,ITq);
    end
end
time_for_interp=toc;

save('d:/fdscov/Spectra_001_q.mat','time_q','u_vel_q','v_vel_q','vel_q',...
    'CHLORINE_VOLUME_FRACTION');




%% Begin to simulate CA dispersion
t=0;            % time start
Dt=0.1;         % time interval
T_end=360;      % time of simulation

% Initial condition
if t == 0

C=zeros(IX,IY);      % Array initialization
Source=zeros(IX,IY);      % Source Character
Source(179:180,120,121)=10;
end
step=0;

% Begin Loop
while t<=T_end
    t=t+Dt;
    step=step+1;
    
    % 
    u_vel_t=reshape(u_vel_q(:,:,1,step),IX,IY);
    v_vel_t=reshape(v_vel_q(:,:,1,step),IX,IY);
    vel_t=reshape(vel_q(:,:,1,step),IX,IY);
    
New_C=fun_update2D(C,parameter);


end
