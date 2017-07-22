%% Script Description:
% Each block is marked as essential or un-necessary for continious execuation
% Essential block is marked as #<Exec: %i>#, where %i indicates the sequence
% Un-necessary block is marked as <optional@%i> where %i indicates the block that should be
%       running before the optional block could run
%
% Regular variable that should be kept in workspace
clear
clc
arch='win'; % | 'linux', 'mac'
fid='055';
leak=[96,97,80,81];   % source location
% leak=[299,300,80,81];   % source location
strength=1;  % source strength
lk_start=100;   % time to leak
lk_end=200;     % time leak ends
layer=1;        % z layer
T_end=400;      % Total simulation time

if strcmp(arch,'win') 
    disp('Working on ''Windows'' platform')
    addpath('D:/fdsmat')
elseif strcmp(arch,'linux')
addpath('~/fdscov')
elseif strcmp(arch,'mac')
    addpath('~/fdscov')
end
% chid=['Spectra_',fid];
chid=['Facility_',fid];

load('mycolormap.mat')

%% #<Exec:01>#   Load fdscoverted data from hard disk
% Windows PATH
% addpath('d:/fdscov')

load([chid,'.mat'])
u_vel=U_VELOCITY;
v_vel=V_VELOCITY;
w_vel=W_VELOCITY;
vel=VELOCITY;
clear U_VELOCITY V_VELOCITY W_VELOCITY
[IX,IY,IZ,IT]=size(vel);
IT=IT-1;
time=time(1:IT);

%% <optional@01>
% Visulization of wind velocity vector to determine the actual direction

[Y0,X0]=meshgrid(-60:120,-60:275);
figure(1)
for i=1:IT
    u_vel_qver=reshape(u_vel(:,:,layer,i),IX,IY);
    v_vel_qver=reshape(v_vel(:,:,layer,i),IX,IY);
    quiver(X0,Y0,u_vel_qver,v_vel_qver)
    title(['Time = ',num2str(round(time(i))),' (s)'])
    drawnow
end
close(1)
pause(2)

%% <optional@02>    Plot Original Chlorine Concentration 

fun_plot2D_facility(CHLORINE_VOLUME_FRACTION,time,3,mycmap,chid,layer,arch)

%% #<Exec:02>#    Some statistics
% time space analysis
Dtimearray=time(2:IT)-time(1:IT-1);
plot(Dtimearray)
Dt_test=mean(Dtimearray);
pause(3)
% mean wind velocity at time t
for i=IT:-1:1
    mean_vel_array(i)=mean(mean(vel(:,:,layer,i)));
end
plot(mean_vel_array)
vel_mean_time=mean(mean_vel_array);

%% #<Exec:03>#    Determine Dt and Dspace given specific Dispersion coefficient
% Rule: mean_vel*Dt/Dspace <= 1/2
Dspace_test=1.0;
% 1. Given Dspace
Dt_upperbound=Dspace_test/2/vel_mean_time;
disp(['Maximum time interval is: <',num2str(Dt_upperbound),' given space interval = ',num2str(Dspace_test)])
% 2. Given Dtime
Dspace_lowbound=2*vel_mean_time*Dt_test;
disp(['Minimum space interval is: >',num2str(Dspace_lowbound),' given time interval = ',num2str(Dt_test)])


%% #<Exec:04>#     Velocity field time interpolation
tic;

% ------------ Dt & Dspace could change according to Exec:03
Dt=0.05;
Dspace=1;
% ------------------------------------------------------------


i=1;j=1;
time_q=0:Dt:T_end;
ITq=length(time_q);

textprogressbar('Interpolating... : ');
for i=IX:-1:1
    for j=IY:-1:1
        u_vel_temp=reshape(u_vel(i, j, layer, 1:IT),1,IT);
        v_vel_temp=reshape(v_vel(i, j, layer, 1:IT),1,IT);
        
        u_vel_new=spline(time,u_vel_temp,time_q);
        v_vel_new=spline(time,v_vel_temp,time_q);
        vel_new=sqrt(u_vel_new.^2+v_vel_new.^2);
        
        u_vel_q(i,j,layer, :)=reshape(u_vel_new,1,1,1,ITq);
        v_vel_q(i,j,layer, :)=reshape(v_vel_new,1,1,1,ITq);
        vel_q(i,j,layer,:)=reshape(vel_new,1,1,1,ITq);
        
        textprogressbar(round((IX-i)/IX*100));
    end
end
time_for_interp=toc;   % about two minites
textprogressbar([' Time Elapse: ',num2str(time_for_interp)]);

if strcmp(arch,'linux')
save(['~/fdscov/',chid,'_q.mat'],'time_q','u_vel_q','v_vel_q','vel_q',...
    'CHLORINE_VOLUME_FRACTION');
elseif strcmp(arch,'win')
save(['D:\fdsmat\',chid,'_q.mat'],'time_q','u_vel_q','v_vel_q','vel_q',...
    'CHLORINE_VOLUME_FRACTION');
end
%% #<Spectial:01>#  Load saved interpolation data if exist
%------------------------------
%
% Go to initial then go back to reload
%
% -----------------------------

if strcmp(arch, 'linux')    
    file_exist=system(['test -e ~/fdscov/',chid,'_q.mat']);
    if (file_exist == 2)
        load(['~/fdscov/',chid,'_q.mat'])
    else
        disp(['File ~/fdscov/',chid,'_q.mat doesn''t exist!'])
    end
elseif strcmp(arch, 'win')
    file_exist=exist(['D:\fdsmat\',chid,'_q.mat'],'file');
    if (file_exist == 2)
        load(['D:\fdsmat\',chid,'_q.mat'])
    else
        disp(['File D:\fdsmat\',chid,'_q.mat doesn''t exist!'])
    end
end
[IX,IY,IZ,IT]=size(u_vel_q);
Dspace=1.0;
idum=randi(IT);
Dt=time_q(idum)-time_q(idum-1);
disp(['Check, Dt is: ',num2str(Dt)]);

%% #<Exec:05># / #<Spectial:02#    Begin to simulate CA dispersion
t=0;            % time start

% Initial condition
if t == 0
    C=zeros(IX,IY);      % Array initialization
    Source=zeros(IX,IY);      % Source Character
    con=zeros(IX,IY,IT);
end
step=0;

% Begin Loop
tic;
textprogressbar('Executing CA: ')
while t<=T_end
    t=t+Dt;
    step=step+1;
    
    % Extract one layer at time step from 3-D velocity field
    u_vel_t=reshape(u_vel_q(:,:,layer,step),IX,IY);
    v_vel_t=reshape(v_vel_q(:,:,layer,step),IX,IY);
    vel_t=reshape(vel_q(:,:,layer,step),IX,IY);
    
    New_C=fun_update2D(C,u_vel_t,v_vel_t,Dt,Dspace,Source);
    
    C=New_C;
    
    if t >= lk_start && t < lk_end
        Source(leak(1):leak(2),leak(3):leak(4))=strength;
    else
        Source=zeros(IX,IY);
    end
    
    con(:,:,step)=New_C;
    
    textprogressbar(round(t/T_end*1000)/10)
        
end
CA_loop_time=toc;
textprogressbar(['Time Elapse: ',num2str(CA_loop_time),' s'])

% if strcmp(arch, 'linux')
%  save(['~/fdscov/',chid,'_con.mat'],'con','time_q');
% elseif strcmp(arch,'win')
%     save(['D:\fdsmat\',chid,'_con.mat'],'con','time_q');
% end

%% #<Exec:06>#  Results visulization
% -------------------------
%   Go to initial to clear
% ---------------------------
if strcmp(arch,'linux')
   load(['~/fdscov/',chid,'_con.mat'],'con','time_q');
elseif strcmp(arch,'win')
    load(['D:\fdsmat\',chid,'_con.mat'],'con','time_q');
end
     
% Generate CA results 
fun_plot2D_facility(con,time_q,2,mycmap,chid,layer,arch)


%% Debug purpose
iskp=100;
k=0;
seqns=1:iskp:size(con,3);
for i=seqns
    k=k+1;
    CMax(k)=max(max(con(:,:,i)));
end
plot(time_q(seqns),CMax)



