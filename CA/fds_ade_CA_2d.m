%% Script Description:
% Each block is marked as essential or un-necessary for continious execuation
% Essential block is marked as #<Exec: %i>#, where %i indicates the sequence
% Un-necessary block is marked as <optional@%i> where %i indicates the block that should be
%       running before the optional block could run
%% #<Exec:01>#   Load fdscoverted data from hard disk
clear
clc

% Windows PATH
% addpath('d:/fdscov')

% Linux PATH
addpath('~/fdscov')
fid='029';
load(['Spectra_',fid,'.mat'])
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
    u_vel_qver=reshape(u_vel(:,:,1,i),IX,IY);
    v_vel_qver=reshape(v_vel(:,:,1,i),IX,IY);
    quiver(X0,Y0,u_vel_qver,v_vel_qver)
    title(['Time = ',num2str(round(time(i))),' (s)'])
    drawnow
end
close(1)
pause(2)

% Plot Original Chlorine Concentration 
fun_plot2D_facility(CHLORINE_VOLUME_FRACTION,time,3)

for i=1:IT
    imagesc(CHLORINE_VOLUME_FRACTION(:,:,1,i)')
    title(['Progress ',num2str(round(i/IT*100)),'%'])
    ax=gca;
    ax.YDir='normal';
    colormap(jet)
    drawnow
end

%% #<Exec:02>#    Some statistics
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
tic;Illustration
Dt=0.06;
Dspace=1;
i=1;j=1;
time_q=0:Dt:360;
ITq=length(time_q);

textprogressbar('Interpolating : ');
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
        
        textprogressbar(round((IX-i)/IX*100));
    end
end
textprogressbar('  Done');

time_for_interp=toc;   % about two minites

save(['~/fdscov/Spectra_',fid,'_q.mat'],'time_q','u_vel_q','v_vel_q','vel_q',...
    'CHLORINE_VOLUME_FRACTION');


%% #<Spectial:01>#  Load saved interpolation data if exist
clear
clc
addpath('~/fdscov')
fid='029';

file_exist=system(['test -e ~/fdscov/Spectra_',fid,'_q.mat']);
if (file_exist == 0)
    load(['~/fdscov/Spectra_',fid,'_q.mat'])
else
    disp(['File ~/fdscov/Spectra_',fid,'_q.mat not exist!'])
    return
end
[IX,IY,IT]=size(u_vel_q);
Dspace=1.0;
idum=randi(IT);
Dt=time_q(idum)-time_q(idum-1);

%% #<Exec:05># / #<Spectial:01#    Begin to simulate CA dispersion
t=0;            % time start
T_end=360;      % time of simulation

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
    
    %
    u_vel_t=reshape(u_vel_q(:,:,1,step),IX,IY);
    v_vel_t=reshape(v_vel_q(:,:,1,step),IX,IY);
    vel_t=reshape(vel_q(:,:,1,step),IX,IY);
    
    New_C=fun_update2D(C,u_vel_t,v_vel_t,Dt,Dspace,Source);
    
    C=New_C;
    
    if t >= 60 
        Source(235:240,90:95)=5;
    end
    
    if t >= 280
        Source=zeros(IX,IY);
    end
    
    con(:,:,step)=New_C;
    
    textprogressbar(round(t/T_end*1000)/10)
        
end
textprogressbar(' Finished')
CA_loop_time=toc;

 save(['~/fdscov/Spectra_',fid,'_con.mat'],'con','time_q');

%% #<Exec:06>#  Results visulization
clear
fid='029';
load(['~/fdscov/Spectra_',fid,'_con.mat'],'con');
     
     
fun_plot2D_facility(con,time_q,2);
% [IX,IY,IT]=size(con);
% v=[5,1,0.5,0.1,0.05];  % contour level
% for k=4000:IT
%     imagesc(X_cor,Y_cor,con(:,:,k)');
%     title(['Time eclapse: ',num2str(round(time_q(k))),' s'])
%     colormap(jet);
%         ax = gca;
%         ax.YDir='normal';
%     set(gca,'YDir','normal');
%     fun_plot2D_facility
%     drawnow
% end


%% #<optional@non>#  Require adding source to path
% demo_textprogressbar
% This a demo for textprogressbar script
textprogressbar('calculating outputs: ');
for i=1:100
    textprogressbar(i);
    pause(0.1);
end
textprogressbar('done');


textprogressbar('saving data:         ');
for i=1:0.5:80
    textprogressbar(i);
    pause(0.05);
end
textprogressbar('terminated');

%% #<optional@non># Require Adding source to path
pb = CmdLineProgressBar('Doing stuff...');
for k = 1 : 100
    pb.print(k,100)
    pause(0.1)
    % do stuff
end
