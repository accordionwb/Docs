% Dimensions
I=600; %x
J=500; %y


C=zeros(I,J); % Array initialization
C(410:420,250:260)=ones(11,11);
New_C=C;

Data=cell(1,2); % Build data cell containing time step and concentration values. 
Data{1}=0;
Data{2}=C;

wc=ones(2,1); % main weight cofficient wa, wb

% direction based cofficient: | south | north | east | west | SE | NW | NE | SW 
U=5; %m/s
theta=deg2rad(135); % degree, positive wind speed value point to center 
utemp=[-U*cos(theta),U*cos(theta),U*sin(theta),-U*sin(theta),U*sin(theta-pi/4),-U*sin(theta-pi/4),...
    U*cos(theta-pi/4),-U*cos(theta-pi/4)];
ud=(utemp+abs(utemp))/2;

kz=ones(2,1); % vertical turbulence cofficient kz', kz''
delta=zeros(I,J);   % deposition 
lambda=zeros(I,J);  % reaction

% Update rules
ii=2:I-1;
jj=2:J-1;


% Set control variables
pause=0;
stop=0;
run=0;
exit=0;

% Time step
Dt=0.01;
Dspace=0.13;
T=1;
t=0;
step=0;
while t<T
    t=t+Dt;
    step=step+1;
New_C(ii,jj)=C(ii,jj)+wc(1)*(ud(1)*(C(ii,jj-1)-C(ii,jj))+ud(2)*(C(ii,jj+1)-C(ii,jj))+ ...
         ud(3)*(C(ii+1,jj)-C(ii,jj))+ud(4)*(C(ii-1,jj)-C(ii,jj)))*Dt/Dspace + ...
    wc(2)*(ud(5)*(C(ii+1,jj-1)-C(ii,jj))+ud(6)*(C(ii-1,jj+1)-C(ii,jj))+ ...
         ud(7)*(C(ii+1,jj+1)-C(ii,jj))+ud(8)*(C(ii-1,jj-1)-C(ii,jj)))*Dt/Dspace/sqrt(2)+ ...
    wc(1)*(kz(1)*(C(ii,jj-1)-C(ii,jj))+kz(2)*(C(ii,jj+1)-C(ii,jj)))*Dt/Dspace - ...
    delta(ii,jj).*C(ii,jj)-lambda(ii,jj).*C(ii,jj);
% Results(step).C=New_C;
C=New_C;

hold on
drawnow 
[CTR,h]=contour(New_C');
delete(h);
title(['Time: ',num2str(t),' s']);
% pause(1);
end
