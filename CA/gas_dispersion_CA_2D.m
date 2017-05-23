% Dimensions
I=60; %x
J=50; %y

% Initial and Source conditions
C=zeros(I,J); % Array initialization
S=zeros(I,J);  % Source Character
S(26,26)=10;
New_C=C;

Data=cell(100,2); % Build data cell containing time step and concentration values.
Data{1,1}=0;
Data{1,2}=C;

wc=ones(2,1); % main weight cofficient wa, wb

% direction based cofficient: | south | north | east | west | SE | NW | NE | SW
U=0; %m/s
theta=deg2rad(225); % degree, positive wind speed value point to center
utemp=[-U*cos(theta),U*cos(theta),U*sin(theta),-U*sin(theta),U*sin(theta-pi/4),-U*sin(theta-pi/4),...
    U*cos(theta-pi/4),-U*cos(theta-pi/4)];
ud=(utemp+abs(utemp))/2;

Kz=ones(2,1); % vertical turbulence cofficient kz', kz''
Kx=1e-4;
Ky=1e-4;
Kxy=1e-4;
delta=zeros(I,J);   % deposition
lambda=zeros(I,J);  % reaction

% Update rules inner index
ii=2:I-1;
jj=2:J-1;


% Set control variables
pause=0;
stop=0;
run=0;
exit=0;
dispersion=1;  % or 1

% Time and space step
Dt=0.01;
Dspace=0.01;
T=600;
t=0;
step=0;

% Plot options
v=[5,1,0.5,0.1,0.05];


% Numerical Stability Check
UFL0=0.35;
UFLU=U*Dt/Dspace;
UFLK=Kx*Dt/Dspace;
if max(UFLU,UFLK)>UFL0
    disp(['Current UFLU=',num2str(UFLU),', UFLK=',num2str(UFLK),' Limit is ',num2str(UFL0)]);
    DspaceU=U*Dt/UFLU;
    DspaceK=Kx*Dt/UFLK;
    error(['Dspace should be min(UFLU, UFLK): ',num2str(DspaceU),', ',num2str(DspaceK)])
else
    while t<T
        t=t+Dt;
        step=step+1;
        %%% Without Dispersion
        if dispersion == 0
            New_C(ii,jj)=C(ii,jj)+...
                wc(1)*Dt/Dspace*(ud(1)*(C(ii,jj-1)-C(ii,jj))+...
                ud(2)*(C(ii,jj+1)-C(ii,jj))+ud(3)*(C(ii+1,jj)-C(ii,jj))+...
                ud(4)*(C(ii-1,jj)-C(ii,jj))) + ...
                wc(2)*Dt/Dspace/sqrt(2)*(ud(5)*(C(ii+1,jj-1)-C(ii,jj))+...
                ud(6)*(C(ii-1,jj+1)-C(ii,jj))+ ...
                ud(7)*(C(ii+1,jj+1)-C(ii,jj))+...
                ud(8)*(C(ii-1,jj-1)-C(ii,jj)))+ ...
                wc(1)*Dt/Dspace*(Kz(1)*(C(ii,jj-1)-C(ii,jj))+...
                Kz(2)*(C(ii,jj+1)-C(ii,jj))) - ...
                delta(ii,jj).*C(ii,jj)-lambda(ii,jj).*C(ii,jj)+S(ii,jj)*Dt;
        end
        %%% With Dispersion
        if dispersion ==1
            %%% Main body
            New_C(ii,jj)=C(ii,jj)+...
                Kx*Dt/Dspace^2*(C(ii+1,jj)-2*C(ii,jj)+C(ii-1,jj))+...
                Ky*Dt/Dspace^2*(C(ii,jj+1)-2*C(ii,jj)+C(ii,jj-1))+...
                Kxy*Dt/2/Dspace^2*(C(ii-1,jj+1)-2*C(ii,jj)+C(ii+1,jj-1))+...
                Kxy*Dt/2/Dspace^2*(C(ii-1,jj-1)-2*C(ii,jj)+C(ii+1,jj+1))+...
                wc(1)*Dt/Dspace*(ud(1)*(C(ii,jj-1)-C(ii,jj))+...
                ud(2)*(C(ii,jj+1)-C(ii,jj))+ud(3)*(C(ii+1,jj)-C(ii,jj))+...
                ud(4)*(C(ii-1,jj)-C(ii,jj))) + ...
                wc(2)*Dt/Dspace/sqrt(2)*(ud(5)*(C(ii+1,jj-1)-C(ii,jj))+...
                ud(6)*(C(ii-1,jj+1)-C(ii,jj))+ud(7)*(C(ii+1,jj+1)-C(ii,jj))+...
                ud(8)*(C(ii-1,jj-1)-C(ii,jj)))- ...
                delta(ii,jj).*C(ii,jj)*Dt-lambda(ii,jj).*C(ii,jj)*Dt+S(ii,jj)*Dt;
            %%% Boundary Conditions
            % Edge I
            % Edge II
            % Edge III
            % Edge IV
            % Corner A
            % Corner B
            % Corner C
            % Corner D
        end
        % Results(step).C=New_C;
        C=New_C;
        
        
        %%% Plot Contour
        % hold on
        if rem(step,100)==0
            drawnow
            [CTR,obj] = contour(New_C',v);
            obj.LevelListMode='manual';
            obj.LevelStepMode='manual';
            obj.ShowText='on';
            title(['Time: ',num2str(t),' s']);
        end
        % pause(1);
    end
end
