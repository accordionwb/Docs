%% Parameters 
parameter.IJ=[700,500];  %x,y,z Demisions
I=parameter.IJ(1);
J=parameter.IJ(2);
% wind direction based cofficient: | south | north | east | west | SE | NW | NE | SW
parameter.U=0.5; %m/s wind speed
parameter.theta=deg2rad(225); % degree, positive wind speed value point to center
% Dispersion coefficient
parameter.Kz=ones(2,1); % vertical turbulence coefficient kz', kz''
parameter.Kx=1e-3;  % horizental dispersion coefficient on x axis
parameter.Ky=1e-3;
parameter.Kxy=1e-3;
% Reaction and deposition
parameter.delta=zeros(I,J);   % deposition
parameter.lambda=zeros(I,J);  % reaction
% Adjoint coefficient
parameter.wc=ones(2,1); % main weight cofficient wa, wb
% Time and space division
parameter.Dt=0.01;  % time step
parameter.Dspace=0.1;  % space step
% control variables
parameter.dispersion=1;  % or 1
% Initial and Source conditions
C=zeros(I,J); % Array initialization
S=zeros(I,J);  % Source Character
S(355:360,255:260)=10*ones(6,6);
parameter.Source=S; %source 

%% Loop settings
% Time and space step
T=10; % time end. (start from t=0)
t=0;
step=0;

% Plot options
v=[5,1,0.5,0.1,0.05];  % contour level


% Numerical Stability Check
UFL0=0.35;
UFLU=parameter.U*parameter.Dt/parameter.Dspace;
UFLK=parameter.Kx*parameter.Dt/parameter.Dspace;
if max(UFLU,UFLK)>UFL0
    disp(['Current UFLU=',num2str(UFLU),', UFLK=',num2str(UFLK),' Limit is ',num2str(UFL0)]);
    DspaceU=parameter.U*parameter.Dt/UFL0;
    DspaceK=parameter.Kx*parameter.Dt/UFL0;
    error(['Dspace should be max(UFLU, UFLK): ',num2str(DspaceU),', ',num2str(DspaceK)])
else
    while t<T
        t=t+parameter.Dt;
        step=step+1;
        New_C=fun_rule2D(C,parameter);
        %%% Without Dispersion
        
        % Results(step).C=New_C;
        C=New_C;
        
        
        %%% Plot Contour
        
        if rem(step,100)==0
            drawnow

            im = image(C','CDataMapping','scaled');
%             colorbar
            colormap(jet);
            ax = gca;
            ax.YDir='normal';
            
            hold on
            
            [CTR,obj] = contour(C',v);
            obj.ShowText='on';            
            title(['Time: ',num2str(t),' s']);
            
            hold off
        end
        % pause(1);
    end
end
