%% Lattice Boltzmann Free-boundary Air flow model
% Developed by BIng Wang, ECUST, Shanghai, China
% Basic model contains 2-D D2Q9 LBM-SRT Configuration
% Side View, Building flow
% D2Q9 configuration is illustrated as follows:
% 6   2   5
%   \ | /
% 3 - 9 - 1
%   / | \
% 7   4   8

%% Define Physical Parameters
clear
clc

L = 30; % length m
H = 20; % length m
h = 5; % Obst height m, Characteristic length.
l = 1;  % Obst width, m
Tf = 8; % Obst SW corner, m
U0 = 5; % Mean wind velocity, m/s
rho0 = 1; % mean density, kg/m3
mu0 = 1/600;  % Mean kinematic viscosity
% Re = U0*h/mu0;
Re = 200;
obst = [8,8+l,0,h]; % obstruction xx and yy

% Main program
% Lattice Unity
% Reynold number should be retained
% Lattice velocity and lattice viscosity could be determined arbitrarily.
Uf = 0.1; % Lattice fluid velocity
Ma = Uf * sqrt(3)
if Ma > 0.3
    error('Mach number is too large, try smaller U0');
end
Nh = 50; % 100 lattice on characteristic length, h
mstep = 10000; % Maximum iteration
% Calculated values
mu = Uf*Nh/Re;
NX = Nh*L/h;
NY = Nh*H/h+1;
xx = 1:NX;
yy = 2:NY-1;
omega = 1 / (3*mu+0.5);
iOBST = [obst(1)/L*NX,obst(2)/L*NX,obst(3)/H*NY,obst(4)/H*NY];
tPlot=10;

[Y,X] = meshgrid(1:NY,1:NX);
OBST = X>=iOBST(1) & X<=iOBST(2) & Y>=iOBST(3) & Y<=iOBST(4);
% OBST = zeros(NX,NY);
OBST(:,[1,NY]) = 1;

bbRegion=find(OBST);

% D2Q9 LATTICE CONSTANTS
w  = [ 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36, 4/9];
ex = [ 1,  0, -1,  0,    1,  -1,  -1,   1, 0];
ey = [ 0,  1,  0, -1,    1,   1,  -1,  -1, 0];
origin=[1   2   3   4   5   6   7   8   9];  %
opp = [ 3,   4,  1,  2,  7,  8,  5,  6,   9];

% Initial condition (rho = 0)
fIn = reshape( w' * ones(1,NX*NY), 9, NX, NY);

% Boundary type control
type = 1    % 0 -> unknown |  1 -> regular

%==========================================================================
%% Main Loop
it=1;
while it <= mstep
    
    % MACROSCOPIC VARIABLES
    rho = sum(fIn);
    UX  = reshape ( ...
        (ex * reshape(fIn,9,NX*NY)), 1,NX,NY) ./rho;
    UY  = reshape ( ...
        (ey * reshape(fIn,9,NX*NY)), 1,NX,NY) ./rho;
    
    % MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS
    % Inlet: Poiseuille profile
    %     L = NY-2;
    %     UX(:,1,yy) = 4 * Uf / (L*L) * (yy.*L-yy.*yy);;
    UX(:,1,yy) = Uf;
    UY(:,1,yy) = 0;
    rho(:,1,yy) = 1 ./ (1-UX(:,1,yy)) .* ( ...
        sum(fIn([2,9,4],1,yy)) + ...
        2*sum(fIn([3,6,7],1,yy)) );
    
    % Outlet: Second order extrapolation
%         rho(:,NX,yy) = rho(:,NX-1,yy);
%         UY(:,NX,yy)  = 0;
%         UX(:,NX,yy)  = UX(:,NX-1,yy);
    
    
    % COLLISION STEP
    for i=1:9
        cu = 3*(ex(i)*UX+ey(i)*UY);
        fEq(i,:,:)  = rho .* w(i) .* ...
            ( 1 + cu + 1/2*(cu.*cu) ...
            - 3/2*(UX.^2+UY.^2) );
        fOut(i,:,:) = omega .* fEq(i,:,:) + ...
            ( 1 - omega ) .* fIn(i,:,:);
    end
    
    
    % MICROSCOPIC BOUNDARY CONDITIONS
    
    if type ==0
        for i=1:9
            % Left boundary (West)
            fOut(i,1,yy) = fEq(i,1,yy) + ...
                18*w(i)*ex(i)*ey(i)* ( fIn(7,1,yy) - ...
                fIn(6,1,yy)-fEq(7,1,yy)+fEq(6,1,yy) );
            
            % Right boundary
            fOut(i,NX,yy) = fEq(i,NX,yy) + ...
                18*w(i)*ex(i)*ey(i)* ( fIn(5,NX,yy) - ...
                fIn(8,NX,yy)-fEq(5,NX,yy)+fEq(8,NX,yy) );
            
            
            % Bounce back region
            fOut(i,bbRegion) = fIn(opp(i),bbRegion);  % opp as Rotation
        end
    end
    
    if type ==1
        for i=1:9
            fOut(i,bbRegion) = fIn(opp(i),bbRegion);  % opp as Rotation
        end
    end
    
    
    % STREAMING STEP
    for i=1:9
        fIn(i,:,:) = circshift(fOut(i,:,:), [0,ex(i),ey(i)]);
    end
    
    if type == 1
        % Left boundary, 2nd Order
        fOut(1,1,yy)=fIn(3,1,yy)+2*rho(:,1,yy).*Uf/3;
        fOut(5,1,yy)=fIn(7,1,yy)+rho(:,1,yy).*Uf/6;
        fOut(8,1,yy)=fIn(6,1,yy)+rho(:,1,yy).*Uf/6;
        
        % Right boundary, 2nd Order
        fOut(7,NX,yy)=fIn(7,NX-1,yy);%-fIn(7,NX-2,yy);
        fOut(3,NX,yy)=fIn(3,NX-1,yy);%-fIn(3,NX-2,yy);
        fOut(6,NX,yy)=fIn(6,NX-1,yy);%-fIn(6,NX-2,yy);
    end
    
    
    % VISUALIZATION
    if (mod(it,tPlot)==0)
        U = reshape(sqrt(UX.^2+UY.^2),NX,NY);
        U(bbRegion) = nan;
        imagesc(U');
        ax = gca;
        ax.YDir='normal';
        title(['iteration ',num2str(it),' / ',num2str(mstep)])
                rectangle('Position',[iOBST(1),iOBST(3),iOBST(2)-iOBST(1),...
                    iOBST(4)-iOBST(3)],'FaceColor','cyan')
        axis equal off; drawnow
    end
    it=it+1;
end

%% Velocity field analysis
iX=floor(linspace(1,NX,3));
UX0=reshape(UX,NX,NY);
UY0=reshape(UY,NX,NY);

for k=1:length(iX)
    plot(UX0(iX(k),:))
    hold on
end
hold off

