%% Lattice Boltzmann Free-boundary Air flow model
% Developed by BIng Wang, ECUST, Shanghai, China
% Basic model contains 2-D D2Q9 LBM-SRT Configuration
% Side View, Building flow 

% Define Physical Parameters
L = 80; % length m
H = 40; % length m
h = 10; % Obst height m, Characteristic length. 
l = 1;  % Obst width, m
Tf = 8; % Obst SW corner, m
U0 = 10; % Mean wind velocity, m/s
rho0 = 1; % mean density, kg/m3
mu0 = 1/600;  % Mean kinematic viscosity
Re = U0*h/mu0;
obst = [8,8+l,0,h]; % obstruction xx and yy 

% Main program
% Lattice Unity
% Reynold number should be retained
% Lattice velocity and lattice viscosity could be determined arbitrarily. 
Uf = 0.1; % Lattice fluid velocity 
Nh = 100; % 100 lattice on characteristic length, h
mstep = 10000; % Maximum iteration
% Calculated values
mu = Uf*Nh/Re;
NX = Nh*L/h;
NY = Nh*H/h+1;
xx = 1:NX;
yy = 2:NY;
omega = 1 / (3*mu+0.5);
iOBST = [obst(1)/L*NX,obst(2)/L*NX,obst(3)/H*NY,obst(4)/H*NY];

[X,Y] = meshgrid(1:NX,1:NY);
OBST = X>=iOBST(1) & X<=iOBST(2) & Y>=iOBST(3) & Y<=iOBST(4);
OBST(:,1) = 1;
bbregion=find(OBST);

% D2Q9 LATTICE CONSTANTS
w  = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
ex = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
ey = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
opp = [ 1,   4,  5,  2,  3,    8,   9,   6,   7];

% Initial condition (rho = 0)
fIn = reshape( w' * ones(1,NX*NY), 9, NX, NY);
rho = sum(fIn);
UX  = reshape ( ...
          (ex * reshape(fIn,9,NX*NY)), 1,NX,NY) ./rho;
UY  = reshape ( ...
          (ey * reshape(fIn,9,NX*NY)), 1,NX,NY) ./rho;
      
%==========================================================================
% Main Loop
it=1;
while it <= mstep
    
    % MACROSCOPIC VARIABLES
    rho = sum(fIn);
UX  = reshape ( ...
          (ex * reshape(fIn,9,NX*NY)), 1,NX,NY) ./rho;
UY  = reshape ( ...
          (ey * reshape(fIn,9,NX*NY)), 1,NX,NY) ./rho;
    
    % BOUNDARY CONDITIONS
    % West Inlet: Plug flow
    UX(:,1,yy) = Uf;
    UY(:,1,yy) = 0;
    rhow = 1 ./ (1-UX(:,1,yy)) .* ( ...
        sum(fIn([1,3,5],1,yy)) + ...
        2*sum(fIn([4,7,8],1,yy)) );
    rho(:,1,yy)=rhow;
    % East Outlet: Open boundary
    rho(:,NX,yy) = rho(:,NX-1,yy);
    UY(:,NX,yy)  = 0;
    UX(:,NX,yy)  = UX(:,NX-1,yy);
    % North boundary: Open boundary
    % South Boundary: bounce back
    fIn(1,NX,yy)=2.*fIn(1,NX-1,yy)-fIn(1,NX-2,yy);
    
    
    % COLLISION STEP
    for i=1:9
       cu = 3*(ex(i)*UX+ey(i)*UY);
       fEq(i,:,:)  = rho .* w(i) .* ...
         ( 1 + cu + 1/2*(cu.*cu) ...
           - 3/2*(UX.^2+UY.^2) );
       fOut(i,:,:) = omega .* fEq(i,:,:) + ...
         ( 1 - omega ) .* fIn(i,:,:);
    end
    
