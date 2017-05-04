function fun_obst_sphere(Sphere_Data,general_conf)
% McDermott
% 6-16-08
% fds_cyl_obst.m
% Modified by A Bova
% 28Mar10

%%%% Input Defination %%%%%
% Sphere_Data(n,k): contains k sphere configuration of n groups
% Sphere_Data(n,1): OBST ID number, int
% Sphere_Data(k,2);  % Center Coordinate on X axis  (m)
% Sphere_Data(k,3);  % Center Coordinate on Y axis (m)
% Sphere_Data(k,4);  % Center Coordinate on Z axis (m)
% Sphere_Data(k,5);  % Mesh Cell Size, Assuming uniform grid resolution (m) 
% Sphere_Data(k,6);  % Outer Raduis (m)

%%% General Configuration %%%
% general_conf is a structure containing color, surf_id and fid and other
% parameter
% By default: 
if nargin==1
    general_conf.color='SILVER';
    general_conf.surf_id='INERT';
    general_conf.fid=fopen('obst.txt','w');
end


% constructs cylinder by rotating rectangular obstructions
color=general_conf.color;
surf_id=general_conf.surf_id;
fid=general_conf.fid;
N=size(Sphere_Data,1);

for k=1:N

index=Sphere_Data(k,1);  % OBST ID Number
xc=Sphere_Data(k,2);  % Center Coordinate on X axis  (m)
yc=Sphere_Data(k,3);  % Center Coordinate on Y axis (m)
zc=Sphere_Data(k,4);  % Center Coordinate on Z axis (m)
dx=Sphere_Data(k,5);  % Mesh Cell Size, Assuming uniform grid resolution (m) 
ro=Sphere_Data(k,6);  % Outer Raduis (m)
% r_shape=Sphere_Data(7);  % Sphere Shape 0-pi
%axis=char(Sphere_Data(9)); % Cylinder axis direction parallel to x 120, y 121 or z 122 (ASCII)


xo = ro/sqrt(3); 

% While loop will continue if cylinder wall is not at least one mesh cell thick



%l_cyl=input('Cylinder length (m): ');
% disp(['Cylinder length (m): ',num2str(l_cyl)]);


% disp(['Cylinder axis direction parallel to x, y or z?: ',axis,'s']);

x0 = [xc,yc,zc]';   % center coordinates (m)

% Define vectors from center to corner of obstruction (o_vec).
% and corner of hole (i_vec).

R_vec=[xo,xo,xo]';
  

% Change as necessary

% filename = ['cylinder_',num2str(index,'%02d'),'.obst'];


d_theta=asin(dx/ro); %minimum angle or rotation
d_phi=asin(dx/ro);
n_theta = round(pi/d_theta); % Turn through pi radians (half-circle).
n_phi=round(pi/d_phi);
% fid = fopen(filename,'w');

info=['/ Sphere No.',num2str(index),', XYZ = ',num2str(xc),',',num2str(yc),',',num2str(zc),...
    ', Radius: ',num2str(ro),' m.'];
fprintf(fid,'%s\n',info);

% fprintf('Writing to file: %s \n',filename);
% fprintf('Writing to file');
% fprintf('Progress.');
for theta = 0:n_theta
    for phi=0:n_phi
    
   
    % Define rotation matrix based on shpere axis.
 M=[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta);
     cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta);
     -sin(phi),cos(phi),0];

    
    % Rotate upper & lower corner coords
    
    R_vec = M*R_vec;
    
    lc_plot = x0 - R_vec; % Translate coords relative to center coords.
    uc_plot = x0 + R_vec;
    
    xmin=num2str(lc_plot(1));
    ymin=num2str(lc_plot(2));
    zmin=num2str(lc_plot(3));
    xmax=num2str(uc_plot(1));
    ymax=num2str(uc_plot(2));
    zmax=num2str(uc_plot(3));

    obst = ['&OBST XB=',xmin,',',xmax,',',ymin,',',ymax,',',zmin,',',zmax,',',...
            'COLOR=''',color,''',','SURF_ID=''',surf_id,'''/'];

    fprintf(fid,'%s\n',obst);
    
    end
end
disp('');
disp(['Sphere ',num2str(index,'%02d'),'/',num2str(N,'%02d'),' Complete!']);
fprintf(fid,'\n');
end
% fclose(fid);