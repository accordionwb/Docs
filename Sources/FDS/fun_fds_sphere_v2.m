% function fun_fds_Sphere_con(fid,sol,xc,yc,zc,dx,ro,l_cyl,axis,index,surf_id1,color)
function fun_fds_sphere_v2(xc,yc,zc,dx,ro)
% McDermott
% 6-16-08
% fds_cyl_obst.m
% Modified by A Bova
% 28Mar10

% constructs Sphere by rotating rectangular obstructions


%index=Sphere_Data(1);  % OBST ID Number
%xc=Sphere_Data(2);  % Center Coordinate on X axis  (m)
%yc=Sphere_Data(3);  % Center Coordinate on Y axis (m)
%zc=Sphere_Data(4);  % Center Coordinate on Z axis (m)
%dx=Sphere_Data(5);  % Mesh Cell Size, Assuming uniform grid resolution (m) 
%ro=Sphere_Data(6);  % Outer Raduis (m)


color='SILVER';
surf_id='INERT';

xo=ro/sqrt(3); 

% While loop will continue if Sphere wall is not at least one mesh cell thick

x0 = [xc,yc,zc]';   % center coordinates (m)
R_vec=[xo,xo,xo]';
  

% Change as necessary

filename = ['Sphere.obst'];


d_theta=asin(dx/ro); %minimum angle or rotation
d_phi=asin(dx/ro);
n_theta = round(pi/d_theta); % Turn through pi radians (half-circle).
n_phi=round(pi/d_phi);
fid = fopen(filename,'w');

info=['/ Obstruction Sphere Center, XYZ = ',num2str(xc),',',num2str(yc),',',num2str(zc),...
    ', Radius: ',num2str(ro),' m.'];
fprintf(fid,'%s\n',info);
zmax0=zc;
zmin0=zc;

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
    zmax0=max(zmax0,uc_plot(3));
    zmin0=min(zmin0,lc_plot(3));
    obst = ['&OBST XB=',xmin,',',xmax,',',ymin,',',ymax,',',zmin,',',zmax,',',...
            'COLOR=''',color,''',','SURF_ID=''',surf_id,'''/'];

    fprintf(fid,'%s\n',obst);
    
    end    
end
% disp('');
% disp(['Sphere ',num2str(index,'%02d'),'Complete!']);
fprintf(fid,'\n');
fclose(fid);
disp(['Zmax is :',num2str(zmax0),'\n'])
disp(['Zmin is :',num2str(zmin0),'\n'])