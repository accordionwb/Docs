function fun_fds_cylinder_v2(sol,xc,yc,zc,dx,ro,l_cyl,axis,index)
% McDermott
% 6-16-08
% fds_cyl_obst.m
% Modified by A Bova
% 28Mar10

% constructs cylinder by rotating rectangular obstructions


disp(['(1) Solid, (2) hollow or (3) "hole only" cylinder?: ',num2str(sol)] );
disp(['Input center coordinates of cylinder using your domain coordinate system.']);
disp(['Center x-coord: ',num2str(xc)]);
disp(['Center y-coord: ',num2str(yc)]);
disp(['Center z-coord: ',num2str(zc)]);
disp(['Mesh cell size (m): ',num2str(dx)]);    % ***** Assumes uniform grid resolution (m)****
disp(['Outer Radius (m): ',num2str(ro)]);
xo = sqrt(2)*ro;
ri=ro;
xi=xo;

% While loop will continue if cylinder wall is not at least one mesh cell thick

if (sol==2)
    while ((ro-ri) < dx);
        ri = input('Inner radius (m): ');
        xi = sqrt(2)*ri;
    end;
end;

%l_cyl=input('Cylinder length (m): ');
disp(['Cylinder length (m): ',num2str(l_cyl)]);

l_hol=l_cyl;

if (sol==2)
    disp('Make hole longer than cylinder to ensure that it penetrates,');
    disp('make shorter than cylinder if only a hollow inner chamber is desired.');
    l_hol = input('Hole length: ');
end;

disp(['Cylinder axis direction parallel to x, y or z?: ',axis,'s']);

x0 = [xc,yc,zc]';   % center coordinates (m)

% Define vectors from center to corner of obstruction (o_vec).
% and corner of hole (i_vec).

if (axis=='x')
    o_vec=0.5*[l_cyl,xo,xo]';
    i_vec=0.5*[l_hol,xi,xi]';
else
    if (axis=='y')
        o_vec=0.5*[xo,l_cyl,xo]';
        i_vec=0.5*[xi,l_hol,xi]';
    else o_vec=0.5*[xo,xo,l_cyl]';
        i_vec=0.5*[xi,xi,l_hol]';
    end;
    
end;

% Change as necessary
surf_id1='''INERT''';
color ='''SLATE BLUE''';
filename = ['cylinder_',num2str(index),'.obst'];


da=asin(dx/ro); %minimum angle or rotation
rmax = round(pi/da); % Turn through pi radians (half-circle).

fid = fopen(filename,'w');

info=['/ Obstruction No.',num2str(index),', Diameter: ',num2str(2*ro),' m.'];
fprintf(fid,'%s\r',info);

fprintf('Writing to file: %s \r',filename);
fprintf('Progress.');
cntr=0;
for theta = 0:rmax
    
    %progress bar
    if (cntr==10);
        fprintf('. \b');
        cntr=1;
    end;
    
    % Define rotation matrix based on cylinder axis.
    
    if (axis=='x')
        M=[1., 0, 0; 0, cos(theta),-sin(theta); 0, sin(theta), cos(theta)];
    else
        if (axis=='y')
            M=[cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
        else M=[cos(theta),-sin(theta),0; sin(theta), cos(theta),0; 0, 0, 1];
        end;
    end;
    
    % Rotate upper & lower corner coords
    
    o_vec = M*o_vec;
    
    lc_plot = x0 - o_vec; % Translate coords relative to center coords.
    uc_plot = x0 + o_vec;
    
    xmin=num2str(lc_plot(1));
    ymin=num2str(lc_plot(2));
    zmin=num2str(lc_plot(3));
    xmax=num2str(uc_plot(1));
    ymax=num2str(uc_plot(2));
    zmax=num2str(uc_plot(3));
    
    if (sol==3)
        obst=['&HOLE XB=',xmin,',',xmax,',',ymin,',',ymax,',',zmin,',',zmax,'/'];
    else
        obst = ['&OBST XB=',xmin,',',xmax,',',ymin,',',ymax,',',zmin,',',zmax,',',...
            'COLOR=',color,',','SURF_ID=',surf_id1,'/'];
    end;
    
    fprintf(fid,'%s\r',obst);
    
    % Do hollow cylinder calculations only if necessary
    
    if (sol==2)
        i_vec = M*i_vec;
        
        lh_plot = x0 - i_vec; % Translate coords relative to center coords.
        uh_plot = x0 + i_vec;
        
        xhmin=num2str(lh_plot(1));
        yhmin=num2str(lh_plot(2));
        zhmin=num2str(lh_plot(3));
        xhmax=num2str(uh_plot(1));
        yhmax=num2str(uh_plot(2));
        zhmax=num2str(uh_plot(3));
        hole = ['&HOLE XB=',xhmin,',',xhmax,',',yhmin,',',yhmax,',',zhmin,',',zhmax,' /'];
        
        fprintf(fid,'%s\r',hole);
    end;
    
    cntr=cntr+1;
end
disp('');
disp('Complete');
fclose(fid);