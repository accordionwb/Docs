% function fun_fds_cylinder_con(fid,sol,xc,yc,zc,dx,ro,l_cyl,axis,index,surf_id1,color)
%function fun_fds_cylinder_con(fid,Cylinder_Data,color,surf_id)
function fun_obst_cylinder(Cylinder_Data,general_conf)
% McDermott
% 6-16-08
% fds_cyl_obst.m
% Modified by A Bova
% 28Mar10

% constructs cylinder by rotating rectangular obstructions

color=general_conf.color;
surf_id=general_conf.surf_id;
fid=general_conf.fid;
N=size(Cylinder_Data,1);

cntr=0;  % Progress

for k=1:N  % Begin Data Set
    index=Cylinder_Data(k,1);  % OBST ID Number
    sol=Cylinder_Data(k,2);  % (1) Solid, (2) hollow or (3) "hole only"
    xc=Cylinder_Data(k,3);  % Center Coordinate on X axis  (m)
    yc=Cylinder_Data(k,4);  % Center Coordinate on Y axis (m)
    zc=Cylinder_Data(k,5);  % Center Coordinate on Z axis (m)
    dx=Cylinder_Data(k,6);  % Mesh Cell Size, Assuming uniform grid resolution (m)
    ro=Cylinder_Data(k,7);  % Outer Raduis (m)
    l_cyl=0.5*Cylinder_Data(k,8);  % Cylinder length (m)
    axis=char(Cylinder_Data(k,9)); % Cylinder axis direction parallel to x 120, y 121 or z 122 (ASCII)
    
    
    xo = ro/sqrt(2);
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
    % disp(['Cylinder length (m): ',num2str(l_cyl)]);
    
    l_hol=l_cyl;
    
    if (sol==2)
        % disp('Make hole longer than cylinder to ensure that it penetrates,');
        % disp('make shorter than cylinder if only a hollow inner chamber is desired.');
        l_hol = input('Hole length: ');
    end;
    
    % disp(['Cylinder axis direction parallel to x, y or z?: ',axis,'s']);
    
    x0 = [xc,yc,zc]';   % center coordinates (m)
    
    % Define vectors from center to corner of obstruction (o_vec).
    % and corner of hole (i_vec).
    
    if (axis=='x')
        o_vec=[l_cyl,xo,xo]';
        i_vec=[l_hol,xi,xi]';
    else
        if (axis=='y')
            o_vec=[xo,l_cyl,xo]';
            i_vec=[xi,l_hol,xi]';
        else o_vec=[xo,xo,l_cyl]';
            i_vec=[xi,xi,l_hol]';
        end;
        
    end;
    
    % Change as necessary
    
    % filename = ['cylinder_',num2str(index,'%02d'),'.obst'];
    
    
    da=asin(dx/ro); %minimum angle or rotation
    rmax = round(pi/da); % Turn through pi radians (half-circle).
    
    % fid = fopen(filename,'w');
    
    info=['/ Cylinder No.',num2str(index),', XYH = ',num2str(xc),',',num2str(yc),',',num2str(l_cyl),...
        ', Diameter: ',num2str(2*ro),' m.'];
    fprintf(fid,'%s\n',info);
    
    % fprintf('Writing to file: %s \n',filename);
    % fprintf('Writing to file');
    % fprintf('Progress.');
    %progress bar
    if (cntr==10);
                fprintf('. \b');
        cntr=1;
    end;
    
    for theta = 0:rmax
        
        
        
        % Define rotation matrix based on cylinder axis.
        
        if (axis=='x')
            M=[1, 0, 0; 0, cos(theta),-sin(theta); 0, sin(theta), cos(theta)];
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
                'COLOR=''',color,''',','SURF_ID=''',surf_id,'''/'];
        end;
        
        fprintf(fid,'%s\n',obst);
        
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
            
            fprintf(fid,'%s\n',hole);
        end;
        
        cntr=cntr+1;
    end
    disp('');
    disp(['Cylinder ',num2str(index,'%02d'),'/',num2str(N,'%02d'),' Complete!']);
    fprintf(fid,'\n');
end
% fclose(fid);