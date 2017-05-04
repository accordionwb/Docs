function fun_obst_cubic(cubic,general_conf)
% Ground centered cubic localization from cubic
% General configuration from generl_config
% including:
%        filename/fid
%        color
%        surf_id
fid=general_conf.fid;
color=general_conf.color;
surf_id=general_conf.surf_id;



M=size(cubic,1);
for i=1:M
    xmin=cubic(i,1);
    xmax=cubic(i,2);
    ymin=cubic(i,3);
    ymax=cubic(i,4);
    zmin=cubic(i,5);
    zmax=cubic(i,6);
    
    xc=mean(xmin,xmax);
    yc=mean(ymin,ymax);
    zc=mean(zmin,zmax);
    xl=xmax-xmin;
    yl=ymax-ymin;
    zl=zmax-zmin;
    
    
    
    info=['/ Cubic No.',num2str(i),', XYZ = ',num2str(xc),',',num2str(yc),',',num2str(zc),...
        ', LWH =: ',num2str(xl),',',num2str(yl),',',num2str(zl),' m.'];
    fprintf(fid,'%s\n',info);
    
    obst = ['&OBST XB=',num2str(xmin),',',num2str(xmax),',',num2str(ymin),',',num2str(ymax),...
        ',',num2str(zmin),',',num2str(zmax),',',...
        'COLOR=''',color,''',','SURF_ID=''',surf_id,'''/'];
    
    fprintf(fid,'%s\n',obst);
end
fprintf(fid,'\n');
