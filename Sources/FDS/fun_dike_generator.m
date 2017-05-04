function fun_dike_generator(conf)
%% This function generate one single domain file with specific file name no return 
% The format of the following parameter is pre-specified. 
%
% Usage: fun_domain_generator(conf) where configure is a structure containing everything needed
%
%    STRUCTURE        TYPE    EXAMPLE
% conf. 
%  |-- head.
%  |    |-- chid       str   'Domain_001'
%  |    |-- title      str   'Bing FDS simulation on domain xxx '
%  |
%  |-- mesh.  
%  |-- misc. 
%  |    |-- restart    str  '.TRUE.' | '.FALSE.'
%  |
%  |-- time            flt   
%  |-- specid          str  'METHANE' | 'CHLORINE'
%  |-- rlse.
%  |    |-- surfid     str  'LEAK'
%  |    |-- mass       flt  2.0 (kg/s)
%  |    |-- rampid     str  'leak_ramp'
%  |    |-- ramp(:,II)  flt  T=%1   F=%2
%  |    |-- ventXB     flt   6 coordinates
%  |    |-- color      str  'RED' 
%  |    
%  |-- wind.
%  |    |-- surfid()   cell  {'Wind_X','Wind_Y','OPEN'}
%  |    |-- temp       flt   15.0 
%  |    |-- VX(III)    flt  VEL=%1   VEL_T=%2, %3
%  |    |-- VY(III)    flt  VEL=%1   VEL_T=%2, %3 
%  |    |-- profile    str  'ATMOSPHERIC'
%  |    |-- Z0         flt   Z0=1.0
%  |    |-- rampid     str   'WindRamp'
%  |    |-- ramp(:,II) flt  T=%1, F=%2
%  |    |-- vent(V)    int  1:WIND_X | 2:WIND_Y | 3:OPEN  [1 3 2 3 3];
%  |   
%  |-- dump.
%  |    |-- massfile   str   '.TRUE.'
%  |
%  |-- slif().
%  |    |-- surface     str    'PBX' | 'PBY' | 'PBZ'
%  |    |-- PB          flt    1.0 | 2.0 | ...
%  |    |-- quantity    str    'VOLUME FRACTION'
%  |
%  |-- devc.
%  |    |-- xyz_config  str    'MAT_Devc_config.mat'
%  |    |-- quantity    str    'VOLUME FRACTION' | 'MASS FRACTION'
%  |
%  |-- obst.
%       |-- rec_config  str    'MAT_Rec_config.mat'
%       |-- rec_color   str    'OLIVE DRAB'
%       |-- rec_surfid  str    'INERT'
%       |-- cyl_config  str    'MAT_Cylinder_config.mat'
%       |-- cyl_color   str    'SILVER'
%       |-- cyl_surfid  str    'INERT'

%% Read configuration
if exist(['Source/',conf.head.chid],'dir')
    disp(['Source/',conf.head.chid,' already exist']);
else
    mkdir('Source/',conf.head.chid);
end
filename=['Source/',conf.head.chid,'/',conf.head.chid,'.fds'];
fid=fopen(filename,'w');


%% Writing configure into FDS script

% 1st entry &HEAD
fprintf(fid,'&HEAD CHID = ''%s'', TITLE=''%s'' /\n',conf.head.chid, conf.head.title);
fprintf(fid,'\n');

% 2nd entry &MESH
fprintf(fid,'%s\n','// Start mesh configuration');
% meshXB(1).line='&MESH ID=''mesh_CoreX'',  IJK=100, 100, 30, XB= -50.0, 50.0,  -50.0, 50.0, 0.0, 30.0 /';  % Zone 1
% meshXB(2).line='&MESH ID=''mesh_CoreY'',  IJK=100, 100, 30, XB=  -50.0, 50.0,    -50.0, 50.0, 0.0, 30.0 /';  % Zone 2
% meshXB(3).line='&MESH ID=''mesh_ExtendX'',  IJK=50, 50, 30, XB= 50.0, 150,   -50.0, 50.0, 0.0, 30.0 /';  % Zone 3
% meshXB(4).line='&MESH ID=''mesh_ExtendY'',  IJK=50, 50, 30, XB=   -50.0, 50.0, 50.0, 150.0, 0.0, 30.0 /';  % Zone 4
% meshXB(5).line='&MESH ID=''mesh_center'', IJK=80, 100,  60,  XB = -20, 20, -25, 25, 0, 15 /';  % Zone 5
% meshXB(6).line='&MESH ID=''mesh_tank'', IJK=96, 192,  80,  XB = -6, 6, -12, 12, 0, 10 /';  % Zone 6

meshXB(1).line='/MESH ID=''mesh_CoreX'',  IJK=100, 100, 30, XB= -50.0, 50.0,  -50.0, 50.0, 0.0, 30.0 /';  % Zone 1
meshXB(2).line='/MESH ID=''mesh_CoreY'',  IJK=100, 100, 30, XB=  -50.0, 50.0,    -50.0, 50.0, 0.0, 30.0 /';  % Zone 2
meshXB(3).line='/MESH ID=''mesh_ExtendX'',  IJK=50, 50, 30, XB= 50.0, 150,   -50.0, 50.0, 0.0, 30.0 /';  % Zone 3
meshXB(4).line='/MESH ID=''mesh_ExtendY'',  IJK=50, 50, 30, XB=   -50.0, 50.0, 50.0, 150.0, 0.0, 30.0 /';  % Zone 4
meshXB(5).line='&MESH ID=''mesh_center'', IJK=150, 200,  100,  XB = -15, 15, -20, 20, 0, 20 /';  % Zone 5
meshXB(6).line='&MESH ID=''mesh_tank'', IJK=100, 200,  100,  XB = -5, 5, -10, 10, 0, 10 /';  % Zone 6

for m=1:numel(conf.mesh)
    if conf.mesh(m)~=0
        fprintf(fid, '%s\n',meshXB(m).line);
    end
end
fprintf(fid,'\n');

% 3rd entry & MISC
fprintf(fid,'&MISC RESTART= %s / # It could be true \n',conf.misc.restart);
fprintf(fid,'&MISC MEAN_FORCING(1:2)=.TRUE.,.TRUE., U0=%2.2f, V0=%2.2f, DT_MEAN_FORCING=0.1 / \n',...
    conf.misc.wind);
fprintf(fid,'\n');

% 4th entry &TIME
fprintf(fid,'&TIME T_END=%3.1f, /\n',conf.time);
fprintf(fid,'\n');

% 5th entry &SPEC
fprintf(fid,'&SPEC ID=''%s'' /\n',conf.specid);
fprintf(fid,'\n');

% 6th entry rlse &SURF
fprintf(fid,'&SURF ID=''%s'', SPEC_ID=''%s'', MASS_FLUX(1)=%2.4f, RAMP_MF(1)=''%s'', TMP_BACK=-161.55, TMP_INNER=-161.55, TMP_FRONT=-161.55 /\n', ...
    conf.rlse.surfid, conf.specid, conf.rlse.mass, conf.rlse.rampid);

[M,N]=size(conf.rlse.ramp);
for i=1:M;
    fprintf(fid,'&RAMP ID=''%s'', T=%3.1f, F=%1.3f, /\n',conf.rlse.rampid, conf.rlse.ramp(i,:));
end
fprintf(fid,'\n');

fprintf(fid,'&VENT XB= %4.1f, %4.1f, %4.1f, %4.1f, %4.1f, %4.1f, SURF_ID=''%s'', COLOR=''%s'' /\n',...
    conf.rlse.ventXB,conf.rlse.surfid,conf.rlse.color);
fprintf(fid,'\n');


axis_cfg={'XMIN','XMAX','YMIN','YMAX','ZMAX'};
for i=1:5   %conf.wind.vent is integer numbered  1|2|3 
    fprintf(fid,'&VENT MB=''%s'', SURF_ID=''%s'' /\n',char(axis_cfg(i)),'OPEN');
end
fprintf(fid,'\n');

% 7th entry &DUMP
fprintf(fid,'&DUMP MASS_FILE=%s  /\n',conf.dump.massfile);
fprintf(fid,'&DUMP PLOT3D_SPEC_ID(5)=''%s'', PLOT3D_QUANTITY(1:5)=''TEMPERATURE'',''U-VELOCITY'',''V-VELOCITY'',''W-VELOCITY'',''%s'', DT_PL3D=%3.1f /\n', ...
    conf.specid,conf.dump.pl3d_quantity,conf.dump.dt_pl3d);
fprintf(fid,'\n');

% 8th entry &ISOF
fprintf(fid,'&ISOF SPEC_ID=''%s'', QUANTITY=''%s'', VALUE(1)=%2.2e, VALUE(2)=%2.2e, VALUE(3)=%2.2e   /\n',...
    conf.specid, conf.isof.quantity, conf.isof.value);
fprintf(fid,'\n');


% 9th entry &SLCF 
fprintf(fid,'%s\n','&SLCF PBZ=1.0, QUANTITY=''VELOCITY'', VECTOR=.TRUE. /');
fprintf(fid,'\n');
for i=1:numel(conf.slif)
    fprintf(fid,'&SLCF %s =%3.1f, QUANTITY=''%s'', SPEC_ID=''%s'' /\n',conf.slif(i).surface, conf.slif(i).PB,...
        conf.slif(i).quantity, conf.specid);
end
fprintf(fid,'\n');

% 10th entry: &DEVC
devc=conf.devc.xyz_config;
for i=1:numel(devc)
    if conf.mesh(devc(i).zone) ==1
%        fprintf(fid,'&DEVC ID=''%s'',XYZ=%3.1f,%3.1f,%3.1f, QUANTITY=''%s'', SPEC_ID=''%s'' /\n',...
 %           devc(i).id, devc(i).xyz, conf.devc.quantity, conf.specid);
    end
end
fprintf(fid,'\n');

% 11th entry: &OBST

% Writing Rectangle Data
Rec_Data=conf.obst.rec_config;
% surf_id1='''INERT''';
% color1 ='''OLIVE DRAB''';
fprintf(fid,'%s\n','/Starting writing rectangular OBST');
[N1,M1]=size(Rec_Data);
for i=1:N1
    fprintf(fid,'&OBST XB=%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f, COLOR=''%s'', SURF_ID=''%s''  /\n',...
        Rec_Data(i,:),conf.obst.rec_color, conf.obst.rec_surfid);
end
fprintf(fid,'\n');

% Writing Cylinder Data
Cylinder_Data=conf.obst.cyl_config;
[N2,M2]=size(Cylinder_Data);
for i=1:N2
    fun_fds_cylinder_con(fid,Cylinder_Data(i,:),conf.obst.cyl_color,conf.obst.cyl_surfid)  
end
fprintf(fid,'\n');

% Last entry: &TAIL
fprintf(fid,'&TAIL /\n');

fclose(fid);
disp(['FDS Script ',conf.head.chid,'.fds has been generated']);
end