function fun_spectra_generator(conf)
%% This function generate one single spectra analysis file with specific file name no return
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
%  |-- mesh.      [0-9], indicating mesh part
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
%  |-- slcf().
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
if exist(['script/',conf.head.chid],'dir')
    disp(['script/',conf.head.chid,' already exist']);
else
    mkdir('script/',conf.head.chid);
end
filename=['script/',conf.head.chid,'/',conf.head.chid,'.fds'];
fid=fopen(filename,'w');
dike_flag=conf.obst.dike_flag;


%% Writing configure into FDS script

%% 1st entry &HEAD
fprintf(fid,'&HEAD CHID = ''%s'', TITLE=''%s'' /\n',conf.head.chid, conf.head.title);
fprintf(fid,'\n');

%% 2nd entry &MESH
fprintf(fid,'%s\n','// Start mesh configuration');
% meshXB(9).line='&MESH ID=''mesh00'',  IJK=215, 60, 50, XB=   0.0, 215.0,    -0.0, 60.0, 0.0, 50.0 /';  % core area
% meshXB(1).line='&MESH ID=''mesh01'',  IJK=50, 50, 50, XB= -50.0, 0.0,  -50.0, 0.0, 0.0, 50.0 /';
% meshXB(2).line='&MESH ID=''mesh02'',  IJK=50, 60, 50, XB=  -50.0, 0.0,    0.0, 60.0, 0.0, 50.0 /';
% meshXB(3).line='&MESH ID=''mesh03'',  IJK=50, 50, 50, XB= -50.0, 0.0,   60.0, 110.0, 0.0, 50.0 /';
% meshXB(4).line='&MESH ID=''mesh04'',  IJK=215, 50, 50, XB=   0.0, 215.0, 60.0, 110.0, 0.0, 50.0 /';
% meshXB(5).line='&MESH ID=''mesh05'',  IJK=50, 50, 50, XB=   215.0, 265.0, 60.0, 110.0, 0.0, 50.0 /';
% meshXB(6).line='&MESH ID=''mesh06'',  IJK=50, 60, 50, XB=  215.0, 265.0, 0.0, 60.0, 0.0, 50.0 /';
% meshXB(7).line='&MESH ID=''mesh07'',  IJK=50, 50, 50, XB= 215.0, 265.0, -50.0, 0.0, 0.0, 50.0 /';
% meshXB(8).line='&MESH ID=''mesh08'',  IJK=215, 50, 50, XB=   0.0, 215.0, -50.0, 0.0, 0.0, 50.0 /';
% for m=1:numel(conf.mesh)
%     if conf.mesh(m)~=0
%         fprintf(fid, '%s\n',meshXB(m).line);
%     end
% end

meshXB(1).line='&MESH ID=''mesh01'',  IJK=60, 180, 50, XB= -60.0, 0.0,  -60.0, 120.0, 0.0, 50.0 /';
meshXB(2).line='&MESH ID=''mesh03'',  IJK=215, 60, 50, XB= 0.0, 215.0,   -60.0, 0.0, 0.0, 50.0 /';
meshXB(3).line='&MESH ID=''mesh00'',  IJK=215, 60, 50, XB=   0.0, 215.0,  0.0, 60.0, 0.0, 50.0 /';
meshXB(4).line='&MESH ID=''mesh02'',  IJK=215, 60, 50, XB=  0.0, 215.0,  60.0, 120.0, 0.0, 50.0 /';
meshXB(5).line='&MESH ID=''mesh04'',  IJK=60, 180, 50, XB=   215.0, 275.0, -60.0, 120.0, 0.0, 50.0 /';
for m=1:numel(meshXB)
    fprintf(fid, '%s\n',meshXB(m).line);
end

fprintf(fid,'\n');

%% 3rd entry & MISC
fprintf(fid,'&MISC RESTART= %s / # It could be true \n',conf.misc.restart);
fprintf(fid,'&MISC MEAN_FORCING(1:2)=.TRUE.,.TRUE., U0=%2.2f, V0=%2.2f, DT_MEAN_FORCING=0.1 / \n',...
    conf.misc.wind);
fprintf(fid,'\n');

%% 4th entry &TIME
fprintf(fid,'&TIME T_END=%3.1f, /\n',conf.time);
fprintf(fid,'\n');

%% 5th entry &SPEC
fprintf(fid,'&SPEC ID=''%s'' /\n',conf.specid);
fprintf(fid,'\n');

%% 6th entry rlse &SURF
fprintf(fid,'&SURF ID=''%s'', SPEC_ID=''%s'', MASS_FLUX(1)=%2.4f, RAMP_MF(1)=''%s'' /\n', ...
    conf.rlse.surfid, conf.specid, conf.rlse.mass, conf.rlse.rampid);

M=size(conf.rlse.ramp,1);
for i=1:M;
    fprintf(fid,'&RAMP ID=''%s'', T=%3.1f, F=%1.3f, /\n',conf.rlse.rampid, conf.rlse.ramp(i,:));
end
fprintf(fid,'\n');

fprintf(fid,'&VENT XB= %4.1f, %4.1f, %4.1f, %4.1f, %4.1f, %4.1f, SURF_ID=''%s'', COLOR=''%s'' /\n',...
    conf.rlse.ventXB,conf.rlse.surfid,conf.rlse.color);
fprintf(fid,'\n');

%% 7th entry wind &SURF
% fprintf(fid,'&SURF ID=''%s'', TMP_FRONT=%2.1f , VEL=%2.1f, VEL_T=%2.1f, %2.1f, PROFILE=''%s'', Z0=%2.1f, PLE=0.094, RAMP_V=''%s'' /\n', ...
%     char(conf.wind.surfid(1)),conf.wind.temp, conf.wind.VX, conf.wind.profile, conf.wind.Z0, conf.wind.rampid);
% fprintf(fid,'&SURF ID=''%s'', TMP_FRONT=%2.1f , VEL=%2.1f, VEL_T=%2.1f, %2.1f, PROFILE=''%s'', Z0=%2.1f, PLE=0.094, RAMP_V=''%s'' /\n', ...
%     char(conf.wind.surfid(2)),conf.wind.temp, conf.wind.VY, conf.wind.profile, conf.wind.Z0, conf.wind.rampid);
%
% [M,N]=size(conf.rlse.ramp);
% for i=1:M
%     fprintf(fid,'&RAMP ID=''%s'', T=%3.1f, F=%1.3f, /\n',conf.wind.rampid,conf.wind.ramp(i,:));
% end
% fprintf(fid,'\n');
axis_cfg={'XMIN','XMAX','YMIN','YMAX','ZMAX'};
for i=1:5   %conf.wind.vent is integer numbered  1|2|3
    fprintf(fid,'&VENT MB=''%s'', SURF_ID=''%s'' /\n',char(axis_cfg(i)),'OPEN');
end
fprintf(fid,'\n');

%% 7th entry &DUMP
fprintf(fid,'&DUMP MASS_FILE=%s  /\n',conf.dump.massfile);
% fprintf(fid,'&DUMP PLOT3D_SPEC_ID(5)=''%s'', PLOT3D_QUANTITY(1:5)=''TEMPERATURE'',''U-VELOCITY'',''V-VELOCITY'',''W-VELOCITY'',''%s'', DT_PL3D=%3.1f /\n', ...
%     conf.specid,conf.dump.pl3d_quantity,conf.dump.dt_pl3d);
fprintf(fid,'\n');

%% 8th entry &ISOF
% fprintf(fid,'&ISOF SPEC_ID=''%s'', QUANTITY=''%s'', VALUE(1)=%2.2e, VALUE(2)=%2.2e, VALUE(3)=%2.2e   /\n',...
%     conf.specid, conf.isof.quantity, conf.isof.value);
% fprintf(fid,'\n');


%% 9th entry &SLCF
fprintf(fid,'\n');
for i=1:numel(conf.slcf)
    fprintf(fid,'&SLCF %s =%3.1f, QUANTITY=''VELOCITY'', VECTOR=.TRUE. /\n', conf.slcf(i).surface, conf.slcf(i).PB);
end
fprintf(fid,'\n');
for i=1:numel(conf.slcf)
    fprintf(fid,'&SLCF %s =%3.1f, QUANTITY=''%s'', SPEC_ID=''%s'' /\n',conf.slcf(i).surface, conf.slcf(i).PB,...
        conf.slcf(i).quantity, conf.specid);
end
fprintf(fid,'\n');

%% 10th entry: &DEVC
load(conf.devc.xyz_config)
for i=1:numel(Devc_config)
    if conf.mesh(Devc_config(i).zone) ==1
        fprintf(fid,'&DEVC ID=''%s'',XYZ=%3.1f,%3.1f,%3.1f, QUANTITY=''%s'', SPEC_ID=''%s'' /\n',...
            Devc_config(i).id, Devc_config(i).xyz, conf.devc.quantity, conf.specid);
    end
end
fprintf(fid,'\n');

%% 11th entry: &OBST

% Writing Rectangle Data
load(conf.obst.rec_config);
% surf_id1='''INERT''';
% color1 ='''OLIVE DRAB''';
fprintf(fid,'%s\n','/Starting writing rectangular OBST');
N=size(Rec_Data,1);
if dike_flag==1
    for i=1:N
        fprintf(fid,'&OBST XB=%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f, COLOR=''%s'', SURF_ID=''%s''  /\n',...
            Rec_Data(i,:),conf.obst.rec_color, conf.obst.rec_surfid);
    end
elseif dike_flag==0
    for i=1:10
        fprintf(fid,'&OBST XB=%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f, COLOR=''%s'', SURF_ID=''%s''  /\n',...
            Rec_Data(i,:),conf.obst.rec_color, conf.obst.rec_surfid);
    end
end
fprintf(fid,'\n');

% Writing Cylinder Data
load(conf.obst.cyl_config);
% N=size(Cylinder_Data,1);
% for i=1:N
general_config.fid=fid;
general_config.color=conf.obst.cyl_color;
general_config.surf_id=conf.obst.cyl_surfid;
fun_obst_cylinder(Cylinder_Data,general_config)
% end
fprintf(fid,'\n');

%% Last entry: &TAIL
fprintf(fid,'&TAIL /\n');

fclose(fid);
disp(['FDS Script ',conf.head.chid,'.fds has been generated']);
end