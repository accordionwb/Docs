function fun_case_report(conf)
%% This function read FDS configure file and generate a formatted text report file as csv format.
filename='Spectra_configuration.csv';
fid=fopen(filename,'w');

fprintf(fid,'Spectra simulation cases. Specie=''%s'', Simulation time=''%d''\n',conf(1).specid,conf(1).time);
fprintf(fid,'CHID,Wind U,Wind V, Release Rate (kg/s),SourceX,SourceY,SourceZ\n');
for i=1:length(conf)
    fprintf(fid,'%s,',conf(i).head.chid);
    fprintf(fid,'%2.4f,',conf(i).misc.wind(1));
    fprintf(fid,'%2.4f,',conf(i).misc.wind(2));
    fprintf(fid,'%2.4f,',conf(i).rlse.mass);
    fprintf(fid,'%3.1f,',mean(conf(i).rlse.ventXB(1:2)));
    fprintf(fid,'%3.1f,',mean(conf(i).rlse.ventXB(3:4)));
    fprintf(fid,'%3.1f',mean(conf(i).rlse.ventXB(5:6)));
    fprintf(fid,'\n');
end
fclose(fid);
