function Obst=fun_obstcode_2D(parameter)
% This file calculate the obst code for 2D domain
IJ=parameter.IJ;
I=IJ(1);
J=IJ(2);
ii=2:I-1;
jj=2:J-1;
Obstcfg=parameter.Obstcfg;
% Initialization
Obstarea=zeros(I,J);
Obstcode=zeros(I,J);

[M,N]=size(Obstcfg);
for i=1:M
    for j=1:N/2
        Obstidx{i,j}=Obstcfg(i,2*j-1):Obstcfg(i,2*j);
    end    
    Obstarea(Obstidx{i,1},Obstidx{i,2})=1;
end

for l=1:M
    iii=Obstcfg(l,1)-1:Obstcfg(l,2)+1;
    jjj=Obstcfg(l,3)-1:Obstcfg(l,4)+1;
    for i=iii
        for j=jjj
            codstr=[num2str(Obstarea(i-1,j+1)),...
            num2str(Obstarea(i,j+1)),...
            num2str(Obstarea(i+1,j+1)),...
            num2str(Obstarea(i-1,j)),...
            num2str(Obstarea(i,j)),...
            num2str(Obstarea(i+1,j)),...
            num2str(Obstarea(i-1,j-1)),...
            num2str(Obstarea(i,j-1)),...
            num2str(Obstarea(i+1,j-1))];
        Obstcode(i,j)=bin2dec(codstr);
        end
    end
end

Obst=Obstcode;

