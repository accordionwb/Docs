function fun_ANN_Sim_Frqc(Check,simout,M)
% This function plot the relative frequency density of the relaive error of concentration and time
% estimated by neural network.
% Input:
%   Check and simout are Phast & ANN calculation results
%   M: the total number of divided on x-axes. default is 51
% Return:
%   plot of concentration and time frequency density figure.
%% Data preparation
X1=Check(1,:)';X2=Check(2,:)';
Y1=simout(1,:)';Y2=simout(2,:)';
C_r=(Y1-X1)./X1;
T_r=(Y2-X2)./X2;
N=length(X1);
% Range
Range_C=[1.1*min(C_r),1.1*max(C_r)];
Range_T=[1.1*min(T_r),1.1*max(T_r)];
% Division
if nargin<3
    M=[101,51];
end
Xc=linspace(Range_C(1),Range_C(2),M(1));
Xt=linspace(Range_T(1),Range_T(2),M(2));
for i=1:length(Xc)-1
    yc(i)=length(find( (C_r)>=Xc(i) & C_r<Xc(i+1) ))/N;
    xc(i)=mean([Xc(i),Xc(i+1)]);
end

for j=1:length(Xt)-1
    yt(j)=length(find( (T_r)>=Xt(j) & T_r<Xt(j+1) ))/N;
    xt(j)=mean([Xt(j),Xt(j+1)]);
end
% Plot figure
figure(2)    % Second figure for concentration and time distribution plots
set(gcf,'position',[100,300,1500,600],'color','w')
subplot(121)
plot(xc,yc)
%  title('(a)  Concentration accuracy distribution');xlabel('Relative deviation of concentration');ylabel('Frequency Density')
title('(a)  Concentration accuracy distribution');xlabel('(C_{net}-C_{phast})/C_{phast}');ylabel('Frequency Density')
subplot(122)
plot(xt,yt)
title('(b)  Time accuracy distribution');xlabel('(T_{net}-T_{phast})/T_{phast}');ylabel('Frequency Density')