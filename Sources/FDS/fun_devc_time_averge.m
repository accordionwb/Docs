function outdata=fun_devc_time_averge(timestep,indata)
% This function reads individual DEVC data and calculate averaged
% concentration with in limited timestep
% Indata is single DEVC time series concentration with time segamentations
T_start=min(indata(:,1));
T_end=max(indata(:,1));
seg=floor((T_end-T_start)/timestep);
k=1; % Original Data index
N=size(indata,1);
for i=1:seg+1
    n=0;
    csum=0;
    while indata(k,1)<=T_start+i*timestep
        csum=csum+indata(k,2);
        n=n+1;
        if k<N
            k=k+1;
        else
            break
        end
    end
    out_t(i)=T_start+(i-1)*timestep+timestep/2;
    out_c(i)=csum/n;
    disp(['Time Block: ',num2str(i),' | Data point ',num2str(k),' Finished'])
    if k>=N
       disp(['Data index ',num2str(k),' exceed maximum data length ',num2str(N)])
       break
    end
end
outdata=[out_t;out_c]';


% Plot Results
set(gcf,'color','w')
FTS1=18; % FontSize 1
TXS1=14;
TXS2=12;
plot(out_t,out_c,'k-','LineWidth',2);
% Label and Font
title([num2str(timestep), 's Time Average Concentraion'],'FontSize',FTS1)
ylabel('Chlorine (ppm)','FontSize',FTS1)
xlabel('Time (s)','FontSize',FTS1)
