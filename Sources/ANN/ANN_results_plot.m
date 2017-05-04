function [Static]=ANN_results_plot(Check,simout,Info,key)
%% Introduction
% Check and simout is are the default input of this program.
% Info contains the information about the input parameter.
% key is the mode of post processing and plot.
%    key:  plot (give the plot of concentration, time W1 W2 L_max L_F
%    key:  static  (give the relative error of concentration and time, and the possiablity analysis of W1/W2 and L_F/L_max

%% Parameter Preparation
status=Info.status;
stability=Info.stability;


%% Part I.  Concentration | Time | W1/W2 | L1/L2
% Direct plot against phast output, to show the
if (strcmp(key,'plot'))
    switch status
        case 0  %Inputs,1:11 and Output, 5:6, ( concentration + all_pressure + all_stability )
            % Plot configuration: plot concentration and arrival time against phast results.
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            time_net=[datestr(now,30),'plot_11-2_CT_STA-',stability,'.mat'];
            figure(1)    % The first figure of concentration and time parity plots
            set(gcf,'position',[100,300,1500,600],'color','w')
            subplot(121)
            plotAdded(mdl1)
            title('(a)  Concentration (ppm)');xlabel('Phast output');ylabel('Neural network output')
            subplot(122)
            plotAdded(mdl2)
            title('(b)  Arrival time (s)');xlabel('Phast output');ylabel('Neural network output')
            Static.MSE_concentration=mdl1.MSE;
            Static.Rsquared_concentration=mdl1.Rsquared;
            Static.MSE_arrival_time=mdl2.MSE;
            Static.Rsuared_arrival_time=mdl2.Rsquared;
            fun_ANN_Sim_Frqc(Check,simout)  % plot frequency density
            save(['../Run_Save/',time_net],'Static','X1','X2','Y1','Y2','mdl1','mdl2')
            
        case 1  % Inputs,1:11 and Output, 1:4, ( L1/L2 & W1/W2 + all_pressure + all_stability )
            % Plot: calculate L1/L2, W1/W2
            X1=Check(1,:)';X2=Check(2,:)';
            X3=Check(3,:)';X4=Check(4,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            Y3=simout(3,:)';Y4=simout(4,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            mdl3=fitlm(X3,Y3);
            mdl4=fitlm(X4,Y4);
            time_net=[datestr(now,30),'_plot_11-4_LW_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,100,1500,900],'color','w')
            subplot(221)
            plotAdded(mdl1)
            title('(a)   W_1 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(222)
            plotAdded(mdl2)
            title('(b)  W_2 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(223)
            plotAdded(mdl3)
            title('(c)  L_F (m)');xlabel('Phast output L_target');ylabel('Neural network output')
            subplot(224)
            plotAdded(mdl4)
            title('(d)  L_max (m)');xlabel('Phast output L_max');ylabel('Neural network output')
            Static.MSE_W1=mdl1.MSE;
            Static.Rsquared_W1=mdl1.Rsquared;
            Static.MSE_W2=mdl2.MSE;
            Static.Rsuared_W2=mdl2.Rsquared;
            Static.MSE_Lf=mdl3.MSE;
            Static.Rsuared_Lf=mdl3.Rsquared;
            Static.MSE_Lmax=mdl4.MSE;
            Static.Rsuared_Lmax=mdl4.Rsquared;
            save(['../Run_Save/',time_net],'Static','X1','X2','X3','X4','Y1','Y2','Y3','Y4','mdl1','mdl2','mdl3','mdl4')
            
        case 2    % Inputs,1:11 and Output, 1:2, ( W1/W2 + all_pressure + all_stability )
            % Plot W1 and W2 each
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            time_net=[datestr(now,30),'_plot_11-2_W1W2_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,300,1500,600],'color','w')
            subplot(121)
            plotAdded(mdl1)
            title('(a)   W_1 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(122)
            plotAdded(mdl2)
            title('(b)  W_2 (m)');xlabel('Phast output');ylabel('Neural network output')
            Static.MSE_W1=mdl1.MSE;
            Static.Rsquared_W1=mdl1.Rsquared;
            Static.MSE_W2=mdl2.MSE;
            Static.Rsuared_W2=mdl2.Rsquared;
            save(['../Run_Save/',time_net],'Static','X1','X2','Y1','Y2','mdl1','mdl2')
            
        case 3 % Inputs,[1:9 11] and Output, 5:6, ( concentration + single_pressure + all_stability )
            % Plot  Concentration and stability
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            time_net=[datestr(now,30),'_plot_9.11-2_CT_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,300,1500,600],'color','w')
            subplot(121)
            plotAdded(mdl1)
            title('(a)  Concentration (ppm)');xlabel('Phast output');ylabel('Neural network output')
            subplot(122)
            plotAdded(mdl2)
            title('(b)  Arrival time (s)');xlabel('Phast output');ylabel('Neural network output')
            Static.MSE_concentration=mdl1.MSE;
            Static.Rsquared_concentration=mdl1.Rsquared;
            Static.MSE_arrival_time=mdl2.MSE;
            Static.Rsuared_arrival_time=mdl2.Rsquared;
           fun_ANN_Sim_Frqc(Check,simout)  % plot frequency density
            save(['../Run_Save/',time_net],'Static','X1','X2','Y1','Y2','mdl1','mdl2')
            
        case 4  % Inputs,[1:9 11] and Output, 1:4, ( L1/L2 & W1/W2 + single_pressure + all_stability )
            % Plot W1 W2 L1 and L2 each
            X1=Check(1,:)';X2=Check(2,:)';
            X3=Check(3,:)';X4=Check(4,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            Y3=simout(3,:)';Y4=simout(4,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            mdl3=fitlm(X3,Y3);
            mdl4=fitlm(X4,Y4);
            time_net=[datestr(now,30),'_plot_9.11-4_LW_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,100,1500,900],'color','w')
            subplot(221)
            plotAdded(mdl1)
            title('(a)   W_1 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(222)
            plotAdded(mdl2)
            title('(b)  W_2 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(223)
            plotAdded(mdl3)
            title('(c)  L_F (m)');xlabel('Phast output L_target');ylabel('Neural network output')
            subplot(224)
            plotAdded(mdl4)
            title('(d)  L_max (m)');xlabel('Phast output L_max');ylabel('Neural network output')
            Static.MSE_W1=mdl1.MSE;
            Static.Rsquared_W1=mdl1.Rsquared;
            Static.MSE_W2=mdl2.MSE;
            Static.Rsuared_W2=mdl2.Rsquared;
            Static.MSE_Lf=mdl3.MSE;
            Static.Rsuared_Lf=mdl3.Rsquared;
            Static.MSE_Lmax=mdl4.MSE;
            Static.Rsuared_Lmax=mdl4.Rsquared;
            save(['../Run_Save/',time_net],'Static','X1','X2','X3','X4','Y1','Y2','Y3','Y4','mdl1','mdl2','mdl3','mdl4')
            
        case 5  % Inputs,1:11 and Output, 1:2, ( W1/W2 + single_pressure + all_stability )
            % Plot W1 and W2
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            time_net=[datestr(now,30),'_plot_11-2_W1W2_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,300,1500,600],'color','w')
            subplot(121)
            plotAdded(mdl1)
            title('(a)   W_1 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(122)
            plotAdded(mdl2)
            title('(b)  W_2 (m)');xlabel('Phast output');ylabel('Neural network output')
            Static.MSE_W1=mdl1.MSE;
            Static.Rsquared_W1=mdl1.Rsquared;
            Static.MSE_W2=mdl2.MSE;
            Static.Rsuared_W2=mdl2.Rsquared;
            save(['../Run_Save/',time_net],'Static','X1','X2','Y1','Y2','mdl1','mdl2')
            
        case 6 % Inputs,1:10 and Output, 5:6, ( concentration + all_pressure + single_stability )
            % Plot  Concentration and stability
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            time_net=[datestr(now,30),'_plot_10-2_CT_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,300,1500,600],'color','w')
            subplot(121)
            plotAdded(mdl1)
            title('(a)  Concentration (ppm)');xlabel('Phast output');ylabel('Neural network output')
            subplot(122)
            plotAdded(mdl2)
            title('(b)  Arrival time (s)');xlabel('Phast output');ylabel('Neural network output')
            Static.MSE_concentration=mdl1.MSE;
            Static.Rsquared_concentration=mdl1.Rsquared;
            Static.MSE_arrival_time=mdl2.MSE;
            Static.Rsuared_arrival_time=mdl2.Rsquared;
            fun_ANN_Sim_Frqc(Check,simout)  % plot frequency density
            save(['../Run_Save/',time_net],'Static','X1','X2','Y1','Y2','mdl1','mdl2')
            
        case 7  % Inputs,1:10 and Output, 1:4, ( L1/L2 & W1/W2 + all_pressure + single_stability )
            % Plot W1 W2 L1 and L2 each
            X1=Check(1,:)';X2=Check(2,:)';
            X3=Check(3,:)';X4=Check(4,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            Y3=simout(3,:)';Y4=simout(4,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            mdl3=fitlm(X3,Y3);
            mdl4=fitlm(X4,Y4);
            time_net=[datestr(now,30),'_plot_10-4_LW_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,100,1500,900],'color','w')
            subplot(221)
            plotAdded(mdl1)
            title('(a)   W_1 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(222)
            plotAdded(mdl2)
            title('(b)  W_2 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(223)
            plotAdded(mdl3)
            title('(c)  L_F (m)');xlabel('Phast output L_target');ylabel('Neural network output')
            subplot(224)
            plotAdded(mdl4)
            title('(d)  L_max (m)');xlabel('Phast output L_max');ylabel('Neural network output')
            Static.MSE_W1=mdl1.MSE;
            Static.Rsquared_W1=mdl1.Rsquared;
            Static.MSE_W2=mdl2.MSE;
            Static.Rsuared_W2=mdl2.Rsquared;
            Static.MSE_Lf=mdl3.MSE;
            Static.Rsuared_Lf=mdl3.Rsquared;
            Static.MSE_Lmax=mdl4.MSE;
            Static.Rsuared_Lmax=mdl4.Rsquared;
            save(['../Run_Save/',time_net],'Static','X1','X2','X3','X4','Y1','Y2','Y3','Y4','mdl1','mdl2','mdl3','mdl4')
            
        case 8  % Inputs,1:10 and Output, 1:2, ( W1/W2 + all_pressure + single_stability )
            % Plot
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            time_net=[datestr(now,30),'_plot_10-2_W1W2_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,300,1500,600],'color','w')
            subplot(121)
            plotAdded(mdl1)
            title('(a)   W_1 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(122)
            plotAdded(mdl2)
            title('(b)  W_2 (m)');xlabel('Phast output');ylabel('Neural network output')
            Static.MSE_W1=mdl1.MSE;
            Static.Rsquared_W1=mdl1.Rsquared;
            Static.MSE_W2=mdl2.MSE;
            Static.Rsuared_W2=mdl2.Rsquared;
            save(['../Run_Save/',time_net],'Static','X1','X2','Y1','Y2','mdl1','mdl2')
            
        case 9 % Inputs,1:9 and Output, 5:6, ( concentration + single_pressure + single_stability )
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            time_net=[datestr(now,30),'_plot_9-2_CT_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,300,1500,600],'color','w')
            subplot(121)
            plotAdded(mdl1)
            title('(a)  Concentration (ppm)');xlabel('Phast output');ylabel('Neural network output')
            subplot(122)
            plotAdded(mdl2)
            title('(b)  Arrival time (s)');xlabel('Phast output');ylabel('Neural network output')
            Static.MSE_concentration=mdl1.MSE;
            Static.Rsquared_concentration=mdl1.Rsquared;
            Static.MSE_arrival_time=mdl2.MSE;
            Static.Rsuared_arrival_time=mdl2.Rsquared;
            fun_ANN_Sim_Frqc(Check,simout)  % plot frequency density
            save(['../Run_Save/',time_net],'Static','X1','X2','Y1','Y2','mdl1','mdl2')
            
        case 10  % Inputs,1:9 and Output, 1:4, ( L1/L2 & W1/W2 + single_pressure + single_stability )
            % Plot W1 W2 L1 and L2 each
            X1=Check(1,:)';X2=Check(2,:)';
            X3=Check(3,:)';X4=Check(4,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            Y3=simout(3,:)';Y4=simout(4,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            mdl3=fitlm(X3,Y3);
            mdl4=fitlm(X4,Y4);
            time_net=[datestr(now,30),'_plot_9-4_LW_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,100,1500,900],'color','w')
            subplot(221)
            plotAdded(mdl1)
            title('(a)   W_1 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(222)
            plotAdded(mdl2)
            title('(b)  W_2 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(223)
            plotAdded(mdl3)
            title('(c)  L_F (m)');xlabel('Phast output L_target');ylabel('Neural network output')
            subplot(224)
            plotAdded(mdl4)
            title('(d)  L_max (m)');xlabel('Phast output L_max');ylabel('Neural network output')
            Static.MSE_W1=mdl1.MSE;
            Static.Rsquared_W1=mdl1.Rsquared;
            Static.MSE_W2=mdl2.MSE;
            Static.Rsuared_W2=mdl2.Rsquared;
            Static.MSE_Lf=mdl3.MSE;
            Static.Rsuared_Lf=mdl3.Rsquared;
            Static.MSE_Lmax=mdl4.MSE;
            Static.Rsuared_Lmax=mdl4.Rsquared;
            save(['../Run_Save/',time_net],'Static','X1','X2','X3','X4','Y1','Y2','Y3','Y4','mdl1','mdl2','mdl3','mdl4')
            
        case 11 % Inputs,1:9 and Output, 1:2, ( W1/W2 + single_pressure + single_stability )
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            mdl1=fitlm(X1,Y1);
            mdl2=fitlm(X2,Y2);
            time_net=[datestr(now,30),'_plot_9-2_W1W2_STA-',stability,'.mat'];
            figure(1)
            set(gcf,'position',[100,300,1500,600],'color','w')
            subplot(121)
            plotAdded(mdl1)
            title('(a)   W_1 (m)');xlabel('Phast output');ylabel('Neural network output')
            subplot(122)
            plotAdded(mdl2)
            title('(b)  W_2 (m)');xlabel('Phast output');ylabel('Neural network output')
            Static.MSE_W1=mdl1.MSE;
            Static.Rsquared_W1=mdl1.Rsquared;
            Static.MSE_W2=mdl2.MSE;
            Static.Rsuared_W2=mdl2.Rsquared;
            save(['../Run_Save/',time_net],'Static','X1','X2','Y1','Y2','mdl1','mdl2')
    end
end
%%  Part II.  Concentration Analysis, Calculating the relative error.

if (strcmp(key, 'static'))
    switch status
        case 0 %Inputs,1:11 and Output, 5:6, ( concentration + all_pressure + all_stability )
            % Plot configuration: plot concentration and arrival time against phast results.
            X1=Check(1,:)';Y1=simout(1,:)';   % concentration
            X2=Check(2,:)';Y2=simout(2,:)';   % arrival time
            relative_error1=(Y1-X1)./X1*100;  mean_r1=mean(relative_error1); max_r1=max(abs(relative_error1));
            relative_error2=(Y2-X2)./X2*100;  mean_r2=mean(relative_error2); max_r2=max(abs(relative_error2));
            figure(1)
            set(gcf,'position',[100,200,1400,900],'color','w')
            subplot(2,1,1)
            bar(relative_error1)
            title('(a)  Relative error of ANN calculated chlorine concentration against phast results')
            xlabel(['Number of cases, total ',num2str(length(X1))])
            ylabel('Relative error, %')
            subplot(2,1,2)
            bar(relative_error2)
            title('(a)  Relative error of ANN calculated chlorine arrival time against phast results')
            xlabel(['Number of cases, total ',num2str(length(X1))])
            ylabel('Relative error, %')
            Static.mean_error_concentration=mean_r1;
            Static.max_error_concentration=max_r1;
            Static.mean_error_arrival_time=mean_r2;
            Static.max_error_arrival_time=max_r2;
            fun_ANN_Sim_Frqc(Check,simout)  % plot frequency density 
            time_str=[datastr(now,30),'static_11-2_CT_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','Y1','Y2','relative_error1','relative_error2')
            
        case 1  % Inputs,1:11 and Output, 1:4, ( L1/L2 & W1/W2 + all_pressure + all_stability )
            X1=Check(1,:)';X2=Check(2,:)';            X3=Check(3,:)';X4=Check(4,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';            Y3=simout(3,:)';Y4=simout(4,:)';
            net_width_ratio=Y1./Y2;
            net_distance_ratio=Y3./Y4;
            phast_width_ratio=X1./X2;
            phast_distance_ratio=X3./X4;
            YN_width=boolean(floor(net_width_ratio));
            YF_width=boolean(floor(phast_width_ratio));
            YN_distance=boolean(floor(net_distance_ratio));
            YF_distance=boolean(floor(phast_distance_ratio));
            compare_width=YN_width-YF_width;  % results: =0, equal; >0 Miss; <0, False
            compare_distance=YN_distance-YF_distance;
            figure(1)    %  plot of ratios
            set(gcf,'position',[100,50,1300,900],'color','w')
            subplot(1,2,1)
            plot(net_distance_ratio,net_width_ratio,'kx')
            title('(a)  Neural network estimation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            subplot(1,2,2)
            plot(phast_distance_ratio,phast_width_ratio,'kx')
            title('(b)  Phast calculation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            figure(2)  % plot of statistic results
            set(gcf,'position',[100,100,1300,500],'color','w')
            subplot(2,1,1)
            bar(compare_distance)
            title('(a)   [(L_F/L_{max})_{net}]_f-[(L_F/L_{max})_{phast}]_f')
            xlabel([' Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            subplot(2,1,2)
            bar(compare_width)
            title('(b)   [(W_1/W_2)_{net}]_f-[(W_1/W_2)_{phast}]_f')
            xlabel(['Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            Static.FN_d=find(compare_distance>0);  % false negative
            Static.FP_d=find(compare_distance<0);
            Static.FN_w=find(compare_width>0);
            Static.FP_w=find(compare_width<0);
            Static.Co_w=find(compare_width==0);  % Correct width ratio
            Static.Co_d=find(compare_distance==0);   % Correct distance ratio
            lens=length(compare_distance);
            text(0.2*lens,-0.5,['False prediction: ',num2str(length(false))])
            text(0.2*lens,0.5,['Miss Prediction: ',num2str(length(miss))])
            axis([0,1.1*length(compare_width),-1,1])
            time_str=[datastr(now,30),'static_11-4_LW_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','X3','X4','Y1','Y2','Y3','Y4','YN_*','FN_*','compare_*')
            
        case 2   % Inputs,1:11 and Output, 1:2, ( W1/W2 + all_pressure + all_stability )
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            net_width_ratio=Y1./Y2;
            phast_width_ratio=X1./X2;
            YN_width=boolean(floor(net_width_ratio));
            YF_width=boolean(floor(phast_width_ratio));
            compare_width=YN_width-YF_width;  % results: =0, equal; >0 Miss; <0, False
            figure(1)
            set(gcf,'position',[100,50,1300,900],'color','w')
            plot(net_distance_ratio,net_width_ratio,'kx')
            title('(a)  Neural network estimation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            figure(2)
            set(gcf,'position',[100,100,1300,500],'color','w')
            bar(compare_width)
            title('(b)   [(W_1/W_2)_{net}]_f-[(W_1/W_2)_{phast}]_f')
            xlabel(['Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            Static.FN_w=find(compare_width>0);
            Static.FP_w=find(compare_width<0);
            Static.Co_w=find(compare_width==0);  % Correct width ratio
            lens=length(compare_width);
            text(0.2*lens,-0.5,['False prediction: ',num2str(length(false))])
            text(0.2*lens,0.5,['Miss Prediction: ',num2str(length(miss))])
            axis([0,1.05*length(compare_width),-1,1])
            time_str=[datastr(now,30),'static_11-2_W1W2_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','Y1','Y2','YN_width','FN_width','compare_width')
            
        case 3  % Inputs,[1:9 11] and Output, 5:6, ( concentration + single_pressure + all_stability )
            % Plot  Concentration and stability
            X1=Check(1,:)';Y1=simout(1,:)';   % concentration
            X2=Check(2,:)';Y2=simout(2,:)';   % arrival time
            relative_error1=(Y1-X1)./X1*100;  mean_r1=mean(relative_error1); max_r1=max(abs(relative_error1));
            relative_error2=(Y2-X2)./X2*100;  mean_r2=mean(relative_error2); max_r2=max(abs(relative_error2));
            figure(1)
            set(gcf,'position',[100,200,1400,900],'color','w')
            subplot(2,1,1)
            bar(relative_error1)
            title('(a)  Relative error of ANN calculated chlorine concentration against phast results')
            xlabel(['Number of cases, total ',num2str(length(X1))])
            ylabel('Relative error, %')
            subplot(2,1,2)
            bar(relative_error2)
            title('(a)  Relative error of ANN calculated chlorine arrival time against phast results')
            xlabel(['Number of cases, total ',num2str(length(X1))])
            ylabel('Relative error, %')
            Static.mean_error_concentration=mean_r1;
            Static.max_error_concentration=max_r1;
            Static.mean_error_arrival_time=mean_r2;
            Static.max_error_arrival_time=max_r2;
            fun_ANN_Sim_Frqc(Check,simout)  % plot frequency density
            time_str=[datastr(now,30),'static_9.11-2_CT_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','Y1','Y2','relative_error1','relative_error2')
            
        case 4  % Inputs,[1:9 11] and Output, 1:4, ( L1/L2 & W1/W2 + single_pressure + all_stability )
            X1=Check(1,:)';X2=Check(2,:)';            X3=Check(3,:)';X4=Check(4,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';            Y3=simout(3,:)';Y4=simout(4,:)';
            net_width_ratio=Y1./Y2;
            net_distance_ratio=Y3./Y4;
            phast_width_ratio=X1./X2;
            phast_distance_ratio=X3./X4;
            YN_width=boolean(floor(net_width_ratio));
            YF_width=boolean(floor(phast_width_ratio));
            YN_distance=boolean(floor(net_distance_ratio));
            YF_distance=boolean(floor(phast_distance_ratio));
            compare_width=YN_width-YF_width;  % results: =0, equal; >0 Miss; <0, False
            compare_distance=YN_distance-YF_distance;
            figure(1)    %  plot of ratios
            set(gcf,'position',[100,50,1300,900],'color','w')
            subplot(1,2,1)
            plot(net_distance_ratio,net_width_ratio,'kx')
            title('(a)  Neural network estimation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            subplot(1,2,2)
            plot(phast_distance_ratio,phast_width_ratio,'kx')
            title('(b)  Phast calculation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            figure(2)  % plot of statistic results
            set(gcf,'position',[100,100,1300,500],'color','w')
            subplot(2,1,1)
            bar(compare_distance)
            title('(a)   [(L_F/L_{max})_{net}]_f-[(L_F/L_{max})_{phast}]_f')
            xlabel([' Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            subplot(2,1,2)
            bar(compare_width)
            title('(b)   [(W_1/W_2)_{net}]_f-[(W_1/W_2)_{phast}]_f')
            xlabel(['Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            Static.FN_d=find(compare_distance>0);  % false negative
            Static.FP_d=find(compare_distance<0);
            Static.FN_w=find(compare_width>0);
            Static.FP_w=find(compare_width<0);
            Static.Co_w=find(compare_width==0);  % Correct width ratio
            Static.Co_d=find(compare_distance==0);   % Correct distance ratio
            lens=length(compare_distance);
            text(0.2*lens,-0.5,['False prediction: ',num2str(length(false))])
            text(0.2*lens,0.5,['Miss Prediction: ',num2str(length(miss))])
            axis([0,1.1*length(compare_width),-1,1])
            time_str=[datastr(now,30),'static_9.11-4_LW_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','X3','X4','Y1','Y2','Y3','Y4','YN_*','FN_*','compare_*')
            
        case 5 % Inputs,1:11 and Output, 1:2, ( W1/W2 + single_pressure + all_stability )
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            net_width_ratio=Y1./Y2;
            phast_width_ratio=X1./X2;
            YN_width=boolean(floor(net_width_ratio));
            YF_width=boolean(floor(phast_width_ratio));
            compare_width=YN_width-YF_width;  % results: =0, equal; >0 Miss; <0, False
            figure(1)
            set(gcf,'position',[100,50,1300,900],'color','w')
            plot(net_distance_ratio,net_width_ratio,'kx')
            title('(a)  Neural network estimation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            figure(2)
            set(gcf,'position',[100,100,1300,500],'color','w')
            bar(compare_width)
            title('(b)   [(W_1/W_2)_{net}]_f-[(W_1/W_2)_{phast}]_f')
            xlabel(['Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            Static.FN_w=find(compare_width>0);
            Static.FP_w=find(compare_width<0);
            Static.Co_w=find(compare_width==0);  % Correct width ratio
            lens=length(compare_width);
            text(0.2*lens,-0.5,['False prediction: ',num2str(length(false))])
            text(0.2*lens,0.5,['Miss Prediction: ',num2str(length(miss))])
            axis([0,1.05*length(compare_width),-1,1])
            time_str=[datastr(now,30),'static_11-2_W1W2_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','Y1','Y2','YN_width','FN_width','compare_width')
            
        case 6  % Inputs,1:10 and Output, 5:6, ( concentration + all_pressure + single_stability )
            X1=Check(1,:)';Y1=simout(1,:)';   % concentration
            X2=Check(2,:)';Y2=simout(2,:)';   % arrival time
            relative_error1=(Y1-X1)./X1*100;  mean_r1=mean(relative_error1); max_r1=max(abs(relative_error1));
            relative_error2=(Y2-X2)./X2*100;  mean_r2=mean(relative_error2); max_r2=max(abs(relative_error2));
            figure(1)
            set(gcf,'position',[100,200,1400,900],'color','w')
            subplot(2,1,1)
            bar(relative_error1)
            title('(a)  Relative error of ANN calculated chlorine concentration against phast results')
            xlabel(['Number of cases, total ',num2str(length(X1))])
            ylabel('Relative error, %')
            subplot(2,1,2)
            bar(relative_error2)
            title('(a)  Relative error of ANN calculated chlorine arrival time against phast results')
            xlabel(['Number of cases, total ',num2str(length(X1))])
            ylabel('Relative error, %')
            Static.mean_error_concentration=mean_r1;
            Static.max_error_concentration=max_r1;
            Static.mean_error_arrival_time=mean_r2;
            Static.max_error_arrival_time=max_r2;
            fun_ANN_Sim_Frqc(Check,simout)   % plot frequency density
            time_str=[datastr(now,30),'static_10-2_CT_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','Y1','Y2','relative_error1','relative_error2')
            
        case 7  % Inputs,1:10 and Output, 1:4, ( L1/L2 & W1/W2 + all_pressure + single_stability )
            X1=Check(1,:)';X2=Check(2,:)';            X3=Check(3,:)';X4=Check(4,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';            Y3=simout(3,:)';Y4=simout(4,:)';
            net_width_ratio=Y1./Y2;
            net_distance_ratio=Y3./Y4;
            phast_width_ratio=X1./X2;
            phast_distance_ratio=X3./X4;
            YN_width=boolean(floor(net_width_ratio));
            YF_width=boolean(floor(phast_width_ratio));
            YN_distance=boolean(floor(net_distance_ratio));
            YF_distance=boolean(floor(phast_distance_ratio));
            compare_width=YN_width-YF_width;  % results: =0, equal; >0 Miss; <0, False
            compare_distance=YN_distance-YF_distance;
            figure(1)    %  plot of ratios
            set(gcf,'position',[100,50,1300,900],'color','w')
            subplot(1,2,1)
            plot(net_distance_ratio,net_width_ratio,'kx')
            title('(a)  Neural network estimation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            subplot(1,2,2)
            plot(phast_distance_ratio,phast_width_ratio,'kx')
            title('(b)  Phast calculation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            figure(2)  % plot of statistic results
            set(gcf,'position',[100,100,1300,500],'color','w')
            subplot(2,1,1)
            bar(compare_distance)
            title('(a)   [(L_F/L_{max})_{net}]_f-[(L_F/L_{max})_{phast}]_f')
            xlabel([' Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            subplot(2,1,2)
            bar(compare_width)
            title('(b)   [(W_1/W_2)_{net}]_f-[(W_1/W_2)_{phast}]_f')
            xlabel(['Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            Static.FN_d=find(compare_distance>0);  % false negative
            Static.FP_d=find(compare_distance<0);
            Static.FN_w=find(compare_width>0);
            Static.FP_w=find(compare_width<0);
            Static.Co_w=find(compare_width==0);  % Correct width ratio
            Static.Co_d=find(compare_distance==0);   % Correct distance ratio
            lens=length(compare_distance);
            text(0.2*lens,-0.5,['False prediction: ',num2str(length(false))])
            text(0.2*lens,0.5,['Miss Prediction: ',num2str(length(miss))])
            axis([0,1.1*length(compare_width),-1,1])
            time_str=[datastr(now,30),'static_10-4_LW_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','X3','X4','Y1','Y2','Y3','Y4','YN_*','FN_*','compare_*')
            
        case 8  % Inputs,1:10 and Output, 1:2, ( W1/W2 + all_pressure + single_stability )
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            net_width_ratio=Y1./Y2;
            phast_width_ratio=X1./X2;
            YN_width=boolean(floor(net_width_ratio));
            YF_width=boolean(floor(phast_width_ratio));
            compare_width=YN_width-YF_width;  % results: =0, equal; >0 Miss; <0, False
            figure(1)
            set(gcf,'position',[100,50,1300,900],'color','w')
            plot(net_distance_ratio,net_width_ratio,'kx')
            title('(a)  Neural network estimation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            figure(2)
            set(gcf,'position',[100,100,1300,500],'color','w')
            bar(compare_width)
            title('(b)   [(W_1/W_2)_{net}]_f-[(W_1/W_2)_{phast}]_f')
            xlabel(['Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            Static.FN_w=find(compare_width>0);
            Static.FP_w=find(compare_width<0);
            Static.Co_w=find(compare_width==0);  % Correct width ratio
            lens=length(compare_width);
            text(0.2*lens,-0.5,['False prediction: ',num2str(length(false))])
            text(0.2*lens,0.5,['Miss Prediction: ',num2str(length(miss))])
            axis([0,1.05*length(compare_width),-1,1])
            time_str=[datastr(now,30),'static_10-2_W1W2_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','Y1','Y2','YN_width','FN_width','compare_width')
            
        case 9 % Inputs,1:9 and Output, 5:6, ( concentration + single_pressure + single_stability )
            X1=Check(1,:)';Y1=simout(1,:)';   % concentration
            X2=Check(2,:)';Y2=simout(2,:)';   % arrival time
            relative_error1=(Y1-X1)./X1*100;  mean_r1=mean(relative_error1); max_r1=max(abs(relative_error1));
            relative_error2=(Y2-X2)./X2*100;  mean_r2=mean(relative_error2); max_r2=max(abs(relative_error2));
            figure(1)
            set(gcf,'position',[100,200,1400,900],'color','w')
            subplot(2,1,1)
            bar(relative_error1)
            title('(a)  Relative error of ANN calculated chlorine concentration against phast results')
            xlabel(['Number of cases, total ',num2str(length(X1))])
            ylabel('Relative error, %')
            subplot(2,1,2)
            bar(relative_error2)
            title('(a)  Relative error of ANN calculated chlorine arrival time against phast results')
            xlabel(['Number of cases, total ',num2str(length(X1))])
            ylabel('Relative error, %')
            Static.mean_error_concentration=mean_r1;
            Static.max_error_concentration=max_r1;
            Static.mean_error_arrival_time=mean_r2;
            Static.max_error_arrival_time=max_r2;
                        fun_ANN_Sim_Frqc(Check,simout)
            time_str=[datestr(now,30),'static_9-2_CT_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','Y1','Y2','relative_error1','relative_error2')
            
        case 10   % Inputs,1:9 and Output, 1:4, ( L1/L2 & W1/W2 + single_pressure + single_stability )
            X1=Check(1,:)';X2=Check(2,:)';            X3=Check(3,:)';X4=Check(4,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';            Y3=simout(3,:)';Y4=simout(4,:)';
            net_width_ratio=Y1./Y2;
            net_distance_ratio=Y3./Y4;
            phast_width_ratio=X1./X2;
            phast_distance_ratio=X3./X4;
            YN_width=boolean(floor(net_width_ratio));
            YF_width=boolean(floor(phast_width_ratio));
            YN_distance=boolean(floor(net_distance_ratio));
            YF_distance=boolean(floor(phast_distance_ratio));
            compare_width=YN_width-YF_width;  % results: =0, equal; >0 Miss; <0, False
            compare_distance=YN_distance-YF_distance;
            figure(1)    %  plot of ratios
            set(gcf,'position',[100,50,1300,900],'color','w')
            subplot(1,2,1)
            plot(net_distance_ratio,net_width_ratio,'kx')
            title('(a)  Neural network estimation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            subplot(1,2,2)
            plot(phast_distance_ratio,phast_width_ratio,'kx')
            title('(b)  Phast calculation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            figure(2)  % plot of statistic results
            set(gcf,'position',[100,100,1300,500],'color','w')
            subplot(2,1,1)
            bar(compare_distance)
            title('(a)   [(L_F/L_{max})_{net}]_f-[(L_F/L_{max})_{phast}]_f')
            xlabel([' Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            subplot(2,1,2)
            bar(compare_width)
            title('(b)   [(W_1/W_2)_{net}]_f-[(W_1/W_2)_{phast}]_f')
            xlabel(['Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            Static.FN_d=find(compare_distance>0);  % false negative
            Static.FP_d=find(compare_distance<0);
            Static.FN_w=find(compare_width>0);
            Static.FP_w=find(compare_width<0);
            Static.Co_w=find(compare_width==0);  % Correct width ratio
            Static.Co_d=find(compare_distance==0);   % Correct distance ratio
            lens=length(compare_distance);
            text(0.2*lens,-0.5,['False prediction: ',num2str(length(false))])
            text(0.2*lens,0.5,['Miss Prediction: ',num2str(length(miss))])
            axis([0,1.1*length(compare_width),-1,1])
            time_str=[datestr(now,30),'static_9-4_LW_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','X3','X4','Y1','Y2','Y3','Y4','YN_*','FN_*','compare_*')
            
        case 11  % Inputs,1:9 and Output, 1:2, ( W1/W2 + single_pressure + single_stability )
            X1=Check(1,:)';X2=Check(2,:)';
            Y1=simout(1,:)';Y2=simout(2,:)';
            net_width_ratio=Y1./Y2;
            phast_width_ratio=X1./X2;
            YN_width=boolean(floor(net_width_ratio));
            YF_width=boolean(floor(phast_width_ratio));
            compare_width=YN_width-YF_width;  % results: =0, equal; >0 Miss; <0, False
            figure(1)
            set(gcf,'position',[100,50,1300,900],'color','w')
            plot(net_distance_ratio,net_width_ratio,'kx')
            title('(a)  Neural network estimation');xlabel('L_F/L_{max}');ylabel('W_1/W_2')
            figure(2)
            set(gcf,'position',[100,100,1300,500],'color','w')
            bar(compare_width)
            title('(b)   [(W_1/W_2)_{net}]_f-[(W_1/W_2)_{phast}]_f')
            xlabel(['Number of cases, total ',num2str(length(simout))])
            ylabel('Fault indicator')
            Static.FN_w=find(compare_width>0);
            Static.FP_w=find(compare_width<0);
            Static.Co_w=find(compare_width==0);  % Correct width ratio
            lens=length(compare_width);
            text(0.2*lens,-0.5,['False prediction: ',num2str(length(false))])
            text(0.2*lens,0.5,['Miss Prediction: ',num2str(length(miss))])
            axis([0,1.05*length(compare_width),-1,1])
            time_str=[datestr(now,30),'static_9-2_W1W2_STA-',stability,'.mat'];
            save(['../Run_Save/',time_str],'Static','X1','X2','Y1','Y2','YN_width','FN_width','compare_width')
            
    end
    
end  % end if strcmp



