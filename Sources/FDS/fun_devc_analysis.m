function devc_analysis=fun_devc_analysis(Head, Data, misc)
% IF_PLOT is boolean.
% Function returns 
%% Set path and read files  Static data that are used for domain generation
% addpath('D:\SVN\Code\Matlab\FDS\Spectra')
load('MAT_Rec_config.mat');
load('MAT_Cylinder_config.mat');
load('MAT_Devc_config.mat')

%% Read and plot active devc
N=size(Data,2);
T=Data(:,1);
used_devc=Data(:,2:N)*1e6; %convert to PPM

for j=1:N-1
    N_value(j)=norm(used_devc(:,j));
    devc_id(j)=Head(j+1);
end
index=find(N_value>=1);

devc_analysis(1).id_count=1;

for j=1:length(devc_id)
    if any(index==j)
        id_count=1;
        devc_analysis(j).id_count=1;
    else
        devc_analysis(j).id_count=0;
    end
    devc_analysis(j).id=devc_id(j);
    devc_analysis(j).time=T;
    devc_analysis(j).con=used_devc(:,j);
    for k=1:length(Devc_config)
        if strcmp(devc_id(j),Devc_config(k).id)
            devc_analysis(j).xyz=Devc_config(k).xyz;
        end
    end
end

% plot results
if misc.IF_PLOT
    for j=1:length(index)
        if exist(misc.foldername,'dir')
            disp([misc.foldername,' already exist'])
        else
            mkdir(misc.foldername)
            disp(['Create Folder: ',misc.foldername])
        end
        figurename1=[misc.foldername,char(devc_id(index(j))),'.png'];
        if exist(figurename1,'file')
        else
            figure(1)
            plot(T,used_devc(:,index(j)));
            tty=strrep(char(devc_id(index(j))),'_','\_'); % replace '_' with '\_'
            title(['DEVC "',tty,'" Output'])
            ylabel('ppm')
            xlabel('Time /s')
            saveas(gcf,figurename1);
            close(1)
        end
    end
    
    
    figurename2=[misc.foldername,'Overview'];
    if exist([figurename2,'.png'],'file')
        disp([figurename2,'.png already exist'])
    else
        X_Corner=Rec_Data(:,1);
        Y_Corner=Rec_Data(:,3);
        X_length=Rec_Data(:,2)-Rec_Data(:,1);
        Y_length=Rec_Data(:,4)-Rec_Data(:,3);
        
        % Cylinder Data reading
        
        % user_fds_cylinder_v2(sol,xc,yc,zc,dx,ro,l_cyl,axis,index)
        % Where:
        %   SOL = (1) Solid / (2) Hollow / (3) Hole only
        %   XC = Center x coordinate
        %   YC = Center y coordinate
        %   ZC = Center z coordinate
        %   DX = Mesh Cell size
        %   R = Outer radius
        %   L_CYL = Cylinder length
        %   AXIS = axis of cylinder
        %INDEX=Data(:,1);
        %SOL=Data(:,2);
        XC=Cylinder_Data(:,3);
        YC=Cylinder_Data(:,4);
        %ZC=Data(:,5);
        %DX=Data(:,6);
        RO=Cylinder_Data(:,7);
        %L_CYL=Data(:,8);
        %AXIS=char(Data(:,9));
        
        
        % Figure configuration
        figure(2)
        axis([-50 265 -50 110])
        set(gcf,'position',[40,100,2400,1200],'color','w')
        set(gca,'XTick',-50:5:265)
        set(gca,'YTick',-50:5:110)
        set(gca,'XMinorTick','on')
        set(gca,'YMinorTick','on')
        grid on
        % Ploting OBST domain
        hold on
        for j= 1:length(X_Corner)   % rectangular
            rectangle('Position',[X_Corner(j),Y_Corner(j),X_length(j),Y_length(j)], 'FaceColor','cyan')
            %     text(X_Corner(i)+X_length(i)/2,Y_Corner(i)+Y_length(i)/2,num2str(i))   % Adding text of
        end
        
        for j=1:length(XC)   % cylindeer
            rectangle('Position',[XC(j)-RO(j),YC(j)-RO(j),2*RO(j),2*RO(j)],'Curvature',[1,1],  'FaceColor','cyan')
            %      text(XC(j)-1,YC(j),num2str(j));
        end
        
        for i=1:length(devc)
            for j=1:length(index)
                if strcmp(devc(i).id, char(devc_id(index(j))))
                    plot(devc(i).xyz(1),devc(i).xyz(2),'b*')
                    ttx=strrep(devc(i).id,'_','\_');
                    text(devc(i).xyz(1),devc(i).xyz(2),ttx)
                end
            end
        end
        fig = gcf;
        fig.PaperUnits = 'points';
        fig.PaperPosition = [40,100,2400,1200];
        print(figurename2,'-dpng','-r400')
        %     saveas(gcf,figurename2);
        hold off
        close(2)
    end
end

