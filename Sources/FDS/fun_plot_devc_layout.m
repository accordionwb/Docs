function Devc_Result=fun_plot_devc_layout(Head, Data, IF_PLOT)
% IF_PLOT is boolean.
%% Set path and read files  Static data that are used for domain generation
addpath('D:\SVN\Code\Matlab\FDS\Spectra')
load('MAT_Rec_config.mat');
load('MAT_Cylinder_config.mat');
load('MAT_Devc_config.mat')

%% Read and plot active devc
[M,N]=size(Data);
T=Data(:,1);
used_devc=Data(:,2:N)*1e6;

for i=1:N-1
    N(i)=norm(used_devc(:,i));
    devc_id(i)=Head(i+1);
end
index=find(N>=1);


for i=1:length(index)
    if IF_PLOT
        figure(1)
        plot(T,used_devc(:,index(i)));
        tty=strrep(char(devc_id(index(i))),'_','\_');
        title(['DEVC "',tty,'" Output'])
        ylabel('ppm')
        xlabel('Time /s')
        figurename=['DEVC_Figure/',char(devc_id(index(i))),'.jpeg'];
        saveas(gcf,figurename);
    end
    Devc_Result(i).id=devc_id(index(i));
    Devc_Result(i).time=T;
    Devc_Result(i).con=used_devc(:,index(i));
end


%% Read content and plot domain with Rectangular and cylinder OBST
if IF_PLOT
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
    %axis([0 215 0 60])
    %set(gcf,'position',[00,300,2150,600],'color','w')
    %set(gca,'XTick',0:5:215)
    %set(gca,'YTick',0:5:60)
    set(gca,'XMinorTick','on')
    set(gca,'YMinorTick','on')
    grid on
    % Ploting OBST domain
    hold on
    for i= 1:length(X_Corner)   % rectangular
        rectangle('Position',[X_Corner(i),Y_Corner(i),X_length(i),Y_length(i)], 'FaceColor','cyan')
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
    hold off
end

