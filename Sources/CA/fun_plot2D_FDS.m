function fun_plot2D_FDS(con,time,dim,mycmap,chid,layer,arch)
%% This script plot all domain OBST using matlab function
% Basic function for the plot is rectangle
% For circle: rectangle('Position',[Bottom_left_x,Bottom_left_y,Diameter_x,Diameter_y],'Curvature',[1,1],  'FaceColor','cyan')
% For Rectangle: rectangle('Position',[Bottom_left_x,Bottom_left_y,Length_x,Length_y], 'FaceColor','cyan')
% axis equal is required

% Add domain configuration file to path

if strcmp(arch,'linux')
    addpath('/home/wangbing/SVN/Code/Matlab/FDS/Scenarios/Facility')
    outpath=['/home/wangbing/Videos/',chid];
elseif strcmp(arch,'win')
    addpath('C:\SVN/Code\Matlab\FDS\Scenarios\Facility')
    outpath=['E:\Videos\',chid];
elseif strcmp(arch,'mac')
    
end
% addpath('~/fdscov')
% fid='029';

if ~exist('Rec_Data','var')
    load('MAT_Rec_config.mat');
    load('MAT_Cylinder_config.mat');
end


% open video file

if ~exist(outpath,'dir')
    mkdir(outpath)
end

Tdefault=1000;

% if ~exist('con','var')
%     load(['~/fdscov/Spectra_',fid,'_con.mat'],'con');
% end

%% Adjustable Control variables
dike_flag=0;  % 0 -> no dike, 1 -> dike on
outer_flag=1;

%% Plot time serials 2D field data
IT=length(time);
X_cor=[-60,275];
Y_cor=[-60,120];

% Plot domain base map
if dike_flag==0
    tmp=Rec_Data(1:10,:);
    Rec_Data=tmp;
end

% Rec Data
X_Corner=Rec_Data(:,1);
Y_Corner=Rec_Data(:,3);
X_length=Rec_Data(:,2)-Rec_Data(:,1);
Y_length=Rec_Data(:,4)-Rec_Data(:,3);

% Cylinder Data reading
Data=Cylinder_Data;
XC=Data(:,3);
YC=Data(:,4);
RO=Data(:,7);

% Check input
if dim==2
    data=con;
    fname=[outpath,'/',chid,'_CA.avi'];
    vid=VideoWriter(fname);
    open(vid)
elseif dim==3
    [D1,D2,D3,D4]=size(con);
    
    data_tmp=con(:,:,layer,:);
    data=reshape(data_tmp,D1,D2,D4);
    fname=[outpath,'/',chid,'.avi'];
    vid=VideoWriter(fname);
    open(vid)
else
    error('Error input argument')
end

% Make movie equivlent to each other by setting total time step compariable
if IT> 1.1*Tdefault
    skip_length=floor(IT/Tdefault);
    skip_count=0;
else
    skip_length=0;
    skip_count=0;
end

% Plot Loop
for k=1:IT-1
    
    skip_count=skip_count+1;
    if k==1
        skp_tmp=skip_length;
        skip_length=1000;
    else
        skip_length=skp_tmp;
    end
    
    
    if  k==1 || skip_count>=skip_length
        
        imagesc(X_cor,Y_cor,data(:,:,k)');
        ax=gca;
        colormap(ax,mycmap)        
        colorbar
        
        if outer_flag==0
            axis([0 215 0 60])
            set(gcf,'position',[0,0,1280,625],'color','w')
            set(gca,'XTick',0:5:215)
            set(gca,'YTick',0:5:60)
        elseif outer_flag==1
            axis([X_cor,Y_cor])
            set(gcf,'position',[0,0,1280,625],'color','w')
            set(gca,'XTick',X_cor(1):10:X_cor(2))
            set(gca,'YTick',Y_cor(1):10:Y_cor(2))
        end
        
        
        
        caxis([0,0.005])
        title(['Time elapse: ',num2str(round(time(k))),' s'])
        colormap(jet);
        set(gca,'YDir','normal');
        
        
        set(gca,'XMinorTick','on')
        set(gca,'YMinorTick','on')
        grid on
        
        
        % Ploting OBST domain
        hold on
        for i= 1:length(X_Corner)
            rectangle('Position',[X_Corner(i),Y_Corner(i),X_length(i),Y_length(i)], 'FaceColor','cyan')
            if i<=10
                text(X_Corner(i)+X_length(i)/2,Y_Corner(i)+Y_length(i)/2,num2str(i),'FontSize',12)   % Adding text of
            end
        end
        
        % Cylindical OBST
        for j=1:length(XC)
            rectangle('Position',[XC(j)-RO(j),YC(j)-RO(j),2*RO(j),2*RO(j)],'Curvature',[1,1],  'FaceColor','cyan')
            text(XC(j)-1,YC(j),num2str(j),'FontSize',12);
        end
        
        % Tick and Label
        set(gca,'FontSize',12);
        xlabel('X length (m)','FontSize',12)
        ylabel('Y length (m)','FontSize',12)
        
        drawnow
        
        %%% End plot 2D field data
        
        A=getframe(gcf);
        disp(['Frame size at iteration <',num2str(k),'> is: ',num2str(size(A.cdata))])
        writeVideo(vid,A);
        hold off
        
        % reset skip_count
        
        skip_count=0;
        
        
    end
    
    
    
end
% textprogressbar('  Movie generated!');
close(vid)
