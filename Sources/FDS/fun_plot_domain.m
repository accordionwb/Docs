function fun_plot_domain(config,Cubic_Data,Cylinder_Data,Sphere_Data)
% This function plot horizontal domain layout with OBST configurationd data
% config containts domain boundary setup

Origin=config.origin;
XX=config.X_length;
YY=config.Y_length;
major_label_step=config.m_label_step;

fontsize=12;



% Reading Cubic Data
X_Corner=Cubic_Data(:,1);
Y_Corner=Cubic_Data(:,3);
X_length=Cubic_Data(:,2)-Cubic_Data(:,1);
Y_length=Cubic_Data(:,4)-Cubic_Data(:,3);

% Reading Cylinder Data reading
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
XC_cylinder=Cylinder_Data(:,3);
YC_cylinder=Cylinder_Data(:,4);
%ZC=Data(:,5);
%DX=Data(:,6);
RO_cylinder=Cylinder_Data(:,7);
%L_CYL=Data(:,8);
%AXIS=char(Data(:,9));


% Reading Sphere Data
XC_sphere=Sphere_Data(:,2);
YC_sphere=Sphere_Data(:,3);
RO_sphere=Sphere_Data(:,6);


% % Domain Figure
axis([Origin(1),XX,Origin(2),YY])
set(gcf,'position',[20,20,20+XX,20+YY],'color','w')
set(gca,'XTick',Origin(1):major_label_step:XX)
set(gca,'YTick',Origin(2):major_label_step:YY)

% axis([-50 265 -50 110])
% set(gcf,'position',[40,100,2400,1200],'color','w')
% set(gca,'XTick',-50:5:265)
% set(gca,'YTick',-50:5:110)


set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
grid on



% Ploting OBST domain
hold on
for i= 1:length(X_Corner)
    rectangle('Position',[X_Corner(i),Y_Corner(i),X_length(i),Y_length(i)], 'FaceColor','cyan')
    text(X_Corner(i)+X_length(i)/2,Y_Corner(i)+Y_length(i)/2,num2str(i),'FontSize',fontsize)   % Adding text of
   
end

for j=1:length(XC_cylinder)
    rectangle('Position',[XC_cylinder(j)-RO_cylinder(j),YC_cylinder(j)-RO_cylinder(j),2*RO_cylinder(j),2*RO_cylinder(j)],'Curvature',[1,1],  'FaceColor','cyan')
    text(XC_cylinder(j)-1,YC_cylinder(j),num2str(j),'FontSize',fontsize);
end

for k=1:length(XC_sphere)
    rectangle('Position',[XC_sphere(k)-RO_sphere(k),YC_sphere(k)-RO_sphere(k),2*RO_sphere(k),2*RO_sphere(k)],'Curvature',[1,1],  'FaceColor','cyan')
    text(XC_sphere(k)-1,YC_sphere(k),num2str(k),'FontSize',fontsize);
end

% Tick and Label
set(gca,'FontSize',fontsize);
xlabel('X×ø±ê (m)','FontSize',fontsize)
ylabel('Y×ø±ê (m)','FontSize',fontsize)
% 9 Zone configuration
% plot([-50,265],[0,0],'k-',[-50,265],[60,60],'k-',[0,0],[-50,110],'k-',[215,215],[-50,110],'k-')
% text(-25,-25,'1','FontSize',24);
% text(-25,30,'2','FontSize',24);
% text(-25,90,'3','FontSize',24);
% text(115,90,'4','FontSize',24);
% text(240,90,'5','FontSize',24);
% text(240,30,'6','FontSize',24);
% text(240,-25,'7','FontSize',24);
% text(115,-25,'8','FontSize',24);
% text(115,10,'9','FontSize',24);
% 
% 
% % % Source XB configuration
% rectangle('Position',[34.9, 19.5, 0.5, 2], 'FaceColor','r')
% text(35.5,20,'A','FontSize',20,'Color','r');
% rectangle('Position',[91, 24.9, 2, 0.5], 'FaceColor','r')
% text(93.5,26,'B','FontSize',20,'Color','r');
% rectangle('Position',[161.9, 25.5, 0.5, 2], 'FaceColor','r')
% text(162.5,27,'C','FontSize',20,'Color','r');

    
    
% load('MAT_Devc_config.mat');
% for i=1:length(devc)
%     plot(devc(i).xyz(1),devc(i).xyz(2),'b*')
%     ttx=strrep(devc(i).id,'_','\_');
%     text(devc(i).xyz(1),devc(i).xyz(2),ttx)
% end