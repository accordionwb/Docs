function is_in=fun_is_in_domain(xq,yq,zq,Cubic_Data,Cylinder_Data,Sphere_Data)
% This function must be used in the directory where OBST_Location.mat is located
% The coordinates of Rectangle must be specified
% Input xq, yq is the point queried 
% The Domain is pre-spedified for FDS simulation

% % Old Spectra Data Field % %
% Rec=[ ...
%     255.0, 262.0, 15.0,  22.0;
%     266.0, 273.0, 15.0,  22.0;
%     250.0, 280.0, 32.0,  40.0;
%     250.0, 285.0, 45.0,  95.0;
%     105.0, 165.0, 65.0,  88.0;
%     20.0,  30.0,  65.0,  70.0;
%     20.0,  30.0,  72.0,  77.0;
%     5.0,   10.0,  65.0,  70.0;
%     60.0,  95.0,  65.0,  80.0;
%     170.0, 200.0, 68.0,  78.0];
% load(cylinder);
% load(rectangle);
% Rec=Rec_Data;

% % Cubic Data Preperation 
XR=[Cubic_Data(:,1),Cubic_Data(:,1),Cubic_Data(:,2),Cubic_Data(:,2),Cubic_Data(:,1)];
YR=[Cubic_Data(:,3),Cubic_Data(:,4),Cubic_Data(:,4),Cubic_Data(:,3),Cubic_Data(:,3)];
ZL=Cubic_Data(:,5); % Z low
ZH=Cubic_Data(:,6); % Z High
m1=size(XR,1);

k=1;
in=0; % Parameter Initialization
on=0; % Parameter Initialization

for i=1:m1  % For all Cubic OBST
    if zq>=ZL(i) && zq<=ZH(i)
    xv=XR(i,:);
    yv=YR(i,:);
    [in(k),on(k)]=inpolygon(xq,yq,xv,yv);
    k=k+1;
    else
        continue
    end
end
% % Cylinder Data Preperation
L0=linspace(0, 2.*pi, 360);
X0=Cylinder_Data(:,3);   % X 
Y0=Cylinder_Data(:,4);   % Y
H0=Cylinder_Data(:,8);   % Height
R0=Cylinder_Data(:,7);   % Raduis

for j=1:length(X0)
    if zq<=H0(j)
    xv=X0(j)+R0(j)*cos(L0)';
    yv=Y0(j)+R0(j)*sin(L0)';
    [in(k),on(k)]=inpolygon(xq,yq,xv,yv);
    k=k+1;
    else 
        continue
    end
end
% % Sphere Data Preperation
L1=linspace(0, 2.*pi, 360);
X1=Sphere_Data(:,3);   % X 
Y1=Sphere_Data(:,4);   % Y
Z1=Sphere_Data(:,5);   % Z
R1=Sphere_Data(:,6);   % Raduis

for l=1:length(X1)
    if Z1(l)>R1(l)
        if zq>=Z1(l)-R1(l) && zq<=Z1(l)+R1(l)
            Rq=sqrt(R1(l)^2-(Z1(l)-zq)^2);
            xv=X1(l)+Rq*cos(L1)';
            yv=Y1(l)+Rq*sin(L1)';
            [in(k),on(k)]=inpolygon(xq,yq,xv,yv);
            k=k+1;
        else
            continue
        end
    else 
        if zq<=Z1(l)+R1(l)
            Rq=sqrt(R1(l)^2-(Z1(l)-zq)^2);
            xv=X1(l)+Rq*cos(L1)';
            yv=Y1(l)+Rq*sin(L1)';
            [in(k),on(k)]=inpolygon(xq,yq,xv,yv);
            k=k+1;
        else
            continue
        end
    end
end
% % Summary


is_in=sum(in)+sum(on);


