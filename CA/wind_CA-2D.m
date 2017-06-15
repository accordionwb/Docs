%% Introduction
% This is a script for 2D wind flow using CA
% Modeling variables and parameters
%
% C(i,j) Cellular automata index
% ALT(i,j)  Altitude (Currently 0 1)
% Type(i,j)  Wind type. Riptide, NOrmal, Double, Triple
% Type_r(i,j)  Wind type received from neighborhood
% Speed(i,j)   Discreated wind velocity: Level 0,1,2,3
% Speed_r(i,j)  Wind speed received from neighboring
% alpha(i,j) index of wind flow (seperation numbers)
%
% Wind direction E-N-W-S : 1-2-3-4
% Defination of Obst/Topography ALT=0,1 in 2D

%% Domain defination
IX=100;
IY=100;
ALT=zeros(IX,IY);
ALT(60:80,40:60)=ones(21,21);
TYPE={'Riptide','Normal','Double','Triple'};
Type_r=cell(IX,IY);
New_Type=cell(IX,IY);
SPL=[0,1,2,3];  %Speed level with actual wind value
New_Speed=zeros(IX,IY);
% Example of vertical wind profile
U0=10; % m/s at 100 m
H=100;
z=0:100;
uz=U0*exp((z-H)/H);
plot(uz,z);

% Initial values
alpha=ones(IX,IY);
Ps=0;
Speed=SPL(2)*ones(IX,IY);
Speed_r=zeros(IX,IY);

% Defining alpha
for i=2:IX-1
    for j=2:IY-1
        
        if ALT(i,j+1)>=ALT(i,j) && ALT(i,j-1)>=ALT(i,j)
            alpha(i,j)=1;
        elseif ALT(i,j+1)>=ALT(i,j) || ALT(i,j-1)>=ALT(i,j)
            alpha(i,j)=2;
        elseif ALT(i,j+1)<ALT(i,j) && ALT(i,j-1)<ALT(i,j)
            alpha(i,j)=3;
        end
    end
end

%% Start loop
dt=0.1;
T=0.5;
t=0;
k=0;
while t<T
    k=k+1;
    t=t+dt;

for i=2:IX-1
    for j=2:IY-1
        
        % Exchanging wind between cells
        % East
        if ALT(i+1,j)>=ALT(i,j)
            Type_r{i,j}=TYPE{2};
            Speed_r(i+1,j)=Speed(i,j)/alpha(i,j);
        else
            Type_r{i+1,j}=TYPE{1};
            Speed_r(i+1,j)=Ps;
        end
        % North
        if ALT(i,j+1)>=ALT(i,j)
            Type_r{i,j+1}=TYPE{2};
            Speed_r(i,j+1)=Speed(i,j)/alpha(i,j);
        else
            Type_r{i,j+1}=TYPE{1};
            Speed_r(i,j+1)=Ps;
        end
        % South
        if ALT(i,j-1)>=ALT(i,j)
            Type_r{i,j-1}=TYPE{2};
            Speed_r(i,j-1)=Speed(i,j)/alpha(i,j);
        else
            Type_r{i,j-1}=TYPE{1};
            Speed_r(i,j-1)=Ps;
        end
        % West
        if ALT(i-1,j)>=ALT(i,j)
            Type_r{i-1,j}=TYPE{2};
            Speed_r(i-1,j)=Speed(i,j)/alpha(i,j);
        else
            Type_r{i-1,j}=TYPE{1};
            Speed_r(i-1,j)=Ps;
        end
        
        % Own cell speed
        New_Type{i,j}=TYPE{alpha(i,j)+1};
        New_Speed(i,j)=Speed_r(i+1,j)/alpha(i+1,j)+Speed_r(i,j+1)/alpha(i,j+1)+Speed_r(i,j-1)/alpha(i,j-1);
        
        % Renew values
        
    end
end
Type=New_Type;
Speed=New_Speed;
end

