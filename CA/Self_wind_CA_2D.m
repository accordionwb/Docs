%% Self built method
% D2Q5 neighborhood
%        i,j+1
% i-1,j  i,j    i+1,j
%        i,j-1
%
% Local speed c_0,c_1,c_2,c_3,c_4
%           c2(0,1)
%  c3(-1,0) c0(0,0) c1(1,0)
%           c4(0,-1)
% Weight coefficient w_1,2,3,4
w1=1/4;
w2=1/4;
w3=1/4;
w4=1/4;

% -------- Parameter Initialization ------------
% Wind_Speed=3; %m/s wind speed
% theta=deg2rad(45);
% U0=Wind_Speed*sin(theta);
% V0=Wind_Speed*cos(theta);
U0=3;
V0=-2;

% local speed
c=cell(5,1);
c{1}=[1,0];
c{2}=[0,1];
c{3}=[-1,0];
c{4}=[0,-1];
c{5}=[0,0];

% Domain settings
IX=100;
IY=100;
ii=2:IX-1;
jj=2:IY-1;
ALT=zeros(IX,IY);
ALT(60:80,40:60)=ones(21,21);

SU=U0*ones(IX,IY);
SV=V0*ones(IX,IY);
SU(ii,jj)=0;
SV(ii,jj)=0;
New_SU=SU;
New_SV=SV;

%% Time counting
dt=0.1;
% dx=0.5*Wind_Speed*dt;
t=0;
T=10;
k=0;

while t<=T
    t=t+dt;
    k=k+1;
    
    for i=ii
        for j=jj
    % Local rules
    
%     if ALT(i+1,j)==ALT(i,j) && ALT(i,j+1)==ALT(i,j)
    New_SU(i,j)=SU(i,j)+w1*SU(i+1,j)-w3*SU(i-1,j);    
    New_SV(i,j)=SV(i,j)+w2*SV(i,j+1)-w4*SV(i,j-1);
    
%     elseif ALT(i+1,j)>ALT(i,j)
%         
%          New_SU(i,j)=SU(i,j)-2*w3*SU(i-1,j);    
%     New_SV(i,j)=SV(i,j)+w2*SV(i,j+1)-w4*SV(i,j-1);
%     
%     elseif ALT(i,j+1)>ALT(i,j)
%         
%          New_SU(i,j)=SU(i,j)+w1*SU(i+1,j)-w3*SU(i-1,j);
%     New_SV(i,j)=SV(i,j)-2*w4*SV(i,j-1);
%     
%     elseif ALT(i,j+1)<ALT(i,j)
%         
%          New_SU(i,j)=SU(i,j)+w1*SU(i+1,j)-w3*SU(i-1,j);
%     New_SV(i,j)=SV(i,j)+2*w2*SV(i,j+1);%-w4*SV(i,j-1);
%     
%     elseif ALT(i+1,j)<ALT(i,j)
%         
%   New_SU(i,j)=SU(i,j)+2*w1*SU(i+1,j);%-w3*SU(i-1,j);
%     New_SV(i,j)=SV(i,j)+w2*SV(i,j+1)-w4*SV(i,j-1);
%     end
    
    
    
    
    
%     % Boundary conditions
%     % West Inlet
%     if U0>=0
%         New_SU(1,j)=U0;
%         SU(1,j)=U0;
%     else
%         New_SU(1,j)=SU(1,j)+w1*SU(1+1,j);
%     end
%     New_SV(1,j)=SV(1,j)+w2*SV(1,j+1)-w4*SV(1,j-1);
%     % North
%     if V0>=0
%         New_SV(i,IY)=SV(i,IY)-w4*SV(i,IY-1);
%     else
%         New_SV(i,IY)=V0;
%         SV(i,IY)=V0;
%     end
%     New_SU(i,IY)=SU(i,IY)+w1*SU(i+1,IY)-w3*SU(i-1,IY);
%     
%     % East
%     if U0>=0
%         New_SU(IX,j)=SU(IX,j)-w3*SU(IX-1,j);
%     else
%         New_SU(IX,j)=U0;
%         SU(IX,j)=U0;
%     end
%     
%     New_SV(IX,j)=SV(IX,j)+w2*SV(IX,j+1)-w4*SV(IX,j-1);
%     % South
%     if V0>=0
%         New_SV(i,1)=V0;
%         SV(i,1)=V0;
%     else
%         New_SV(i,1)=SV(i,1)+w2*SV(i,1+1);
%     end
%     New_SU(i,1)=SU(i,1)+w1*SU(i+1,1)-w3*SU(i-1,1);
    
        end
    end
    disp(['Iteration: ',num2str(k,'%5d')])
    SU=New_SU;
    SV=New_SV;
    drawnow
    quiver(SU,SV);
end

