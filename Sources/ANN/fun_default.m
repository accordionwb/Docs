%% default function
function working_case=fun_default(Index, PhastC, parameter)
%% Initialization of Input data
% ----------------- Parameter information --------------------
direction=parameter.direction;% direction
Source_distance=parameter.Source_distance;
D0=Source_distance;% source distance
% sensor parameters
gap=parameter.gap;  % fixed value
sangle=parameter.sangle;   % fixed angle between sensors
% output file name
% fixed working range from source
default_lowbond=5;   % fixed value
default_upbond=200;

%% Select case
N=length(Index);
%------------ Judging whether the boolean "Preditiction" would work -----%       % 4: Is Predition?  Normal/Source/Target
for i=1:N
    Distance=PhastC(Index(i)).Distance;
    Concentration=PhastC(Index(i)).Concentration;
    Width=PhastC(Index(i)).Width;
    Times=PhastC(Index(i)).Times;
    
    
    %  ---------------- Evaluation of interp x based on distance category -----------  % default mode
    
    starts=find((Distance<=default_lowbond & Distance>0),1,'last');
    ends=find(Distance>=default_upbond,1,'first');
    x=[Distance(starts),default_lowbond:0.5:default_upbond,Distance(ends)];
    
    % ----------------------  Data intergration --------------------------------
    x0=Distance(starts:ends);
    c0=Concentration(starts:ends);
    t0=Times(starts:ends);
    w0=Width(starts:ends);
    c=spline(x0,c0,x);
    t=spline(x0,t0,x);
    w=spline(x0,w0,x);  % half width of plume
    angle_smoke=atan(w./x);  % Unit in rad
    
    % Part of Input angel calculation
    angle_sensor1=0;    % Unit is rad  !!!!!!!!!
    angle_sensor2=atan(gap*sin(sangle)/(gap*cos(sangle)+D0));
    angle_sensor3=0;
    angle_sensor4=-atan(gap*sin(sangle)/(gap*cos(sangle)+D0));
    %%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Num_of_direction=length(direction);
    for s=1:Num_of_direction      % Wind Direction is s
        angle_wind=deg2rad(direction(s));  % Rad angle of wind
        angle_sen_wind(1)=abs(angle_wind-angle_sensor1);  % included angel of wind and sensor
        angle_sen_wind(2)=abs(angle_wind-angle_sensor2);
        angle_sen_wind(3)=abs(angle_wind-angle_sensor3);
        angle_sen_wind(4)=abs(angle_wind-angle_sensor4);
        %----------------------------------------------------------------------------------------------------------------%
        sen_wind_pro(1)=D0/cos(angle_sensor1)*cos(angle_sen_wind(1));                  % Projection on wind direction
        sen_wind_pro(2)=(D0+gap*cos(sangle))/cos(angle_sensor2)*cos(angle_sen_wind(2));                                  %
        sen_wind_pro(3)=(D0+2*gap*cos(sangle))/cos(angle_sensor3)*cos(angle_sen_wind(3));                                %
        sen_wind_pro(4)=(D0+gap*cos(sangle))/cos(angle_sensor4)*cos(angle_sen_wind(4));                                  %
        %----------------------------------------------------------------------------------------------------------------%
        
        % Condition classification part: Judging whether three sensors could gain record
        acting_num=boolean([0 0 0 0]);
        
        %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
        
        if (angle_wind>0)  %% First condition, wind_angle >0
            inx1=find(x>=sen_wind_pro(1),1,'first');
            inx2=find(x>=sen_wind_pro(2),1,'first');
            inx3=find(x>=sen_wind_pro(3),1,'first');
            inx4=find(x>=sen_wind_pro(4),1,'first');
            
            if (angle_sen_wind(1)<angle_smoke(inx1))
                X(1)=c(inx1);T0=t(inx1);
                acting_num(1)=true;
            else
                X(1)=0.01; T0=t(inx1);   %first Sensor
            end
            
            if (angle_sen_wind(2)< angle_smoke(inx2))
                X(2)=t(inx2)-T0;X(3)=c(inx2); % Second Sensor
                acting_num(2)=true;
            else
                X(2)=7200;X(3)=0.01; % Second Sensor
            end
            
            if (angle_sen_wind(3)< angle_smoke(inx3))
                X(4)=t(inx3)-T0;X(5)=c(inx3); % Third Sensor
                acting_num(3)=true;
            else
                X(4)=7200;X(5)=0.01; % Third Sensor
            end
            
            if (angle_sen_wind(4)< angle_smoke(inx4))
                X(6)=t(inx4)-T0;X(7)=c(inx4); % Forth Sensor
                acting_num(4)=true;
            else
                X(6)=7200;X(7)=0.01; % Forth Sensor
            end
            %------------------------------------------------------------------
            
            
            Xo(1)=X(1); % 1st con
            Xo(2)=X(2); % 2rd time
            Xo(3)=X(3); % 2rd con
            Xo(4)=X(4); % 3th time
            Xo(5)=X(5); % 3th con
            
            %-------------------------------------------------------------------------------------%
        elseif (angle_wind<0)   %% Condition 2, Wind angle <0
            
            inx1=find(x>=sen_wind_pro(1),1,'first');
            inx2=find(x>=sen_wind_pro(2),1,'first');
            inx3=find(x>=sen_wind_pro(3),1,'first');
            inx4=find(x>=sen_wind_pro(4),1,'first');
            if (angle_sen_wind(1)<angle_smoke(inx1))
                X(1)=c(inx1);T0=t(inx1);
                acting_num(1)=true;
            else
                X(1)=0.01; T0=t(inx1);   %first Sensor
            end
            if (angle_sen_wind(2)< angle_smoke(inx2))
                X(2)=t(inx2)-T0;X(3)=c(inx2); % Second Sensor
                acting_num(2)=true;
            else
                X(2)=7200;X(3)=0.01; % Second Sensor
            end
            if (angle_sen_wind(3)< angle_smoke(inx3))
                X(4)=t(inx3)-T0;X(5)=c(inx3); % Third Sensor
                acting_num(3)=true;
            else
                X(4)=7200;X(5)=0.01; % Third Sensor
            end
            if (angle_sen_wind(4)< angle_smoke(inx4))
                X(6)=t(inx4)-T0;X(7)=c(inx4); % Forth Sensor
                acting_num(4)=true;
            else
                X(6)=7200;X(7)=0.01; % Forth Sensor
            end
            
            %%%%%%%%%%%%%% Judeing of pressure category %%%%%%%%%%%%%%%
            
            Xo(1)=X(1); % first
            Xo(2)=X(6); % forth
            Xo(3)=X(7); % forth
            Xo(4)=X(4); % third
            Xo(5)=X(5); % thifd
            
            %----------------------------------------------------------------------------------%
        else   % wind_angle =0
            inx1=find(x>=sen_wind_pro(1),1,'first');
            inx2=find(x>=sen_wind_pro(2),1,'first');
            inx3=find(x>=sen_wind_pro(3),1,'first');
            inx4=find(x>=sen_wind_pro(4),1,'first');
            if (angle_sen_wind(1)<angle_smoke(inx1))
                X(1)=c(inx1);T0=t(inx1);
                acting_num(1)=true;
            else
                X(1)=0.01; T0=t(inx1);   %first Sensor
            end
            if (angle_sen_wind(2)< angle_smoke(inx2))
                X(2)=t(inx2)-T0;X(3)=c(inx2); % Second Sensor
                acting_num(2)=true;
            else
                X(2)=7200;X(3)=0.01; % Second Sensor
            end
            if (angle_sen_wind(3)< angle_smoke(inx3))
                X(4)=t(inx3)-T0;X(5)=c(inx3); % Third Sensor
                acting_num(3)=true;
            else
                X(4)=7200;X(5)=0.01; % Third Sensor
            end
            if (angle_sen_wind(4)< angle_smoke(inx4))
                X(6)=t(inx4)-T0;X(7)=c(inx4); % Forth Sensor
                acting_num(4)=true;
            else
                X(6)=7200;X(7)=0.01; % Forth Sensor
            end
            %------------------------------------------------------------------
            %%%%%%%%%%%%%% Judeing of Pressure %%%%%%%%%%%%%%
            
            Xo(1)=X(1); % first
            Xo(2)=X(2); % second
            Xo(3)=X(3); % second
            Xo(4)=X(4); % third
            Xo(5)=X(5); % third
            
        end   % Done  % Output is Xo with 5 parameters
        
        %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
        
        if (sum(acting_num)>=3)   % i is the Index
            working_case(i).judgement(s)=true;
            working_case(i).active_sensor(s).acting_num=acting_num;
            working_case(i).outx(s).X=Xo';
            working_case(i).alarm_time(s)=T0;
        else
            working_case(i).judgement(s)=false;
            working_case(i).active_sensor(s).acting_num=acting_num;
            working_case(i).outx(s).X=zeros(5,1);
            working_case(i).alarm_time(s)=0;
            disp(['Case: ',num2str(Index(i)),';  skip wind direction: ',num2str(direction(s))])
            continue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end  % End of direction selection of default mode % End of s
end   % End of i from Index
end