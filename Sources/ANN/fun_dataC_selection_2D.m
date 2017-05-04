%% Main function
function Results=fun_dataC_selection_2D(PhastC,mode,parameter,filter)
%------------------------------- Function Description ----------------------------------------------
% This function calculate all the information needed for neural network training from integrated
% PHAST Combined results. The function works in the following algorithm:
%      1. Filter out the post generated scenarios (including wind direction) that couldn't meet the 
%         requirement of "working case" check. That is at least 3 sensor should have readings.
%      2. If "target mode" is selected, the function will calculate downwind concentration and
%         poisition related parameters. No matter whether the target area will be influenced.
%      3. The Results returned is impeccable. It contains all the cases generated from the original 
%        4886 cases. 
% ---------------------------------------Input Requirements ----------------------------------------
% PhastC:
%      Original case classified by category,containing:
%      type | stability | pressure | bole size | wind speed | distance | concentration | half width | times
%
% mode:
%      1. 'default', normal mode means the only filtration is the "working case", return the working
%          case index, Not used for neural network training.
%      2. 'target', target mode means the program will predict the influence towards downwind areas,
%          doing data extracting at distance and skip the cases that are not influenced.
%      3. 'source', source mode return the information of source, including release rate.
% parameter:
%      1. Real wind direction, default -20:2:20
%      2. downwind distance, the sensitive area fron source, default (500, 0)
%      3. Source distance, distance to the release source
%      4. Sensor gap: distance between each sensro
%      5. Sensor included angle: angle between each sensor
%      6. Low bound. fixed starting length from the source, for data selection
% filter:
%      1. stability category using ['all', or the specific category amaong the list. (A-F)]
%      2. pressure [0 / 3 / 11], 0 means all the storage pressure
%      3. wind [0 or other positive value], indicate the lower bound of the wind.
%      4. category [train / test], two different category of original scenario.
%      5. filename, indicating the file name of output file.
% ------------------------------------ End of notes ------------------------------------------------
%% Defination of Input parameters
if (nargin <= 2)
    error('Input error: insufficient input elements, at least parmeter should be included')
elseif (nargin == 3)
    filter.pressure=0; % bar; unused
    filter.wind=0;
    filter.stability='all';
    filter.category='train';
    if (strcmp(mode , 'default'))
        disp(['Using mode: ',mode,', for calculation'])
    elseif(strcmp(mode ,  'target'))
        disp(['Using mode: ',mode,', for calculation'])
    elseif(strcmp(mode , 'source'))
        disp(['Using mode: ',mode,', for calculation'])
    else
        disp(['Unidentified mode: ',mode,', Please check again'])
    end
else
    if (strcmp(mode ,  'default'))
        disp(['Using mode: ',mode,', for calculation'])
    elseif(strcmp(mode , 'target'))
        disp(['Using mode: ',mode,', for calculation'])
    elseif(strcmp(mode , 'source'))
        disp(['Using mode: ',mode,', for calculation'])
    else
        disp(['Unidentified mode: ',mode,', Please check again'])
    end
end

%% Initialization

% ---------------- Filter information ----------------
stability=filter.stability;
pressure=filter.pressure;  % bar  unused
wind=filter.wind;  % windspeed
% ---------------- parmeter information -------------------
%% Data Process and Selection Process Start
n=1; % active case number is initialized

Num_of_Case=length(PhastC);

%%  Scenarios filtration, activate all the scenarios from input filters.

for k=1:Num_of_Case
    % -------------------- Basic information transition ------------------------%
    Summary(k).type=PhastC(k).type;
    Summary(k).stability=PhastC(k).stability;
    Summary(k).pressure=PhastC(k).pressure;
    Summary(k).bole=PhastC(k).bole;
    Summary(k).wind_speed=PhastC(k).wind_speed;
    % -------------------- Selection of data category: --------% 1.  Categroy selection (type: train/test)
    if(strcmp(PhastC(k).type,filter.category) || strcmp('all',filter.category))
        if(strcmp(PhastC(k).type,filter.category))
            disp(['Using type',PhastC(k).type,', at case number:',num2str(k)])
        else
            disp(['Mode: all, type:',PhastC(k).type,' at case number:',num2str(k)])
        end
        Summary(k).filter=boolean([1,0,0,0]);  %[train/test, all/part, ]
    else
        Summary(k).filter=boolean([0,0,0,0]);
        continue  % k increase
    end
    % -------------------- Filtration of atmospheric stability -------------------------- % 2: stability kept as origin
    if (strcmp(stability,'all') || strcmp(stability,'All') )
        all_stability=boolean(true);
        if (all_stability)
            Summary(k).filter=boolean([1,1,0,0]);
            disp('Using all stability category')
        end
    else
        all_stability=boolean(false);
        if (char(PhastC(k).stability) == stability)
            Summary(k).filter=boolean([1,1,0,0]);
            disp(['Using stability ',stability])
        else
            Summary(k).filter=boolean([1,0,0,0]);
            continue  % k increase
        end
    end
    
    
    %----------------------- Data Selection Judgeing Pressure ------------------------%  3 : Pressure
    
    if (pressure == 0)
        all_pressure=boolean(true);
        Summary(k).filter=boolean([1,1,1,0]);
        disp('Using all pressure category')
    else
        all_pressure=boolean(false);
        if(pressure ==PhastC(k).pressure)
            Summary(k).filter=boolean([1,1,1,0]);
            disp(['Use data of pressure ',num2str(PhastC(k).pressure),'bar'])
        else
            Summary(k).filter=boolean([1,1,0,0]);
            disp(['Skip selecting data of pressure ',num2str(PhastC(k).pressure),'bar'])
            continue
        end
    end
    
    % -------------- Selection of wind based on Order category ----------% 4: Wind speed lowbond
    
    if (wind~=0 && wind >= PhastC(k).wind_speed)  % minimium wind speed  filteration
        Summary(k).filter=boolean([1,1,1,0]);  % Delete the wind speed under the threshold.
        disp(['Skip scenario under wind speed',num2str(wind),'m/s'])
        continue
    else
        Summary(k).filter=boolean([1,1,1,1]);   % meaning: all selections have been accomplished
    end
    if (sum(Summary(k).filter)==4)
        Index(n)=k;
        n=n+1;
    end
end
%% Calculation of input and output
direction=parameter.direction;
if (strcmp (mode,'default'))
    working_case=fun_default(Index, PhastC, parameter);   % Calling default function, calculate the working case index and other information.
    Num_of_direction=length(direction);
    for k=1:Num_of_Case
        if(sum(Summary(k).filter) == 4)
            i=find(Index==k);
            for s=1:Num_of_direction
                Summary(k).direction(s)=direction(s);
                Summary(k).Sensor(s).is_working_case=working_case(i).judgement(s);
                Summary(k).Sensor(s).active_sensor=working_case(i).active_sensor(s).acting_num;
                Summary(k).Sensor(s).alarm_time=working_case(i).alarm_time(s);
                Summary(k).Sensor(s).outx=working_case(i).outx(s);
            end
        else  % do not pass the validation check (one or more selection part is denied )
            for s=1:Num_of_direction
                Summary(k).direction(s)=direction(s);
                Summary(k).Sensor(s).is_working_case=false;
                Summary(k).Sensor(s).active_sensor=0;
                Summary(k).Sensor(s).alarm_time=0;
                Summary(k).Sensor(s).outx=0;
            end
        end
    end
     Results=summary_to_data(Summary,mode,parameter,filter);
elseif(strcmp(mode, 'target'))
    working_case=fun_default(Index, PhastC, parameter);
    target_case=fun_target(Index,PhastC,parameter,working_case);
    Num_of_direction=length(direction);
    for k=1:Num_of_Case
        if(sum(Summary(k).filter) == 4)   % creadiable scenario
            i=find(Index==k);
            Summary(k).Target_is_reach=target_case(i).is_reach;
            if (Summary(k).Target_is_reach)
                for s=1:Num_of_direction
                    Summary(k).direction(s)=direction(s);
                    Summary(k).Sensor(s).is_working_case=working_case(i).judgement(s);
                    Summary(k).Sensor(s).active_sensor=working_case(i).active_sensor(s).acting_num;
                    Summary(k).Sensor(s).alarm_time=working_case(i).alarm_time(s);
                    Summary(k).Sensor(s).outx=working_case(i).outx(s).X;
                    Summary(k).Target(s).outx=target_case(i).outx(s).X;
                    Summary(k).Target(s).outy=target_case(i).outy(s).Y;
                end
            else
                for s=1:Num_of_direction
                    Summary(k).direction(s)=direction(s);
                    Summary(k).Sensor(s).is_working_case=working_case(i).judgement(s);
                    Summary(k).Sensor(s).active_sensor=working_case(i).active_sensor(s).acting_num;
                    Summary(k).Sensor(s).alarm_time=working_case(i).alarm_time(s);
                    Summary(k).Sensor(s).outx=working_case(i).outx(s).X;
                    %                     Summary(k).Target(s).outx=0;
                    %                     Summary(k).Target(s).outy=0;
                end
            end
        else                              % Incrediable scenarios
            for s=1:Num_of_direction
                Summary(k).direction(s)=direction(s);
                Summary(k).Sensor(s).is_working_case=false;
                %                 Summary(k).Sensor(s).active_sensor=0;
                %                 Summary(k).Sensor(s).alarm_time=0;
                %                 Summary(k).Sensor(s).outx=0;
                %                 Summary(k).Target(s).outx=0;
                %                 Summary(k).Target(s).outy=0;
                %             end
            end
        end
    end
     Results=summary_to_data(Summary,mode,parameter,filter);
elseif(strcmp(mode,'source'))
    error('Source trace is not finished')     
end

end % End of primary function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%------------ Judging whether the boolean "Preditiction" would work -----%       % 4: Is Predition?  Default/Source/Target
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



%% Target function
function target_case=fun_target(Index,PhastC,parameter,working_case)

%% Initialization of Input data
% ----------------- Parameter information --------------------
direction=parameter.direction;% direction
% downwind area location
downwindx=parameter.downwindx;
downwindy=parameter.downwindy;
target_angle=atan(downwindy./downwindx);
downwind=sqrt(downwindx.^2+downwindy.^2);
% sensor parameters
% output file name
% fixed working range from source
target_lowbond=min(downwind)-min([50,min(downwind)*0.1]);
target_upbond=max(downwind)+max([50,max(downwind)*0.1]); % up bound
stability_list=cell({'A';'AB';'B';'BC';'C';'CD';'D';'E';'F'});
%% ------------ Preperation of target outputs -----%       %  Target mode
N=length(Index);
for i=1:N
    Distance0=PhastC(Index(i)).Distance;
    Concentration0=PhastC(Index(i)).Concentration;
    Width0=PhastC(Index(i)).Width;
    Times0=PhastC(Index(i)).Times;
    Lmax=max(Distance0);
    % preprocessing of initial Data.
    dend=length(Distance0);
    dstart=find((Distance0<=target_lowbond*0.3 & Distance0>0),1,'last');
    DD0=Distance0(dstart);DDN=Distance0(dend);
    Distance=linspace(DD0,DDN,dend*10);   % 10 times encryption
    Concentration=spline(Distance0(dstart:dend),Concentration0(dstart:dend),Distance);
    Width=spline(Distance0(dstart:dend),Width0(dstart:dend),Distance);
    Times=spline(Distance0(dstart:dend),Times0(dstart:dend),Distance);
    
    % ---------------- Evaluation of interp x based on distance category -----------  % target mode
    if (Lmax>target_lowbond)
        starts=find((Distance<=target_lowbond & Distance>0),1,'last');
        
        if (Lmax<=target_upbond)
            ends=find((Distance<=Lmax & Distance>target_lowbond),1,'last');
            x=[Distance(starts),target_lowbond:0.5:floor(Distance(ends)),Distance(ends)];
        else
            ends=find(Distance>=target_upbond,1,'first');
            x=[Distance(starts),target_lowbond:0.5:round(Distance(ends))];
        end
        target_case(i).is_reach=true;
    else
        target_case(i).is_reach=false;
        continue;
    end
    % ----------------------  Data intergration --------------------------------
    x0=Distance(starts:ends);
    
    c0=Concentration(starts:ends);
    t0=Times(starts:ends);
    w0=Width(starts:ends);
    if (length(c0)<2)
        disp(['Length of c0 is: ',num2str(length(c0))]);
    end
    c=spline(x0,c0,x);
    t=spline(x0,t0,x);
    w=spline(x0,w0,x);  % half width of plume
    
    nt=1;  % target case index
    
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    Num_of_direction=length(direction);
    Num_of_downwind=length(downwind);
    for s=1:Num_of_direction      % Wind Direction is s
        if(working_case(i).judgement(s))
            angle_wind=deg2rad(direction(s));  % Rad angle of wind
            tar_wind_angle=target_angle-angle_wind;
            
            %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
            nout=1;
            for j=1:Num_of_downwind
                L_target=downwind(j)*cos(tar_wind_angle(j));
                Xo(6)=downwindx(j);
                Xo(7)=downwindy(j);
                Xo(8)=PhastC(Index(i)).wind_speed;   % wind speed
                Xo(9)=direction(s);   % wind direction
                Xo(10)=PhastC(Index(i)).pressure;
                Xo(11)=find(strcmp(stability_list, PhastC(Index(i)).stability));  %stability
                if (L_target <= Lmax)
                    iny1= find(x>=downwind(j)*cos(tar_wind_angle(j)),1,'first');
                    Y(5)=c(iny1);                                  % 5 is concentration
                    Y(6)=t(iny1)-working_case(i).alarm_time(s);    % 6 is time
                    Y(2)=max([w(iny1),0.01]);                                  % 2: plume half width
                else
                    Y(5)=0;
                    Y(6)=7200;
                    Y(2)=0.01;
                end
                Y(1)=downwind(j)*sin(abs(tar_wind_angle(j)));     % 1 Target position width
                Y(3)=L_target;                                    % 3 target distance
                Y(4)=Lmax;                                        % 4 Maximum distance
                % Result store
                X_down(:,nout)=Xo';
                Y_down(:,nout)=Y';
                nout=nout+1;
            end
            %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
            target_case(i).outx(s).X=X_down;
            target_case(i).outy(s).Y=Y_down;
        else
            target_case(i).outx(s).X=zeros(6,Num_of_downwind);
            target_case(i).outy(s).Y=zeros(6,Num_of_downwind);
        end  % end if (is working case)
        
    end  % end for number of direction
end  % end for Number of selected scenarios (Index)
end  % End function fun_target

%% Results function
function Results=summary_to_data(Summary,mode,parameter,filter)
N=length(Summary);
M=length(Summary(1).Sensor);
% Define target mode output 
if (strcmp(mode,'target'))
    filename=[datestr(now,30),'_DataC_target_',parameter.filename,'_STA-',filter.stability,'.mat'];
    W=length(parameter.downwindx);
    n=1;
    for k=1:N
        if(sum(Summary(k).filter) == 4)
            if (Summary(k).Target_is_reach)
                for s=1:M
                    if(Summary(k).Sensor(s).is_working_case)
                        for j=1:W
                            Input(1:5,n)=Summary(k).Sensor(s).outx;
                            Input(6:11,n)=Summary(k).Target(s).outx(6:11,j);
                            Output(:,n)=Summary(k).Target(s).outy(:,j);
                            n=n+1;
                        end   % end for j=1:W
                    end % End if  is working case
                end  % end for directions
            else
                continue;
            end
        end  % end if filtrated
    end  % End for all cases
    Results.Summary=Summary;
    Results.Input=Input;
    Results.Output=Output;
    Results.info=['Target mode, total ',num2str(n),' cases.'];
    message1=cell({'X1','1st Con'; ...
        'X2','2nd time'; ...
        'X3','2nd Con';...
        'X4','3rd time';...
        'X5','3rd Con';...
        'X6','Down_X';...
        'X7','Down_Y';...
        'X8','Wind Speed';...
        'X9','Wind direction';...
        'X10','Pressure';...
        'X11','Stability'});
    message2=cell({'Y1','tar_position_width'; ...
        'Y2','half width'; ...
        'Y3','L_target';...
        'Y4','L_max';...
        'Y5','Concentration';...
        'Y6','Time'});
    Results.Input_structure=message1;
    Results.Output_structure=message2;
    Results.saved_file=filename;
    Results.filter=filter;
    Results.parameter=parameter;
    [pathdir,name,ext]=fileparts(filename);
    information={'Information of saved Results from data preperation process';
        ['The file name is: ',filename]};
    output_name=['../Run_Save/',name,'.txt'];
    diary(output_name);
    diary on;
    information
    parameter
    filter
    diary off;
    save(['../Run_Save/',filename],'Results')
else
    direction=parameter.direction;
    filename=[datestr(now,30),'DataC_default_',parameter.filename,'_STA-',filter.stability,'.mat'];
    for s=1:M
        Num_of_working_case=0;
        Eff_case=0;
        Optimal.Sensor(s).direction=direction(s);
        for k=1:N
            if(sum(Summary(k).filter) == 4)
                Num_of_working_case=Num_of_working_case+Summary(k).Sensor(s).is_working_case;
                Eff_case=Eff_case+1;
            end
        end
        Optimal.Sensor(s).working_ratio=Num_of_working_case/Eff_case*100;
    end
    Results.info=['Default mode, total ',num2str(Eff_case),' cases.'];
    Results.Optimal=Optimal;
    Results.Summary=Summary;
    Results.saved_file=filename;
    Results.filter=filter;
    Results.parameter=parameter;

    [pathdir,name,ext]=fileparts(filename);
    information={'Information of saved Results from data preperation process';
        ['The file name is: ',filename]};
    output_name=['../Run_Save/',name,'.txt'];
    diary(output_name);
    diary on;
    information
    parameter
    filter
    diary off;
    save(['../Results/',filename],'Results')
end
end



