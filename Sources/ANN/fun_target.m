
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
