%% Results function
function Results=summary_to_data(Summary,mode,parameter,filter)
N=length(Summary);
M=length(Summary(1).Sensor);
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
        'X2','2st time'; ...
        'X3','2st Con';...
        'X4','3st time';...
        'X5','3st Con';...
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

