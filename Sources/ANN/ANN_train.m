function [net,tr,Input,Target,Test,Check,simout,Info]=ANN_train(Results,Var,dirc)


%% Parameter praparation
filter=Results.filter;
parameter=Results.parameter;
mode=Var.mode;
play_sound=Var.sound;
%%  ANN training process
[Input,Target,Test,Check,Info]=fun_ANN_all_preperation(Results,Var);
Info.type
dim=Var.dim;
if (strcmp(mode, 'concentration'))
    net1=newff(Input,Target(1,:),dim);  % For concentration use only
    net1.trainParam.epochs=Var.epochs;
    net1.divideParam.testRatio=0.10;
    net1.divideParam.trainRatio=0.75;
    net1 = init(net1);
    [net1,tr1]=train(net1,Input,Target(1,:));
    simout(1,:)=sim(net1,Test);
    nntraintool close
    if play_sound
        sound(sin(0:0.1:1000),44100)
    end
    % ---------------------------------------------------------------- %
    net2=newff(Input,Target(2,:),dim);  % For time use only
    net2.trainParam.epochs=Var.epochs;
    net2.divideParam.testRatio=0.10;
    net2.divideParam.trainRatio=0.75;
    net2 = init(net2);
    [net2,tr2]=train(net2,Input,Target(2,:));
    simout(2,:)=sim(net2,Test);
    % ---------------------------------------------------------------- %
    tr.concentrationParam=tr1;
    tr.timeParam=tr2;
    net.net1=net1;
    net.net2=net2;
else
    net=newff(Input,Target,dim);  %## number of neurons in the hidden layer ##%
    net.trainParam.epochs=Var.epochs;
    net.divideParam.testRatio=0.10;
    net.divideParam.trainRatio=0.75;
    net = init(net);
    [net,tr]=train(net,Input,Target);
    simout=sim(net,Test);
end


%% Training Results saving
pwd
file_name=[datestr(now,30),'_network_STA-',filter.stability,'.mat'];
[pathdir,name,ext]=fileparts(file_name);
information={'%%%%  Information of saved network from data preperation process  %%%%';
    ['The file name is: ',file_name]};
output_name=['../Run_Save/',name,'.txt'];
diary(output_name);
diary on;
information
parameter
filter
Var
Info.type
diary off
save(['../Run_Save/',file_name],'net','tr','Input','Target','Test','Check','simout','Info')

if play_sound
    load handel
    sound(y,Fs)
end
end

function [Input,Output,Test,Check,Info]=fun_ANN_all_preperation(Results,Var)
%% Introduction
% This function is used to extract the generated data for neural network training.
% status:
%      0/1/2     'concentration'  /  'both' / 'width'
%      3/0     'pressure'  Y/N
%      6/0     'stability' Y/N

%% initialization
In0=Results.Input;
Out0=Results.Output;

% Variable
mode=Var.mode;
if (Results.filter.pressure==0)
    single_pressure=false;
else
    single_pressure=true;
end

if (strcmp(Results.filter.stability, 'all') || strcmp(Results.filter.stability, 'ALL'))
    single_stability=false;
else
    single_stability=true;
end

% Evaluation
para=[0 0 0];
if (strcmp(mode, 'concentration'))
    para(1)=0;
elseif(strcmp(mode, 'both'))
    para(1)=1;
elseif(strcmp(mode,'width'))
    para(1)=2;
end
if (single_pressure)
    para(2)=3;
else
    para(2)=0;
end
if (single_stability)
    para(3)=6;
else
    para(3)=0;
end
status=sum(para);


%% Siwtch case of status;
switch status
    case 0  % concentration + all_pressure + all_stability
        index_width=find(Out0(1,:)./Out0(2,:)<1);
        index_distance=find(Out0(3,:)./Out0(4,:)<1);
        index_used=intersect(index_width,index_distance);
        In=In0(1:11,index_used);
        Out=Out0(5:6,index_used);
        Info.type=['Status: ',num2str(status),' || Inputs,1:11 and Output, 5:6, ( concentration/time + all_pressure + all_stability )'];
        Info.status=status;
        Info.index1=index_width;
        Info.index2=index_distance;
        Info.index=index_used;
        
    case 1  % both L1/L2 W1/W2 + all_pressure + all_stability
        In=In0(1:11,:);
        Out=Out0(1:4,:);
        Info.type=['Status: ',num2str(status),' || Inputs,1:11 and Output, 1:4, ( L1/L2 & W1/W2 + all_pressure + all_stability )'];
        Info.status=status;
        
    case 2 % W1/W2 + all_pressure + all_stability
        index_distance=find(Out0(3,:)./Out0(4,:)<1);
        index_used=index_distance;
        In=In0(1:11,index_used);
        Out=Out0(1:2,index_used);
        Info.type=['Status: ',num2str(status),' || Inputs,1:11 and Output, 1:2, ( W1/W2 + all_pressure + all_stability )'];
        Info.status=status;
        Info.index2=index_distance;
        Info.index=index_used;
        
    case 3 % concentration + single_pressure + all stability
        index_width=find(Out0(1,:)./Out0(2,:)<1);
        index_distance=find(Out0(3,:)./Out0(4,:)<1);
        index_used=intersect(index_width,index_distance);
        In=In0([1:9 11],index_used);
        Out=Out0(5:6,index_used);
        Info.type=['Status: ',num2str(status),' || Inputs,[1:9 11] and Output, 5:6, ( concentration/time + single_pressure + all_stability )'];
        Info.status=status;
        Info.index1=index_width;
        Info.index2=index_distance;
        Info.index=index_used;
        
    case 4 %both + single_pressure + all_stability
        In=In0([1:9 11],:);
        Out=Out0(1:4,:);
        Info.type=['Status: ',num2str(status),' || Inputs,[1:9 11] and Output, 1:4, ( L1/L2 & W1/W2 + single_pressure + all_stability )'];
        Info.status=status;
        
    case 5 % width (W1/W2) + single_pressure + all_stability
        index_distance=find(Out0(3,:)./Out0(4,:)<1);
        index_used=index_distance;
        In=In0([1:9 11],index_used);
        Out=Out0(1:2,index_used);
        Info.type=['Status: ',num2str(status),' || Inputs,1:11 and Output, 1:2, ( W1/W2 + single_pressure + all_stability )'];
        Info.status=status;
        Info.index2=index_distance;
        Info.index=index_used;
        
    case 6 % concentration + all_pressure + single_stability
        index_width=find(Out0(1,:)./Out0(2,:)<1);
        index_distance=find(Out0(3,:)./Out0(4,:)<1);
        index_used=intersect(index_width,index_distance);
        In=In0(1:10,index_used);
        Out=Out0(5:6,index_used);
        Info.type=['Status: ',num2str(status),' || Inputs,1:10 and Output, 5:6, ( concentration/time + all_pressure + single_stability )'];
        Info.status=status;
        Info.index_width=index_width;
        Info.index_distance=index_distance;
        Info.index=index_used;
        
    case 7  % both + all_pressure + single_stability
        In=In0(1:10,:);
        Out=Out0(1:4,:);
        Info.type=['Status: ',num2str(status),' || Inputs,1:10 and Output, 1:4, ( L1/L2 & W1/W2 + all_pressure + single_stability )'];
        Info.status=status;
        
    case 8  % width + all_pressure + single_stability
        index_distance=find(Out0(3,:)./Out0(4,:)<1);
        index_used=index_distance;
        In=In0(1:10,index_used);
        Out=Out0(1:2,index_used);
        Info.type=['Status: ',num2str(status),' || Inputs,1:10 and Output, 1:2, ( W1/W2 + all_pressure + single_stability )'];
        Info.status=status;
        Info.index_distance=index_distance;
        Info.index=index_used;
        
    case 9 %  concentration + single_pressure + single_stability
        index_width=find(Out0(1,:)./Out0(2,:)<1);
        index_distance=find(Out0(3,:)./Out0(4,:)<1);
        index_used=intersect(index_width,index_distance);
        In=In0(1:9,index_used);
        Out=Out0(5:6,index_used);
        Info.type=['Status: ',num2str(status),' || Inputs,1:9 and Output, 5:6, ( concentration/time + single_pressure + single_stability )'];
        Info.status=status;
        Info.index_width=index_width;
        Info.index_distance=index_distance;
        Info.index=index_used;
        
    case 10 % both + single_pressure + single_stability
        In=In0(1:9,:);
        Out=Out0(1:4,:);
        Info.type=['Status: ',num2str(status),' || Inputs,1:9 and Output, 1:4, ( L1/L2 & W1/W2 + single_pressure + single_stability )'];
        Info.status=status;
        
    case 11 % width + single_pressure + single_stability
        index_distance=find(Out0(3,:)./Out0(4,:)<1);
        index_used=index_distance;
        In=In0(1:9,index_used);
        Out=Out0(1:2,index_used);
        Info.type=['Status: ',num2str(status),' || Inputs,1:9 and Output, 1:2, ( W1/W2 + single_pressure + single_stability )'];
        Info.status=status;
        Info.index_distance=index_distance;
        Info.index_used=index_used;
        
end
Info.stability=Results.filter.stability;
%% Random number of input
if (Var.random)
    Nmax=length(In);
    rtemp=randperm(Nmax);
    Nin1=floor(Nmax*Var.percentage(1));
    Nsn1=Var.percentage(2);
    Nin=max(Nin1,Nsn1);    % number of pairs used
    rdex=rtemp(1:Nin);   % Index of pairs used
    
    r1=Var.train_vs_test(1);r2=Var.train_vs_test(2);
    train_number=floor(Nin*r1/(r1+r2));
    for i=1:train_number
        Input_t(:,i)=In(:,rdex(i));
        Target_t(:,i)=Out(:,rdex(i));
        index_train(i)=rdex(i);
    end
    for i=train_number+1:Nin
        Test_t(:,i-train_number)=In(:,rdex(i));
        Check_t(:,i-train_number)=Out(:,rdex(i));
        index_test(i-train_number)=rdex(i);
    end
    Info.random='yes';
    Info.index_train=index_train;
    Info.index_test=index_test;
    
    Input=Input_t;
    Output=Target_t;
    Test=Test_t;
    Check=Check_t;
end


end


