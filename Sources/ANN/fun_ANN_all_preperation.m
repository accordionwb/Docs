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
        index1=find(Out0(1,:)./Out0(2,:)<1);
        index2=find(Out0(3,:)./Out0(4,:)<1);
        index=intersect(index1,index2);
        In=In0(1:11,index);
        Out=Out0(5:6,index);
        Info.type='Inputs,1:11 and Output, 5:6, ( concentration + all_pressure + all_stability )';
        Info.index1=index1;
        Info.index2=index2;
        Info.index=index;
        
    case 1  % both L1/L2 W1/W2 + all_pressure + all_stability
        In=In0(1:11,:);
        Out=Out0(1:4,:);
        Info.type='Inputs,1:11 and Output, 1:4, ( L1/L2 & W1/W2 + all_pressure + all_stability )';
    case 2 % W1/W2 + all_pressure + all_stability
        index2=find(Out0(3,:)./Out0(4,:)<1);
        index=index2;        
        In=In0(1:11,index);
        Out=Out0(1:2,index);        
        Info.type='Inputs,1:11 and Output, 1:2, ( W1/W2 + all_pressure + all_stability )';
        Info.index2=index2;
        Info.index=index;
        
    case 3 % concentration + single_pressure + all stability
        index1=find(Out0(1,:)./Out0(2,:)<1);
        index2=find(Out0(3,:)./Out0(4,:)<1);
        index=intersect(index1,index2);
        
        In=In0([1:9 11],index);
        Out=Out0(5:6,index);
        
        Info.type='Inputs,[1:9 11] and Output, 5:6, ( concentration + single_pressure + all_stability )';
        Info.index1=index1;
        Info.index2=index2;
        Info.index=index;
        
    case 4 %both + single_pressure + all_stability
        In=In0([1:9 11],:);
        Out=Out0(1:4,:);        
        Info.type='Inputs,[1:9 11] and Output, 1:4, ( L1/L2 & W1/W2 + single_pressure + all_stability )';
        
    case 5 % width (W1/W2) + single_pressure + all_stability
        index2=find(Out0(3,:)./Out0(4,:)<1);
        index=index2;        
        In=In0([1:9 11],index);
        Out=Out0(1:2,index);       
        Info.type='Inputs,1:11 and Output, 1:2, ( W1/W2 + single_pressure + all_stability )';
        Info.index2=index2;
        Info.index=index;
        
    case 6 % concentration + all_pressure + single_stability
        index1=find(Out0(1,:)./Out0(2,:)<1);
        index2=find(Out0(3,:)./Out0(4,:)<1);
        index=intersect(index1,index2);       
        In=In0(1:10,index);
        Out=Out0(5:6,index);        
        Info.type='Inputs,1:10 and Output, 5:6, ( concentration + all_pressure + single_stability )';
        Info.index1=index1;
        Info.index2=index2;
        Info.index=index;
        
    case 7  % both + all_pressure + single_stability
        In=In0(1:10,:);
        Out=Out0(1:4,:);       
        Info.type='Inputs,1:10 and Output, 1:4, ( L1/L2 & W1/W2 + all_pressure + single_stability )';
        
    case 8  % width + all_pressure + single_stability
        index2=find(Out0(3,:)./Out0(4,:)<1);
        index=index2;       
        In=In0(1:10,index);
        Out=Out0(1:2,index);       
        Info.type='Inputs,1:10 and Output, 1:2, ( W1/W2 + all_pressure + single_stability )';
        Info.index2=index2;
        Info.index=index;
       
    case 9 %  concentration + single_pressure + single_stability
        index1=find(Out0(1,:)./Out0(2,:)<1);
        index2=find(Out0(3,:)./Out0(4,:)<1);
        index=intersect(index1,index2);        
        In=In0(1:9,index);
        Out=Out0(5:6,index);        
        Info.type='Inputs,1:9 and Output, 5:6, ( concentration + single_pressure + single_stability )';
        Info.index1=index1;
        Info.index2=index2;
        Info.index=index;
        
    case 10 % both + single_pressure + single_stability
        In=In0(1:9,:);
        Out=Out0(1:4,:);       
        Info.type='Inputs,1:9 and Output, 1:4, ( L1/L2 & W1/W2 + single_pressure + single_stability )';
        
    case 11 % width + single_pressure + single_stability
        index2=find(Out0(3,:)./Out0(4,:)<1);
        index=index2;        
        In=In0(1:9,index);
        Out=Out0(1:2,index);        
        Info.type='Inputs,1:9 and Output, 1:2, ( W1/W2 + single_pressure + single_stability )';
        Info.index2=index2;
        Info.index=index;
      
end

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





