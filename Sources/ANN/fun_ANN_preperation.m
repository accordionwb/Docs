function [Input,Target,Test,Check,Info]=fun_ANN_preperation(Results1,Results2,Var)
%% Introduction
% This function is used to extract the generated data for neural network training.
% status:
%      0/1/2     'concentration'  /  'both' / 'width'
%      3/0     'pressure'  Y/N
%      6/0     'stability' Y/N

%% initialization
Input0=Results1.Input;
Target0=Results1.Output;
Test0=Results2.Input;
Check0=Results2.Output;
% Variable
mode=Var.mode;
single_pressure=Var.single_pressure;
single_stability=Var.single_stability;
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
        index1=find(Target0(1,:)./Target0(2,:)<1);
        index2=find(Target0(3,:)./Target0(4,:)<1);
        index_input=intersect(index1,index2);
        index3=find(Check0(1,:)./Check0(2,:)<1);
        index4=find(Check0(3,:)./Check0(4,:)<1);
        index_test=intersect(index3,index4);
        Input=Input0(1:11,index_input);
        Target=Target0(5:6,index_input);
        Test=Test0(1:11,index_test);
        Check=Check0(5:6,index_test);
        Info.type='Inputs,1:11 and Output, 5:6, ( concentration + all_pressure + all_stability )';
        Info.index1=index1;
        Info.index2=index2;
        Info.index_input=index_input;
        Info.index3=index3;
        Info.index4=index4;
        Info.index_test=index_test;
    case 1  % both L1/L2 W1/W2 + all_pressure + all_stability
        Input=Input0(1:11,:);
        Target=Target0(1:4,:);
        Test=Test0(1:11,:);
        Check=Check0(1:4,:);
        Info.type='Inputs,1:11 and Output, 1:4, ( L1/L2 & W1/W2 + all_pressure + all_stability )';
    case 2 % W1/W2 + all_pressure + all_stability
        index2=find(Target0(3,:)./Target0(4,:)<1);
        index_input=index2;
        index4=find(Check0(3,:)./Check0(4,:)<1);
        index_test=index4;
        Input=Input0(1:11,index_input);
        Target=Target0(1:2,index_input);
        Test=Test0(1:11,index_test);
        Check=Check0(1:2,index_test);
        Info.type='Inputs,1:11 and Output, 1:2, ( W1/W2 + all_pressure + all_stability )';
        Info.index2=index2;
        Info.index_input=index_input;
        Info.index4=index4;
        Info.index_test=index_test;
    case 3 % concentration + single_pressure + all stability
        index1=find(Target0(1,:)./Target0(2,:)<1);
        index2=find(Target0(3,:)./Target0(4,:)<1);
        index_input=intersect(index1,index2);
        index3=find(Check0(1,:)./Check0(2,:)<1);
        index4=find(Check0(3,:)./Check0(4,:)<1);
        index_test=intersect(index3,index4);
        Input=Input0([1:9 11],index_input);
        Target=Target0(5:6,index_input);
        Test=Test0([1:9 11],index_test);
        Check=Check0(5:6,index_test);
        Info.type='Inputs,[1:9 11] and Output, 5:6, ( concentration + single_pressure + all_stability )';
        Info.index1=index1;
        Info.index2=index2;
        Info.index_input=index_input;
        Info.index3=index3;
        Info.index4=index4;
        Info.index_test=index_test;
    case 4 %both + single_pressure + all_stability
        Input=Input0([1:9 11],:);
        Target=Target0(1:4,:);
        Test=Test0([1:9 11],:);
        Check=Check0(1:4,:);
        Info.type='Inputs,[1:9 11] and Output, 1:4, ( L1/L2 & W1/W2 + single_pressure + all_stability )';
    case 5 % width (W1/W2) + single_pressure + all_stability
        index2=find(Target0(3,:)./Target0(4,:)<1);
        index_input=index2;
        index4=find(Check0(3,:)./Check0(4,:)<1);
        index_test=index4;
        Input=Input0([1:9 11],index_input);
        Target=Target0(1:2,index_input);
        Test=Test0([1:9 11],index_test);
        Check=Check0(1:2,index_test);
        Info.type='Inputs,1:11 and Output, 1:2, ( W1/W2 + single_pressure + all_stability )';
        Info.index2=index2;
        Info.index_input=index_input;
        Info.index4=index4;
        Info.index_test=index_test;
    case 6 % concentration + all_pressure + single_stability
        index1=find(Target0(1,:)./Target0(2,:)<1);
        index2=find(Target0(3,:)./Target0(4,:)<1);
        index_input=intersect(index1,index2);
        index3=find(Check0(1,:)./Check0(2,:)<1);
        index4=find(Check0(3,:)./Check0(4,:)<1);
        index_test=intersect(index3,index4);
        Input=Input0(1:10,index_input);
        Target=Target0(5:6,index_input);
        Test=Test0(1:10,index_test);
        Check=Check0(5:6,index_test);
        Info.type='Inputs,1:10 and Output, 5:6, ( concentration + all_pressure + single_stability )';
        Info.index1=index1;
        Info.index2=index2;
        Info.index_input=index_input;
        Info.index3=index3;
        Info.index4=index4;
        Info.index_test=index_test;
    case 7  % both + all_pressure + single_stability
        Input=Input0(1:10,:);
        Target=Target0(1:4,:);
        Test=Test0(1:10,:);
        Check=Check0(1:4,:);
        Info.type='Inputs,1:10 and Output, 1:4, ( L1/L2 & W1/W2 + all_pressure + single_stability )';
    case 8  % width + all_pressure + single_stability
        index2=find(Target0(3,:)./Target0(4,:)<1);
        index_input=index2;
        index4=find(Check0(3,:)./Check0(4,:)<1);
        index_test=index4;
        Input=Input0(1:10,index_input);
        Target=Target0(1:2,index_input);
        Test=Test0(1:10,index_test);
        Check=Check0(1:2,index_test);
        Info.type='Inputs,1:10 and Output, 1:2, ( W1/W2 + all_pressure + single_stability )';
        Info.index2=index2;
        Info.index_input=index_input;
        Info.index4=index4;
        Info.index_test=index_test;
    case 9 %  concentration + single_pressure + single_stability
        index1=find(Target0(1,:)./Target0(2,:)<1);
        index2=find(Target0(3,:)./Target0(4,:)<1);
        index_input=intersect(index1,index2);
        index3=find(Check0(1,:)./Check0(2,:)<1);
        index4=find(Check0(3,:)./Check0(4,:)<1);
        index_test=intersect(index3,index4);
        Input=Input0(1:9,index_input);
        Target=Target0(5:6,index_input);
        Test=Test0(1:9,index_test);
        Check=Check0(5:6,index_test);
        Info.type='Inputs,1:9 and Output, 5:6, ( concentration + single_pressure + single_stability )';
        Info.index1=index1;
        Info.index2=index2;
        Info.index_input=index_input;
        Info.index3=index3;
        Info.index4=index4;
        Info.index_test=index_test;
    case 10 % both + single_pressure + single_stability
        Input=Input0(1:9,:);
        Target=Target0(1:4,:);
        Test=Test0(1:9,:);
        Check=Check0(1:4,:);
        Info.type='Inputs,1:9 and Output, 1:4, ( L1/L2 & W1/W2 + single_pressure + single_stability )';
    case 11 % width + single_pressure + single_stability
        index2=find(Target0(3,:)./Target0(4,:)<1);
        index_input=index2;
        index4=find(Check0(3,:)./Check0(4,:)<1);
        index_test=index4;
        Input=Input0(1:9,index_input);
        Target=Target0(1:2,index_input);
        Test=Test0(1:9,index_test);
        Check=Check0(1:2,index_test);
        Info.type='Inputs,1:9 and Output, 1:2, ( W1/W2 + single_pressure + single_stability )';
        Info.index2=index2;
        Info.index_input=index_input;
        Info.index4=index4;
        Info.index_test=index_test;
end

%% Random number of input
if (Var.random)
    Nmax1=length(Input);
    rtemp=randperm(Nmax1);
    Nin1=floor(Nmax1*Var.train_ratio(1));
    Nsn1=Var.train_ratio(2);
    Nin=max(Nin1,Nsn1);
    rdex=rtemp(1:Nin);
    for i=1:length(rdex)
        Input_t(:,i)=Input(:,rdex(i));
        Target_t(:,i)=Target(:,rdex(i));
    end
    
    Nmax2=length(Test);
    rtemp2=randperm(Nmax2);
    Nin2=floor(Nmax2*Var.test_ratio(1));
    Nsn2=Var.test_ratio(2);
    Ntn=max(Nin2,Nsn2);
    rdex2=rtemp2(1:Ntn);
    for j=1:length(rdex2)
        Test_t(:,j)=Test(:,rdex2(j));
        Check_t(:,j)=Check(:,rdex2(j));
    end
    Info.random='yes';
    Info.r_index1=rdex;
    Info.r_index2=rdex2;
    
    Input=Input_t;
    Target=Target_t;
    Test=Test_t;
    Check=Check_t;
end




