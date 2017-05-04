function [cells,init_state]=fun_init_cells(XX,YY,N,type)
%% This function create all zero init_state for game of life
% and cooresponding cells with either random init state or random state
% Parameter defination:
% 1. XX: length on X axis
% 2. YY: length on Y axis
% 3. ZZ: length on Z axis (optional)
% 4. N : total number of ones
% 5. 'type':  String type of inilization
%       'centralized'
%       'normal'
%       'line'

init_state = zeros(XX,YY);
cells=init_state;
switch type
    case 'centralized'
        xn=XX*sqrt(N/XX/YY);
        yn=YY*sqrt(N/XX/YY);
        ratio=1.25*xn/XX;
        ix=randi([round((0.5-ratio/2)*XX),round((0.5+ratio/2)*XX)],round(xn),1);
        iy=randi([round((0.5-ratio/2)*YY),round((0.5+ratio/2)*YY)],round(yn),1);
        cells(ix,iy)=1;
        
    case 'normal'
        xn=round(XX*sqrt(N/XX/YY));
        yn=round(YY*sqrt(N/XX/YY));
        ix=round(random('norm',XX/2,XX/10,[xn,1]));
        iy=round(random('norm',YY/2,YY/10,[yn,1]));
        cells(ix,iy)=1;
        
    case 'line'
        
        nlines=randi(10); % random number of lines
        for i=1:nlines
            which_xy=randn;
            if which_xy>0  % choose x-axis
                ll=randi(XX);
                start_x=randi(XX-ll);
                start_y=randi(YY);
                cells(start_x:start_x+ll,start_y)=1;
                
            else  % y-axis
                ll=randi(YY);
                start_x=randi(XX);
                start_y=randi(YY-ll);
                cells(start_x,start_y:start_y+ll)=1;
            end
        end
        
        
end


