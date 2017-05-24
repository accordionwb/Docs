%% Game of life 
% 1. User specifications
XX=800;
YY=400;
ZZ=10;

Or_value=3;
And_value=2;

% 2. Initialize Cells
% for 2-d cell
N=6000;
type='line';  
%       'centralized'
%       'normal'
%       'line'
[cells,init_state]=fun_init_cells(XX,YY,N,type);

% 3. INitialize Image
% % im = image(cat(3,cells',cells',cells'));
im = image(cells');
% im = image(cells,'AlphaDataMapping','none');
axis equal
axis tight


% 4. build the GUI
%define the plot button
plotbutton=uicontrol('style','pushbutton',...
   'string','Run', ...
   'fontsize',12, ...
   'position',[20,100,50,20], ...
   'callback', 'run=1;');

%define the pause button
erasebutton=uicontrol('style','pushbutton',...
   'string','Pause', ...
   'fontsize',12, ...
   'position',[20,140,50,20], ...
   'callback','freeze=1;');

%define the Quit button
quitbutton=uicontrol('style','pushbutton',...
   'string','Quit', ...
   'fontsize',12, ...
   'position',[20,180,50,20], ...
   'callback','quit=1;close;');

number = uicontrol('style','text', ...
    'string','1', ...
   'fontsize',12, ...
   'position',[20,220,50,20]);

reloadbutton=uicontrol('style','pushbutton',...
    'string','Reload',...
    'fontsize',12,...
    'position',[20,260,60,20],...
    'callback','reload=1;');


% 5. Launch the GUI 
quit= 0; %wait for a quit button push
run = 0; %wait for a draw 
freeze = 0; %wait for a freeze
reload = 0;
% index changes
ix = 2:XX-1;
iy = 2:YY-1;
while (quit==0) 
       if (run==1)
	      %nearest neighbor sum
          sum=init_state;  % initialize sum state
          sum(ix,iy) = cells(ix,iy-1) + cells(ix,iy+1) + ...
              cells(ix-1, iy) + cells(ix+1,iy) + ...
              cells(ix-1,iy-1) + cells(ix-1,iy+1) + ...
              cells(ix+1,iy-1) + cells(ix+1,iy+1);
          
          sum(1,iy) = cells(1,iy-1)+cells(1,iy+1) + ...
              cells(XX,iy) + cells(2,iy) + ...
              cells(XX,iy-1) + cells(XX,iy+1) + ...
              cells(2,iy-1)+cells(2,iy+1);
          sum(XX,iy) = cells(XX,iy-1)+cells(XX,iy+1) + ...
              cells(XX-1,iy) + cells(1,iy) + ...
              cells(XX-1,iy-1) + cells(XX-1,iy+1) + ...
              cells(1,iy-1)+cells(1,iy+1);
          sum(ix,1)= cells(ix,YY) + cells(ix,2) + ...
              cells(ix-1,1) + cells(ix+1,1) + ...
              cells(ix-1,YY) + cells(ix-1,2) + ...
              cells(ix+1,YY) + cells(ix+1,2);
          sum(ix,YY)= cells(ix,YY-1) + cells(ix,1) + ...
              cells(ix-1,YY) + cells(ix+1,YY) + ...
              cells(ix-1,YY-1) + cells(ix-1,1) + ...
              cells(ix+1,YY-1) + cells(ix+1,1);
	      % The CA rule 
	      cells = (sum==Or_value) | (sum==And_value & cells);       
          
	      %draw the new image
	      set(im, 'cdata', cat(3,cells',cells',cells') )
          
	      %update the step number diaplay
	      stepnumber = 1 + str2num(get(number,'string'));
	      set(number,'string',num2str(stepnumber))
      end
      if (freeze==1)
	     run = 0;
	     freeze = 0;
      end
      drawnow  %need this in the loop for controls to work  
      if reload==1
          [cells,init_state]=fun_init_cells(XX,YY,N,type);
          stepnumber=1;
          set(number,'string',num2str(stepnumber));
          reload=0;
      end     
      
end

%% for function test use
