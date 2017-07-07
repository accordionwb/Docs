% 3D array of nodes to store the vectors that represent the occupied
% cells at each node.
%   0 - Cell unoccupied.
%   1 - Cell occupied.
%
% 1st Index -- Node x-coordinate.
% 2nd Index -- Node y-coordinate.
% 3rd Index -- Cell number
%
% The elements of the occupancy vectors correspond to the cells in the
% following way:
%
%              3   2
%               \ /
%            4 - O - 1
%               / \
%              5   6
%
% Observe that this convention differs slightly from Wolf-Gladrow's.

%% % control Parameter
clear
clc
% FHP method built-in character
Ncell = 6;

% Number of nodes in each direction.
% These must be multiples of 32 for the coarse graining to work.
% Unit: 2^n for coarse graining
% Ratio: M vs. N multipliers
unit=5;
Ratio=[20,8];

% Define fluid region, Logical values, 1 = fluid, 0 = obst/boundary
% fluid = true(NX, NY);
% Insert a flat plate as the obstacle.
R_obstx=[1/6,1/6];
R_obsty=[1/3,2/3];  % line obst
R_obstr=1/6;   % circle obst
type = 1;  % 0 == plate ; 1 == circle

% Number of timesteps over which to simulate.
t_end = 1000;

% plot results ?
is_plot=false;

% For results refinement
grain_size=8;
t_plot=10;

%% Execuate the program
% X and Y dimensions
NX = Ratio(1)*2^unit;
NY = Ratio(2)*2^unit;

% Insert a circular cylinder as the obstacle.
[Y,X]=meshgrid(1:NY,1:NX);
if type == 0
    iOBST = round([NX*R_obstx(1), NX*R_obstx(2), NY*R_obsty(1), NY*R_obsty(2)]);
    fluid = X < iOBST(1) | X > iOBST(2) | Y < iOBST(3) | Y > iOBST(4);
else
    theta = 0:0.001:2*pi;
    xc = round(R_obstx(1)*NX);
    yc = round(1/2*NY);
    rc = round(R_obstr*NY);
    fluid = (X-xc).^2+(Y-yc).^2 > rc.^2;
end

% Top and Bottom solid boundary
fluid(:,[1,NY]) = false;

% Initialize conditions
nodes = false(Ncell,NX,NY);
ffregion=find(fluid);
bbregion=find(~fluid);
% nodes(1,ffregion) = true; % initialize all horizontal speed
nodes(1,1,:)=true;


% Start the main loop over time steps.
result=cell(t_end,2);
tic; % Time program exectution.
for t = 1:t_end
    % Collision
    ac_nodes=FHP_collision(nodes,fluid);
    
    % Streaming
    nodes = FHP_streaming(ac_nodes);
    
    % Averaging
    av_vel=FHP_post(nodes, grain_size);
    
    % Save result
    result{t,1}=t;
    result{t,2}=av_vel;
    
    if mod(t,t_plot) == 0
        disp(['Iteration: ',num2str(t),'/',num2str(t_end)])
        
        % Plot figure
        if is_plot
            mx=size(av_vel,1);
            my=size(av_vel,2);
            
            % Pre-allocate vectors for the averaged velocities.
            [MY, MX] = meshgrid(1:my,  1:mx);
            % Store the velocity components.
            av_vel_x_comps = av_vel(:,:,1);
            av_vel_y_comps = av_vel(:,:,2);
            
            % Plot the average velocity field.
            quiver(MX, MY, av_vel_x_comps, av_vel_y_comps);
            
            % Plot the channel boundaries.
            hold on;
            plot([1, mx], [0.75, 0.75], 'k-');
            plot([1, mx], [my + 0.25, my + .25], 'k-');
            
            % Display the flow obstacle.
            
            k = 1;
            
            for i = 1:NX
                for j = 2:NY-2
                    if  fluid(i, j) == 0
                        obs_x(k) = 0.5 + (NX ./ (grain_size .* (NX - 1))) .* (i - 1);
                        obs_y(k) = 0.5 + (NY ./ (grain_size .* (NY - 1))) .* (j - 1);
                        k = k + 1;
                    end
                end
            end
            
            
            plot(obs_x, obs_y, 'r-');
            hold off
            axis equal;
            drawnow
            
        end
    end
    
end
save('results.mat','result');
toc; % Print the time it took to execute.

