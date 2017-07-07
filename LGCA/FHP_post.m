function av_vel=FHP_post(nodes, grain_size)

[Ncell,NX,NY]=size(nodes);

% Define the lattice velocities.
for k=Ncell:-1:1
    lvc(k,1)=cos(pi*(k-1)/3);
    lvc(k,2)=sin(pi*(k-1)/3);
end
% Subdivide the total domain into subdomains of size 32x32 for the
% purposes of coarse-graining.  See pg.  51.
% grain_size = 8;
grain_x = NX / grain_size;
grain_y = NY / grain_size;


% Iterate over the entire domain, averaging and storing the results as
% we go.

for i = grain_x:-1:1
    % Calculate the lower and upper x-boundaries.
    x_cell = (i - 1)*grain_size + 1: i*grain_size;
    for j = grain_y:-1:1
        % Calculate the lower and upper y-boundaries.
        y_cell = (j - 1)*grain_size + 1: j*grain_size;
        
        % Get the number of particles moving in each direction in the
        % current subdomain.
        for k=1:Ncell
        np(k) = sum(sum(nodes(k,x_cell, y_cell)));
        end
        % Compute the average velocity.
        av_vel(i,j,:) = (1/(grain_size.^2))*np*lvc;
        
    end
end



end
