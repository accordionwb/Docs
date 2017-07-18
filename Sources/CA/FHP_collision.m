function ac_nodes=FHP_collision(nodes,fluid)
% The function returns 3D array contains 2D FHP logical values
% Data Structure:
% nodes(i,j,k)
%       i = 1:NX
%       j = 1:NY
%       k = 1:6    the six cells in a node
%
% Cell configuration:
%              3   2
%               \ /
%            4 - O - 1
%               / \
%              5   6
%
%
% Routine 1. Boolean Algebra to calculate all possible situations
% Collision type described as follows:
% 
% a) 2-particle head-on collision
%                          \ //
%                         - O -
%                         // \
%               \ /   
%              = O =
%               / \
%                         \\ /
%                         - O -
%                          / \\
%
% b) symmetric 3-particle collisions
%               \ //            \\ /
%              = O -    <=>     - O =  
%               / \\            // \
%
% c) 4-particle head-on collision
%                         \\ //
%                         - O -
%                         // \\
%               \ //
%              = O =
%              // \
%                         \\ /
%                         = O =
%                          / \\
%
% d) 2-particle head-on collisions with seperator
%               \ /             \\ /
%              = O =    <=>     - O -  
%              // \             // \\
%%



% Parameter initialization
[Ncell,NX,NY]=size(nodes);

% numparts = [zeros(NX,1),sum(nodes(:),3),zeros(NX,1)];

% Build double nodes on k to iterate over 2*Ncell
db_nodes=nodes;
db_nodes(Ncell+1:2*Ncell,:,:)=nodes;
% Define probability array for collision type (a)
RA = logical(round(rand(Ncell,NX,NY)));
% Probability array for collision type (c)
% RC = logical(round(rand(NX,NY,Ncell)));

% opposite rotation index
rott=[4 5 6 1 2 3];
bbregion = find(~fluid);
% Loop over cells for in-fluid collisions
for k=Ncell:-1:1
    
    % Define type (b) collisions
    LB1=xor(db_nodes(k,:,:),db_nodes(k+1,:,:));
    LB2=xor(db_nodes(k+1,:,:),db_nodes(k+2,:,:));
    LB3=xor(db_nodes(k+2,:,:),db_nodes(k+3,:,:));
    LB4=xor(db_nodes(k+3,:,:),db_nodes(k+4,:,:));
    LB5=xor(db_nodes(k+4,:,:),db_nodes(k+5,:,:));
    LB =  LB1 & LB2 & LB3 & LB4 & LB5;
    
    % Define type (a) collision
    LA1=db_nodes(k,:,:) & db_nodes(k+3,:,:) & ~(db_nodes(k+1,:,:) |...
        db_nodes(k+2,:,:) | db_nodes(k+4,:,:) | db_nodes(k+5,:,:));
    
    LA2=RA(k,:,:) & db_nodes(k+1,:,:) & db_nodes(k+4,:,:) &...
        ~(db_nodes(k,:,:) | db_nodes(k+2,:,:) | db_nodes(k+3,:,:) | db_nodes(k+5,:,:));
    
    LA3=~RA(k,:,:) & db_nodes(k+2,:,:) & db_nodes(k+5,:,:) &...
        ~(db_nodes(k,:,:) | db_nodes(k+1,:,:) | db_nodes(k+3,:,:) | db_nodes(k+4,:,:));
    
    Change= (LB | LA1 | LA2 | LA3) & reshape(fluid,1,NX,NY);
    
    ac_nodes(k,:,:) = xor(nodes(k,:,:), Change);
   
   
    % Bounce-back boundary conditions at fluid boundaries and obstacles
    ac_nodes(k,bbregion) = nodes(rott(k),bbregion);
end

% debug part
% actmp1=reshape(ac_nodes(1,:,:),NX,NY);
% actmp2=reshape(ac_nodes(2,:,:),NX,NY);
% actmp3=reshape(ac_nodes(3,:,:),NX,NY);
% actmp4=reshape(ac_nodes(4,:,:),NX,NY);
% actmp5=reshape(ac_nodes(5,:,:),NX,NY);
% actmp6=reshape(ac_nodes(6,:,:),NX,NY);
% 
% temp1=reshape(nodes(1,:,:),NX,NY);
% 
% disp('finish collision')

end