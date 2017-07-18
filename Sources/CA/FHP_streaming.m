function nodes=FHP_streaming(ac_nodes)
% This function calculate the streaming/propagation of lattice speed
% Iterate over all the nodes, propagating the particles as we go.
[Ncell,NX,NY]=size(ac_nodes);

% Cell coordinate shift vectors:
% Basic configuration
% ck0 <=> mod(*,2) == 0
% ck1 <=> mod(*,2) == 1
c1=[0,1,0];
c20=[0,1,1]; 
c21=[0,0,1]; 
c30=[0,0,1];
c31=[0,-1,1];
c4=[0,-1,0];
c50=[0,0,-1];
c51=[0,-1,-1];
c60=[0,1,-1];
c61=[0,0,-1];

n2region=find(mod(1:NY,2) == 1);  % region mod(*,2) == 1
y2region=find(mod(1:NY,2) == 0);  % region mod(*,2) ~= 0 

% Extend nodes boundary by 1 on both X and Y directions
xx=2:NX+1;
yy=2:NY+1;
ex_nodes=false(Ncell,NX+2,NY+2);
ex_nodes(:,xx,yy)=ac_nodes;
% Perodic boundary copy:
ex_nodes(:,NX+2,yy)=ac_nodes(:,1,:); % left -> extend right
ex_nodes(:,1,yy)=ac_nodes(:,NX,:);   % right -> extend left
ex_nodes(1,2,yy)=true(1,NY);  % Inlet Boundary condition
ex_nodes(:,xx,NY+2)=ac_nodes(:,:,1); % bottom -> extend top
ex_nodes(:,xx,1)=ac_nodes(:,:,NY);   % top -> extend bottom


for k=Ncell:-1:1
    switch k
        case 1
            nodes(k,:,:) = circshift(ex_nodes(k,xx,yy),c1);
        case 2
            nodes(k,:,y2region) = circshift(ex_nodes(k,xx,y2region+1),c20);
            nodes(k,:,n2region) = circshift(ex_nodes(k,xx,n2region+1),c21);            
        case 3
            nodes(k,:,y2region) = circshift(ex_nodes(k,xx,y2region+1),c30);
            nodes(k,:,n2region) = circshift(ex_nodes(k,xx,n2region+1),c31);              
        case 4
            nodes(k,:,:) = circshift(ex_nodes(k,xx,yy),c4);
        case 5
            nodes(k,:,y2region) = circshift(ex_nodes(k,xx,y2region+1),c50);
            nodes(k,:,n2region) = circshift(ex_nodes(k,xx,n2region+1),c51);
        case 6
            nodes(k,:,y2region) = circshift(ex_nodes(k,xx,y2region+1),c60);
            nodes(k,:,n2region) = circshift(ex_nodes(k,xx,n2region+1),c61);
    end
end

% debug part
% stemp10=reshape(nodes(1,:,:),NX,NY);
% disp('finish streaming')

end