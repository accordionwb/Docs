function batch_plot(filelist,mode,output1,output2)
% mode includes: AVI  GIF   ALL
% output1 is the name of GIF file
% output2 is the name of AVI file
%% Jpeg and AVI Making
t0=cputime;
N=length(filelist);
% Jpeg making
if (strcmp(mode,'AVI') || strcmp(mode,'ALL'))  % OR logical
    for i=1:N
        inputname=['MatFig/',filelist(i).name];
        [pathstr,name,ext]=fileparts(inputname);
        outputname=['MatJpg/',name,'.jpg'];
        fid=openfig(inputname);
        set(fid,'position',[100,100,500,400],'color','w')  % 100,100 is the left-bottom corner position; 500,400 is the length and width of the figure
        drawnow
        frame=getframe(fid);
        im=frame2im(frame);
        [A,map]=rgb2ind(im,256);
        imwrite(A,map,outputname,'jpg','Comment','lillian-Jpeg'); % Writing jpeg file
        mov(i)=frame;
        close(fid)
    end
    aviname=output2;
    movie2avi(mov,aviname,'compression','None','fps',5); % writing AVI file, note: fps means frame per second, default is 15.
end
%% GIF making
if (strcmp(mode,'GIF') || strcmp(mode,'ALL')) % OR logical
   
    for i=1:N
        inputname=['MatFig/',filelist(i).name];
        fid=openfig(inputname);
        set(fid,'position',[100,100,1000,800],'color','w')
        drawnow
        frame=getframe(fid);
        im=frame2im(frame);
        [A,map]=rgb2ind(im,256);
        filename=output1;
        if i==1;  % the first frame
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.2); % Delaytime, seconds between two frame.
        else  % the following frame
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.2);
        end
        close(fid)
    end
end
%%
if (~strcmp(mode,'GIF') && ~strcmp(mode,'ALL') && ~strcmp(mode,'AVI')) % AND logical
    disp('Wrong mode, choose from "AVI", "GIF" and "ALL".')
end
t1=cputime-t0;
disp(['Total cpu time is: ',num2str(t1)])
return