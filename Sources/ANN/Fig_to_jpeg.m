%% change dir  lillian modified here
cd f:/lillian-fig
%% rename and sort
sourcename=dir(strcat('MatFig', '/*.fig'));
m=1;n=1;
file_fL(27*(m-1)+1:27*m)=flipud(sourcename(27*(n-1)+1:27*n));
m=2;n=2;
file_fL(27*(m-1)+1:27*m)=sourcename(27*(n-1)+1:27*n);
m=1;n=3;
file_fR(27*(m-1)+1:27*m)=flipud(sourcename(27*(n-1)+1:27*n));
m=2;n=4;
file_fR(27*(m-1)+1:27*m)=sourcename(27*(n-1)+1:27*n);
m=1;n=5;
file_t(27*(m-1)+1:27*m)=flipud(sourcename(27*(n-1)+1:27*n));
m=2;n=6;
file_t(27*(m-1)+1:27*m)=sourcename(27*(n-1)+1:27*n);
 %%  Function
  batch_plot(file_fL,'AVI','fl_figure.gif','fl_movie.avi')
%   batch_plot(file_fR,'ALL','fR_figure.gif','fR_movie.avi')
%   batch_plot(file_t,'ALL','t_figure.gif','t_movie.avi')
