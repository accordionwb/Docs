function h= csvreadh( filename, R1, delim )
%CSVREADH Read a comma separated value file with header.
%   [H,M] = CSVREADH('FILENAME') reads a comma separated value formatted file
%   FILENAME.  The result data is returned in M, the header in H. 
%   The file can only contain numeric values as data and a string for 
%   the header.

% Validate input args
if nargin==0
    narginchk(1,1); 
end

% Get Filename
if ~ischar(filename)
    error('csvreadh:FileNameMustBeString', ...
        'Filename must be a string.'); 
end

% Make sure file exists
if exist(filename,'file') ~= 2 
    error('csvreadh:FileNotFound',...
    'File not found.');
end

if nargin==1
    delim = ',';
    R1=0;
elseif nargin ==2
    delim=',';
end

% open input file
fid = fopen( filename );
read_line_num=0;
while read_line_num~=R1
line = fgetl( fid );
read_line_num=read_line_num+1;
end
h = regexp( line, delim, 'split' );


% m = [];
% % this is not quick for sure, but works
% while 1
%     line = fgetl( fid );
%     if ~ischar(line), break, end
%     m = [m; str2double(regexp( line, ',', 'split' ))];
% end

fclose(fid);