%% Function of random point generation
function [Xr,Yr]=fun_down_point_gen(X_reg,Y_reg,N)

%% TEST01 tests I4_TO_HAMMERSLEY_SEQUENCE.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 May 2008
%
%  Author:
%
%    John Burkardt
%

nmax = N;

dim_num = 2;
n = N;
step = 1;
seed(1:dim_num) = 0;
leap(1:dim_num) = 1;
base(1) = -nmax;
for i = 2 : dim_num
    base(i) = prime ( i - 1 );
end

r = i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base );
Xr=X_reg(1)+(X_reg(2)-X_reg(1))*r(1,:);
Yr=Y_reg(1)+(Y_reg(2)-Y_reg(1))*r(2,:);
end
