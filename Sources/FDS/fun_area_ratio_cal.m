function ratio=fun_area_ratio_cal(X_length,Y_length,Cubic_Data,Cylinder_Data,Sphere_Data)
% This function calculate the empty area propotion amoung the whole
% calcuation domain. 
% Input: X_length, Y_length | Horizantal length and width
% Cubic, Cylinder, Sphere configurations
total_area=X_length*Y_length;

r_cyl=Cylinder_Data(:,7);
r_sph=Sphere_Data(:,6);
x_cub=Cubic_Data(:,2)-Cubic_Data(:,1);
y_cub=Cubic_Data(:,4)-Cubic_Data(:,3);

obst_area=sum(pi*r_cyl.^2)+sum(pi*r_sph.^2)+sum(x_cub.*y_cub);
ratio=obst_area/total_area;