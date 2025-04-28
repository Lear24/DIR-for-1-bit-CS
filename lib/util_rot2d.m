function v2 = util_rot2d(v1,theta)
 M = [cos(theta),-sin(theta);sin(theta),cos(theta)];
v2 = M*v1;

end