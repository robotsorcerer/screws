%% The curvature form
function K = curvature_form(f)
% My implementation follows from Montana's definition
% Kinematics of Contact and Grasp, 1988
% Curvature form defined as
% Given coordinate axes x, y, z, find K. 
%  
% f: an invertible mapping from U \subset R^2 to S \subset R^3, 3x1
% symbolic
% x, y, z: vectors 3x1
% K = [x(u), y(u)]^T * [dz(u)/du)/norm($\partial f/\partial u$),
%                       dz(u)/du)/norm($\partial f/\partial v$)]$

syms u v R real
fu = diff(f, u); 
fv = diff(f, v);

x = simplify(fu/norm(fu))
y = simplify(fv/norm(fv)); 
z = simplify(f/R);

zu = diff(z, u);

K = [x, y].' * [zu/norm(fu), zu/norm(fv)]; 