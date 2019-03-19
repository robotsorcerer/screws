%% Twist
function tw = twist(w, q)
% Computes the twist vector about axis w through a point q on w
% w and q must be 3 x 1
tw = [cross(-w, q); w];
end