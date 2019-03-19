%
function etwist = twistmap(w, v theta)
% Exponential map of a twist
% Corresponds to the relative motion of a point from an initial coordinate
% w.r.t to a final coordinate
% w and v must be 3 x 1
assert(isequal(size(w)), [3,1], 'w must be 3x1');
assert(isequal(size(v)), [3,1], 'v must be 3x1');

eangle = rodrigues(w, theta);
etwist = [eangle, (eye(3) - eangle)*cross(w, v);...
           zeros(1, 3), 1];
end
