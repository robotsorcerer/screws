%% Exponential map of angle
function eangle = rodrigues(w, theta)
% Exponential map of angle
% w must be 3 x 1
assert(isequal([3,1], size(w)), 'w must be 3x1');

what = skewsem(w);

if ((norm(w)==1))
    eangle = eye(3) + what * sin(theta) + what^2 * (1 - cos(theta));
else
    eangle = eye(3) + what/norm(w) * sin(norm(w)*theta) + ...
             (what^2/norm(what)^2) * (1 - cos(norm(w)*theta));
end
end