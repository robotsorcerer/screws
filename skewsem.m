function what = skewsem(w)
    % generates a skew symmetric matrix from a vector w
    what = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end