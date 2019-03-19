function what = skewsem(w)
    % generates a skew symmetric matrix from a vector w
    what = [0 -w3 w2; w3 0 -w1; -w2 w1 0];
end