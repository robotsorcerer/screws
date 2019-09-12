function [ro, delta_volume, percent_volume] = isochoric_check(Ri, Ro, ri)
    %% Given a raference configuration radii and an internal radius in the current configuration
    % find the external radius we need to equilaterally expand the iab 
    % 
    ro = (Ro^3 + ri^3 - Ri^3)^(1/3);
    vol_change_ref = (4*pi/3) * (Ro^3-Ri^3);
    vol_change_cur = (4*pi/3) * (ro^3-ri^3);
    delta_volume = vol_change_cur - vol_change_ref;
    percent_volume = (delta_volume/vol_change_ref)*100;
end