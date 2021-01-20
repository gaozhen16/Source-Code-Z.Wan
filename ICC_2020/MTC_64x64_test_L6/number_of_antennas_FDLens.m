function N = number_of_antennas_FDLens(D, ele_max, azi_max)
    N = 0;
    for me = -floor(D*sin(ele_max/2)):floor(D*sin(ele_max/2))
       N = N + 1 + 2*floor(D*sin(azi_max/2)*cos(asin(me/D)));
    end
end