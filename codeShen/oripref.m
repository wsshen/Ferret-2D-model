function [pref_x] = oripref(x,y)

    global z
    pref = angle(z);
    pref_x = pref(x,y);


end