function  h = gaussian_2d(sigma,x1,y1,x2,y2)
    
    GlobalVariables_orimap
    d = (x1-x2).^2 + (y1-y2).^2;
    
    h = exp(-1 * d / (2 * sigma^2)); % shape of external input

end