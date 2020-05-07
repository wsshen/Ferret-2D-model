function  h = gaussian_(sigma,x,y)

    h = exp(-1 * min(abs(x-y),180-abs(x-y)).^2 / (2 * sigma^2)); % shape of external input

end