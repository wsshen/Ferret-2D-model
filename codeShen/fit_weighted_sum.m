function [fitobject,gof,output] = fit_weighted_sum(X,Y,Z)
    
    weighted_sum = 'a*x+b*y';

    [fitobject,gof,output] = fit([X,Y],Z,weighted_sum,'MaxIter',1000);

end