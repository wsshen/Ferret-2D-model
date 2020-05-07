function [h] = create_inputs_2d(ori,j,c,varargin)
    % j represents either E or I receives the input. 1 is E, and 2 is I. 
    
    GlobalVariables_orimap
%     assign(varargin{:})
    
    h_temp = c * gaussian_(sigma_ff,ori,mod(angle(z)/2,pi)*180/pi);

    h = repmat(reshape(h_temp,num_units,1),2,1);
end