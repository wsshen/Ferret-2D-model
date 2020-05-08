function W = SSN_2d_initconnections(cmark)
    
    global num_units grid_sz
    global J sigma_J
    global z
    global sigma_ori
    global kk
    
    
    a = cmark(1);
    if a =='E'
        a = 1;
    else a = 2;
    end
    
    b = cmark(2);
    if b =='E'
        b=1;
    else b=2;
    end
    W = zeros(num_units,num_units);

    W_temp = 0.25 * randn(num_units,num_units) *J(a,b) + J(a,b);
    W_temp(find(W_temp<0))=0;    

    oripref_vec = reshape(mod(angle(z)/2,pi)*180/pi,num_units,1);
    x_vec = repmat([1:grid_sz],1,grid_sz); % create a row vector for x
    y_vec = reshape(repmat([1:grid_sz],grid_sz,1),1,num_units); % create a row vector for y
    
    p = kk(b) * gaussian_(sigma_ori,repmat(oripref_vec,1,num_units),repmat(oripref_vec',num_units,1))...
        .* gaussian_2d(sigma_J(a,1),repmat(x_vec',1,num_units),repmat(y_vec',1,num_units),repmat(x_vec,num_units,1),repmat(y_vec,num_units,1));
    
    W = p.* W_temp;

end