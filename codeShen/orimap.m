function z = orimap()
    
    global z
    global grid_sz
    n=30;
    kc = 8/16;
    j = 0:1:n-1;
    
    rnd = rand(1,n)*75;
    l= rand(1,n);

    l(find(l>=0.5))=1;
    l(find(l<0.5))= -1;
    z = zeros(grid_sz,grid_sz);
    for i = 1:grid_sz
        for ii = 1:grid_sz
            
            z(i,ii) = sum(exp( sqrt(-1) .* (rnd+l.* kc .* ( cos(j*pi/n)*i + sin(j*pi/n)*ii) )));
        end
    end
     
    pref = rescale(mod(angle(z),2*pi),[0 2*pi],[1 256]);

    figure
    image(pref)
    ctab = fitzlabclut(256);
    colormap(ctab)
    colorbar
    
end