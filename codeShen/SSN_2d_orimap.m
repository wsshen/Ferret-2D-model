function [G_r,I,network_ratio] = SSN_2d_orimap(external_input,varargin)
    
    global grid_sz num_units
    global z
    global t dt
    global tau_E_m tau_I_m
    global k_m n_m
    global G_Wee G_Wei G_Wie G_Wii
    assign(varargin{:})
    G_r = gpuArray(zeros(2,grid_sz,grid_sz,length(t))); % firing rate of E and I units
    G_r_ss = gpuArray(zeros(2,grid_sz,grid_sz,length(t)));

    tau_E = 1./( 0.05*randn(grid_sz,grid_sz) * tau_E_m + tau_E_m);
    tau_I = 1./( 0.05*randn(grid_sz,grid_sz) * tau_I_m + tau_I_m);
    tau(1,:,:) = [tau_E];
    tau(2,:,:) = tau_I;
    
    k_temp = 0.05*randn(grid_sz,grid_sz) * k_m + k_m;
    k = repmat(reshape(k_temp,num_units,1),2,1);
    
    n_E = 0.05*randn(grid_sz,grid_sz) * n_m(1) + n_m(1);
    n_I = 0.05*randn(grid_sz,grid_sz) * n_m(2) + n_m(2);
    n = [reshape(n_E,num_units,1);reshape(n_I,num_units,1)];

    for i = 1:length(t)-1
                
        G_W = [G_Wee -1 * G_Wei ; G_Wie -1 * G_Wii];
        r_temp = [reshape(squeeze(G_r(1,:,:,i)),num_units,1);reshape(squeeze(G_r(2,:,:,i)),num_units,1)];

        input_fb = G_W * r_temp;

        I = external_input + input_fb;

        [neg_ind] = find(I<0);% thresholding the input function
        I(neg_ind) = 0;
        r_ss_temp = k .* (I).^n;
%             r_ss_temp = gather(Gr_ss_temp);
        G_r_ss(1,:,:,i) = reshape(r_ss_temp(1:num_units),grid_sz,grid_sz);
        G_r_ss(2,:,:,i) = reshape(r_ss_temp(num_units+1:2*num_units),grid_sz,grid_sz);
         % updating the firing rate
        dr = (-G_r(:,:,:,i) + G_r_ss(:,:,:,i)) * dt .* tau;

        G_r(:,:,:,i+1) = G_r(:,:,:,i) + dr;

    end
    input_network_E = G_Wee * reshape(squeeze(G_r(1,:,:,i)),num_units,1) + G_Wei * reshape(squeeze(G_r(2,:,:,i)),num_units,1);
    input_network_I = G_Wie * reshape(squeeze(G_r(1,:,:,i)),num_units,1) + G_Wii * reshape(squeeze(G_r(2,:,:,i)),num_units,1);

    total_input_E = input_network_E + external_input(1:num_units);
    total_input_I = input_network_I + external_input(1:num_units);
    network_ratio(1,:,:) = reshape(input_network_E./total_input_E,grid_sz,grid_sz);
    network_ratio(2,:,:) = reshape(input_network_I./total_input_I,grid_sz,grid_sz);
end
