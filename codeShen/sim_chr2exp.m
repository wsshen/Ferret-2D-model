function [weights_1,weights_2,linear_sum,responses_both] = sim_chr2exp(varargin)

    % one call function example is --> sim_chr2exp('SubtractSpon',1,'sigma_ff',145,'cen_ex_ratio',1/3)
    
    GlobalVariables_orimap
    SubtractSpon = 0; % don't subtract spontaneous activities by default
    
    cen_ex_ratio = 0; % no central exclusion by default
    assign(varargin{:})
    % setting up the strings used for eval calls
    if SubtractSpon ==1
        add_noise_str = '+bg_noise;';
        minus_baseline = '-r_bg;';
    else
        add_noise_str = [];
        minus_baseline = [];
    end
    
    if cen_ex_ratio
        central_exclusion = ones(grid_sz,grid_sz);
        central_exclusion(round(grid_sz*cen_ex_ratio+1):round(grid_sz*cen_ex_ratio*2),round(grid_sz*cen_ex_ratio+1):round(grid_sz*cen_ex_ratio*2))=0; % assume that exclusion ratio is 1/3.
        cen_ex = repmat(reshape(central_exclusion,num_units,1),2,1);
        cen_ex_str = '.* cen_ex;';
    else
        cen_ex_str = [];
        
    end
    % simulate chr2 experiments with ortho column stimulation excluding the
    % central region

    
    z = orimap();
    Wee = SSN_2d_initconnections('EE');
    Wei = SSN_2d_initconnections('EI'); 
    Wie = SSN_2d_initconnections('IE');
    Wii = SSN_2d_initconnections('II');

    G_Wee = gpuArray(Wee);
    G_Wei = gpuArray(Wei);
    G_Wie = gpuArray(Wie);
    G_Wii = gpuArray(Wii);

    pref = mod(angle(z)/2,pi)*180/pi;

    ori = [0:3.4:86.6;90:3.4:90+86.6];
    ori_bin = [5:10:175];
    figure
    input_strength = [1:1:50];
    weights_1 = zeros(size(ori,2),length(input_strength));
    weights_2 = zeros(size(ori,2),length(input_strength));
    responses_both = zeros(size(ori,2),length(input_strength));
    linear_sum = zeros(size(ori,2),length(input_strength));

    bg_noise = rand([num_units*2,1])*7;
    [G_r_bg,G_I_bg,dummy] = SSN_2d_orimap(bg_noise);
    r_bg = gather(G_r_bg(:,:,:,end));

    for k = 1:size(ori,2)
        k
        for ii = 1:length(input_strength)
            ii 
            E_s1 = [];
            I_s1 = [];
            E_s2 = [];
            I_s2 = [];
            E_d2 = [];
            I_d2 = [];

            eval(['external_input1 = create_inputs_2d(ori(1,k),1,input_strength(ii))' add_noise_str]) ;
            [G_r1,G_I,G_networkinput_ratios] = SSN_2d_orimap(external_input1) ;
            eval(['r_single1 = gather(G_r1(:,:,:,end))' minus_baseline]);

            E_single1 = squeeze(r_single1(1,:,:,end));
            I_single1 = squeeze(r_single1(2,:,:,end));
            for jj = ori_bin
                [indx,indy] = find(abs(pref - jj)<5);
                E_s1(end+1) = mean(diag(E_single1(indx,indy)));
                I_s1(end+1) = mean(diag(I_single1(indx,indy)));
            end
            eval(['external_input2 = create_inputs_2d(ori(2,k),1,input_strength(ii))' cen_ex_str add_noise_str]) ;
            [G_r2,G_I,G_networkinput_ratios] = SSN_2d_orimap(external_input2);
            eval(['r_single2 = gather(G_r2(:,:,:,end))' minus_baseline]) ;

            E_single2 = squeeze(r_single2(1,:,:,end));
            I_single2 = squeeze(r_single2(2,:,:,end));
            for jj = ori_bin
                [indx,indy] = find(abs(pref - jj)<5);
                E_s2(end+1) = mean(diag(E_single2(indx,indy)));
                I_s2(end+1) = mean(diag(I_single2(indx,indy)));
            end
            external_input3 = external_input1 + external_input2;
            [G_r3,G_I,G_networkinput_ratios] = SSN_2d_orimap(external_input3);
            eval(['r_double = gather(G_r3(:,:,:,end))' minus_baseline]) ;

            E_double2 = squeeze(r_double(1,:,:,end));
            I_double2 = squeeze(r_double(2,:,:,end));

            for jj = ori_bin
                [indx,indy] = find(abs(pref - jj)<5);
                E_d2(end+1) = mean(diag(E_double2(indx,indy)));
                I_d2(end+1) = mean(diag(I_double2(indx,indy)));
            end

            [fitobject_E,gof,output] = fit_weighted_sum(E_s1',E_s2',E_d2');
            [fitobject_I,gof,output] = fit_weighted_sum(I_s1',I_s2',I_d2');

            weights_1(k,ii) = fitobject_E.a;
            weights_2(k,ii) = fitobject_I.a;

            [indx,indy] = find(abs(pref - ori(1,k))<5);
            E_pref_temp = mean(diag(E_single1(indx,indy)));

            [indx,indy] = find(abs(pref - ori(1,k))<5);
            E_ortho_temp = mean(diag(E_single2(indx,indy)));

            [indx,indy] = find(abs(pref - ori(1,k))<5);
            E_both_temp = mean(diag(E_double2(indx,indy)));

            responses_both(k,ii) = E_both_temp;
            linear_sum(k,ii) =  E_pref_temp + E_ortho_temp;
        %     
        %     fb(1,:,:) = reshape(fb_(1:num_units),grid_sz,grid_sz);
        %     fb(2,:,:) = reshape(fb_(num_units+1:2*num_units),grid_sz,grid_sz);

        %     plot(i,mean(diag(squeeze(r(1,indx,indy)))),'ro')
        %     plot(i,100*mean(diag(squeeze(ratios(1,indx,indy)))),'ro','MarkerFaceColor','r')
        %     hold on
        % %     plot(i,mean(diag(squeeze(r(2,indx,indx)))),'b*')
        %     plot(i,100*mean(diag(squeeze(ratios(2,indx,indy)))),'bo','MarkerFaceColor','b')
        %     plot(i,100-100*mean(diag(squeeze(ratios(1,indx,indy)))),'r*')
        %     plot(i,100-100*mean(diag(squeeze(ratios(2,indx,indy)))),'b*')
        % 
        % 
        %     legend('network_E','network_I','external_E','external_I')
        %     xlabel('Strength')
        %     ylabel('Percent of input')
        %     ylim([0 100])
        end
    end

    plot(input_strength,mean(weights_1),'g')
    hold on 
    plot(input_strength,mean(weights_2),'k')
    legend('E','I')
    xlabel('strength')
    ylabel('relative weight')

    save('SSN_2d.mat','weights_1','weights_2','responses_both','linear_sum')

end