% SSN_2d_orimap test
clear all
close all

GlobalVariables_orimap


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

ori = [0:8.5:85;90:8.5:85+90];
ori_bin = [5:10:175];
figure
input_strength = [1:1:50];
weights_1 = zeros(size(ori,2),length(input_strength));
weights_2 = zeros(size(ori,2),length(input_strength));

for k = 1:size(ori,2)
    for i = 1:length(input_strength)
        i 
        E_s1 = [];
        I_s1 = [];
        E_s2 = [];
        I_s2 = [];
        E_d2 = [];
        I_d2 = [];

        external_input1 = create_inputs_2d(ori(1,k),1,input_strength(i));
        [G_r1,G_I,G_networkinput_ratios] = SSN_2d_orimap(external_input1);
        r_single1 = gather(G_r1(:,:,:,end));

        E_single1 = squeeze(r_single1(1,:,:,end));
        I_single1 = squeeze(r_single1(2,:,:,end));
        for j = ori_bin
            [indx,indy] = find(abs(pref - j)<5);
            E_s1(end+1) = mean(diag(E_single1(indx,indy)));
            I_s1(end+1) = mean(diag(I_single1(indx,indy)));
        end
        external_input2 = create_inputs_2d(ori(2,k),1,input_strength(i));
        [G_r2,G_I,G_networkinput_ratios] = SSN_2d_orimap(external_input2);
        r_single2 = gather(G_r2(:,:,:,end));

        E_single2 = squeeze(r_single2(1,:,:,end));
        I_single2 = squeeze(r_single2(2,:,:,end));
        for j = ori_bin
            [indx,indy] = find(abs(pref - j)<5);
            E_s2(end+1) = mean(diag(E_single2(indx,indy)));
            I_s2(end+1) = mean(diag(I_single2(indx,indy)));
        end
        external_input3 = external_input1 + external_input2;
        [G_r3,G_I,G_networkinput_ratios] = SSN_2d_orimap(external_input3);
        r_double = gather(G_r3(:,:,:,end));

        E_double2 = squeeze(r_double(1,:,:,end));
        I_double2 = squeeze(r_double(2,:,:,end));

        for j = ori_bin
            [indx,indy] = find(abs(pref - j)<5);
            E_d2(end+1) = mean(diag(E_double2(indx,indy)));
            I_d2(end+1) = mean(diag(I_double2(indx,indy)));
        end

        [fitobject_E,gof,output] = fit_weighted_sum(E_s1',E_s2',E_d2');
        [fitobject_I,gof,output] = fit_weighted_sum(I_s1',I_s2',I_d2');

        weights_1(k,i) = fitobject_E.a;
        weights_2(k,i) = fitobject_I.a;

        I_ = gather(G_I);
        ratios = gather(G_networkinput_ratios);

        I(1,:,:) = reshape(I_(1:num_units),grid_sz,grid_sz);
        I(2,:,:) = reshape(I_(num_units+1:2*num_units),grid_sz,grid_sz);



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
   