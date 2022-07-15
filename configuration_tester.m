%% Initialization
c = 343;

%% Array-matrix
row_elements = 8;
column_elements = 8;
r_prime = zeros(3,row_elements*column_elements);

%% Evenly spaced matrix initialization
%matrix array in the xy-plane
%uni_distance = 2*pi/2;      % Give length in kd
%uni_distance = 2*pi*lambda_d;
distance_mm = 10^-3 * 20;

max_frequency = c/distance_mm/2;

freq_res = 200;
freq_band = linspace(100,max_frequency,freq_res);
mode_max = 7;
data_output = zeros(freq_res,4);

for freq_ind = 1:freq_res
    freq_ind
    best_directivity = 0;
    data_output(freq_ind,1) = freq_band(freq_ind);

    for mode = 1:mode_max
        mode

        frequency = freq_band(freq_ind);
        uni_distance = 2*pi*frequency/c*distance_mm;
        lambda = c/frequency;
        lambda_rel = 0.02/lambda;
        
        element_index = 1;
        for i = 1:row_elements
            for j = 1:column_elements
                r_prime(1,element_index) = i*uni_distance;
                r_prime(2,element_index) = j*uni_distance;
                element_index = element_index + 1;
            end
        end
        
        %set origin of the matrix to (0,0)
        r_prime(1,:) = r_prime(1,:)-row_elements*uni_distance/2 - uni_distance/2;
        r_prime(2,:) = r_prime(2,:)-column_elements*uni_distance/2 - uni_distance/2;
        
        %figure(1);
        %plot(r_prime(1,:),r_prime(2,:),'linestyle','none','marker','*');
        
        %% Phase and amplitude correction
        theta_listen_deg = 0;
        phi_listen_deg = 0;
        
        theta_l = theta_listen_deg/180*pi;
        phi_l = phi_listen_deg/180*pi;
        phase_corr = zeros(1,row_elements*column_elements);
        
        element_index = 1;
        for i = 1:row_elements*column_elements
            
            phase_corr(element_index) =  exp(-1i*(r_prime(1,i)*sin(theta_l)*cos(phi_l) + ...
                    r_prime(2,i)*sin(theta_l)*sin(phi_l)));
                element_index = element_index + 1;
        end
        
        phase_amplitudes = zeros(1,row_elements*column_elements);
        
        
        %Binomial distribution      Reduces the sidelobe levels drastically
        element_index = 1;
        for i = 0:row_elements-1
            row_amp = nchoosek(row_elements-1,i);
            for j = 0:column_elements-1
                column_amp = nchoosek(column_elements-1,j);
                phase_amplitudes(element_index) = row_amp*column_amp;
                element_index = element_index +1;
            end
        end
        %phase_corr = phase_amplitudes.*phase_corr;
        
        %On-Off amplitude distribution
        phase_amplitudes_OF = zeros(1,row_elements*column_elements);
        row_lim = ceil((row_elements)/mode);
        column_lim = ceil((column_elements)/mode);
        for i = 1:row_lim
            for j = 1:column_lim
                element_index = (mode*(i-1))*row_elements + mode*(j-1) +1;
                phase_amplitudes_OF(element_index) = 1;
            end
        end
        
        phase_corr = phase_amplitudes_OF.*phase_corr;
        %figure(39)
        %plot(r_prime(1,:).*phase_amplitudes_OF,r_prime(2,:).*phase_amplitudes_OF,'linestyle','none','marker','*');
        
        %% Get far-field
        resolution = 200;
        [AF,theta,phi] = array_factor(r_prime,phase_corr,resolution);
        
        %% Radiated power
        antenna_element = 1;                
        power = abs(AF.*antenna_element).^2;    %P = |AF|^2
        integrand = ((sin(theta))').*power;     %Integral to be calculated: S |AF|^2 sin(theta)dthetadphi
        
        temp = cumtrapz(phi,cumtrapz(theta,integrand',2));
        power_rad = (temp(end,end) - temp(1,1))/(2*pi);
        
        AF_dbi = 10*log10(power/power_rad);     %Power in terms of isotropic radiation

        max_dB = max(max(AF_dbi));
        test = 1;
        if max_dB > best_directivity && mode*lambda_rel <= 0.5 
            data_output(freq_ind,2) = mode;
            data_output(freq_ind,3) = max_dB;
            data_output(freq_ind,4) = lambda_rel;
            best_directivity = max_dB;
        end
    end
end

%% Plotting 
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;0.3010 0.7450 0.9330];
figure(1)
%plot(freq_config_mode_20_mm,dbi_config_mode_20_mm,'color',colors(1,:),'LineWidth',0.6);
%hold on
plot(freq_20_mm,dbi_20_mm,'color',colors(2,:),'LineWidth',0.9)
hold on
plot(freq_30_mm,dbi_30_mm,'color',colors(3,:),'LineWidth',0.9)
hold on
plot(freq_40_mm,dbi_40_mm,'color',colors(4,:),'LineWidth',0.9)
xlim([freq_config_mode_20_mm(1) freq_config_mode_20_mm(end)])
ylim([min(dbi_config_mode_20_mm) max(dbi_config_mode_20_mm)])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('Frequency (Hz)','Interpreter','latex','FontSize',18);
ylabel('Broadside gain (dBi)','Interpreter','latex','FontSize',18)
legend('d = 20 mm','d = 30 mm','d = 40 mm','Interpreter','latex','location','southeast')
%xticks([0 30 60 90 120 150 180 210 240 270 300 330 360])
grid on
set(gca,'GridALpha',0.1,'LineWidth',.3);
exportgraphics(gcf,'unifrom_dist_plot.png','Resolution',300)


figure(2)
plot(freq_config_mode_20_mm,dbi_config_mode_20_mm,'color',colors(1,:),'LineWidth',0.9);
hold on
plot(freq_20_mm,dbi_20_mm,'color',colors(2,:),'LineWidth',0.9)

xlim([freq_config_mode_20_mm(1) 3500])
ylim([min(dbi_config_mode_20_mm) 12])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('Frequency (Hz)','Interpreter','latex','FontSize',18);
ylabel('Broadside gain (dBi)','Interpreter','latex','FontSize',18)
legend('adaptive config on','adaptive config off','Interpreter','latex','location','southeast')
%xticks([0 30 60 90 120 150 180 210 240 270 300 330 360])
grid on
set(gca,'GridALpha',0.1,'LineWidth',.3);
exportgraphics(gcf,'config_performance.png','Resolution',300)

%% Plot geometry
colors_alpha = [0 0.4470 0.7410 0.1; 0.8500 0.3250 0.0980 0.5;0.9290 0.6940 0.1250 0.5 ...
    ;0.4940 0.1840 0.5560 0.5; 0.4660 0.6740 0.1880 0.5;0.3010 0.7450 0.9330 0.5];

mode = 1;
r_a1 = [0;0;0];
plot_array = matrix_array(r_a1,0.04,8,8);
array_matrices = [plot_array];
sub_arrays = length(array_matrices);

for array = 1:sub_arrays
    z_array_coord = zeros(length(array_matrices(array).r_prime(1,:)));
    plot3(array_matrices(array).r_prime(1,:),array_matrices(array).r_prime(2,:),...
        z_array_coord,'color',colors(1,:),'linestyle','none','marker','o',...
        'MarkerFaceColor','none','linewidth',1.3,'MarkerSize',9);
    axis square
    xlim([-0.08 0.08])
    ylim([-0.08 0.08])
    zlim([-0.02 0.02])
    hold on
end
for array = 1:sub_arrays
    row_elements = array_matrices(array).row_elements;
    column_elements = array_matrices(array).column_elements;
    phase_amplitudes_OF = zeros(1,row_elements*column_elements);
    row_lim = ceil((row_elements)/mode);
    column_lim = ceil((column_elements)/mode);
    for i = 1:row_lim
        for j = 1:column_lim
            element_index = (mode*(i-1))*row_elements + mode*(j-1) +1;
            phase_amplitudes_OF(element_index) = 1;
            plot3(array_matrices(array).r_prime(1,element_index),array_matrices(array).r_prime(2,element_index),...
            0,'color',colors(array,:),'linestyle','none','marker','o',...
            'MarkerFaceColor',colors(array,:),'linewidth',1.3,'MarkerSize',9);
            hold on
        end
    end
end
% e_factor = 0.00006;
% c_factor = e_factor * (uni_a_d-80);
% an = annotation('doublearrow');
% an.X = [0.491-c_factor, 0.568 + c_factor];
% an.Y =  [0.38-c_factor,  0.421 + c_factor];
% 
% dim = [.55 .1 .3 .3];
% str = 'd_a';
% annotation('textbox',dim,'String','$d$','Interpreter','latex','FontSize',12,'LineStyle','none')
% 


grid on
%axis equal
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xticks([-0.06 -0.04 -0.02 0 0.02 0.04 0.06]);
yticks([-0.06 -0.04 -0.02 0 0.02 0.04 0.06]);
zticks([-0.02 0 0.02])
xlabel('x (m)','Interpreter','latex','FontSize',11);
ylabel('y (m)','Interpreter','latex','FontSize',11)
title(strcat('Configuration mode ' ,int2str(mode)),'Interpreter','latex','FontSize',15);
%exportgraphics(gcf,'mode 7.png','Resolution',300)

%% Temporary plotting
max_directivity = 7.5;
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840];

figure(3)
plot(theta/pi*180,AF_8(:,1),'color',colors(1,:),'linewidth',1);
hold on
plot(theta/pi*180,AF_4(:,1),'color',colors(2,:),'linewidth',1);
hold on
plot(theta/pi*180,AF_2(:,1),'color',colors(3,:),'linewidth',1);
hold on
plot(theta/pi*180,AF_1(:,1),'color',colors(4,:),'linewidth',1);
hold on
plot(theta/pi*180,AF_0_5(:,1),'color',colors(5,:),'linewidth',1);
hold on
plot(theta/pi*180,AF_0_25(:,1),'color',colors(6,:),'linewidth',1);
hold on
plot(theta/pi*180,AF_0_125(:,1),'color',colors(7,:),'linewidth',1);
hold on
plot(-theta/pi*180,AF_8(:,1),'color',colors(1,:),'linewidth',1);
hold on
plot(-theta/pi*180,AF_4(:,1),'color',colors(2,:),'linewidth',1);
hold on
plot(-theta/pi*180,AF_2(:,1),'color',colors(3,:),'linewidth',1);
hold on
plot(-theta/pi*180,AF_1(:,1),'color',colors(4,:),'linewidth',1);
hold on
plot(-theta/pi*180,AF_0_5(:,1),'color',colors(5,:),'linewidth',1);
hold on
plot(-theta/pi*180,AF_0_25(:,1),'color',colors(6,:),'linewidth',1);
hold on
plot(-theta/pi*180,AF_0_125(:,1),'color',colors(7,:),'linewidth',1);
xlim([-90 90])
ylim([-5 max_directivity])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\theta$ (deg)','Interpreter','latex','FontSize',18);
ylabel('$|AF(\theta,0)|^2$ (dBi)','Interpreter','latex','FontSize',18)
legend('config mode 1','config mode 2','config mode 3','config mode 4','config mode 5',...
    'config mode 6','config mode 7','Interpreter','latex','location','southeast','FontSize',11)
xticks([-90 -60 -30 0 30 60 90]);
yticks([-5 -2.5 0 2.5 5 7.5])
%yticks([-20 -10 0 5 10 15 20])
grid on
set(gca,'GridALpha',0.1,'LineWidth',.3);
title('8x8 array @ 1.2 kHz','Interpreter','latex','FontSize',15);
exportgraphics(gcf,'1200hz_config_modes_v2.png','Resolution',300)
