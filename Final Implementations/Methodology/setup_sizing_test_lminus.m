%% Setup sizing
%Code written by Bruno Muzy Hegedus

clc
clear

%% Initial Values

lambda = 632.8e-9;
L_FOV_obj =20e-3;

%sCMOS Camera
px_size = 6.5e-6;
L_SENS = 2160*px_size;
Px_pitch_obj = L_SENS/(L_FOV_obj*px_size); %[px per m]

SR = 432; %L_sens/smin

%smin_max = L_SENS/SR; %This gives the size in the image plane of the smallest distinguishable feature, in mm.
%smin_max_px = round(L_SENS/(px_size*SR));
%smin_max = L_FOV_obj/SR; %This gives the size in the object plane of the smallest distinguishable feature, in mm. This size varies in Pixels with Px Pitch
smin_max = 0.0004;


%%Bobcat camera
% px_size = 4.4e-6;
% L_SENS = 1236*px_size;


%px_size = 13e-6;
%dav = 0.2e-3;
%delta_s = 3*px_size;

f_num_data = [11 16 22];
f_data = [0.035 0.04:0.01:0.100 0.105 0.11:0.01:0.2];
l_data = [0.005 0.01:0.01:0.15];

M_data = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
LFOV_data = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
S_data = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
dav_data = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
Z_data = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
df_data = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
m_data = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
CoC_data = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
CoC_obj_data = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
smin_data = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
Px_pitch = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
%Px_pitch_obj = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));
speckle_size = zeros(size(l_data,2), size(f_data, 2),size(f_num_data,2));

%% Solving the system
for k = 1:size(f_num_data,2)
    f_num = f_num_data(k);
    for j = 1:size(f_data, 2)
        f = f_data(j);
        for i = 1:size(l_data,2)
            l = l_data(i);

            syms  delta_s dav S Z df M CoC
            assume(df > 0)
            assume(Z > 0)
            %assume (M > 0) 

            %Delta_t = 35*px_size; %resolution defined by what k pixel sees

            eq1= S == l*M; %Sensitivity
            eq2= delta_s == 1.22*lambda*f_num*(M+1); %Speckle size in the image plane
            eq3= dav == l*M/(f_num*(1+M)) + CoC/M*(1+l/df); %Spatial resolution
            eq4 = 1/f == 1/(Z) + 1/df;
            eq5 = M == Z/df;
            eq6 = CoC == l*f*M/(f_num*(df + l)); %blur
            eq7 = L_FOV_obj == L_SENS*(df+l)/(f*(M+1));
            %eq7 = f == L_SENS*(df+l)/(L_FOV_obj + L_SENS);
            %eq7 = L_FOV_obj == l*M/(f_num*(1+M)) + L_SENS/M*(1+l/df);
            
            
            sol = struct2array(solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7], [S, dav, delta_s, Z, df, M, CoC]));
            
            if size(sol,1) ~= 0 %Only write the results if there is a solution
                S_data(i,j,k) = double(sol(1));
                dav_data(i,j,k) = double(sol(2));
                delta_s = double(sol(3));
    
                M_data(i,j,k) = double(sol(6));
                %f = f_data(i,j,k);
    
                Z = double(sol(4));
                Z_data(i,j,k) = Z;
                
    
                df_data(i,j,k) = double(sol(5));
                df = df_data(i,j,k);

                CoC_data(i,j,k) = double(sol(7));

                m = df + l;
                m_data(i,j,k) = m;

                CoC_obj_data(i,j,k) = CoC_data(i,j,k).*m./Z;

                smin_data(i,j,k) = CoC_obj_data(i,j,k)/2;
                
                % Focal length & object distance
                %f_num = f_num_data(i);
                % m_data(i,j,k) = f*(L_SENS + L_FOV)/L_SENS;
                % m = m_data(i,j,k);
                % D_ap(i,j,k) = f/f_num;
                % 
                % CoC_data(i,j,k) = l*f^2/(f_num*m*(df - f));
                % " The minimum flow feature size to be still resolvable by Bos,it is found to correspond to approximately 50% of the CoC size in the object area."
                
                D_ap(i,j,k) = f/f_num;

                speckle_size(i,j,k) = delta_s/px_size;

                LFOV_data(i,j,k) = L_SENS/M_data(i,j,k);

                Px_pitch(i,j,k) = M_data(i,j,k)/px_size;
                %Px_pitch_obj(i,j,k) = (Z/m)/px_size;
            end

        end
    
    end

end
%% Plotting
[X, Y] = meshgrid(l_data, f_data);
% M_data(abs(M_data) > 3) = 0;
% S_data(abs(S_data) > 0.2) = 0;
% Z_data(abs(Z_data) > 2) = 0;

% lb = 0.01;   % lower bounds
% ub = 0.05;  % upper bounds
% 
% P1 = lhsdesign(3, 1);
% P1_scaled = lb + (ub - lb) .* P1;
% 
% test_points1 = [lb, P1_scaled', ub];
% 
% P2 = lhsdesign(5, 1);
% P2_scaled = lb + (ub - lb) .* P2;
% 
% test_points2 = [P2_scaled'];
% 
% P3 = lhsdesign(3, 1);
% P3_scaled = lb + (ub - lb) .* P3;
% 
% test_points3 = [lb, P3_scaled', ub];


for k = 1:size(f_num_data,2)

    % P1 = lhsdesign(3, 1);
    % P1_scaled = lb + (ub - lb) .* P1;
    % 
    % test_points1 = [lb, P1_scaled', ub];
    % 
    % P2 = lhsdesign(5, 1);
    % P2_scaled = lb + (ub - lb) .* P2;
    % 
    % test_points2 = [P2_scaled'];
    % 
    % P3 = lhsdesign(3, 1);
    % P3_scaled = lb + (ub - lb) .* P3;
    % 
    % test_points3 = [lb, P3_scaled', ub];

    figure(1)
    subplot(size(f_num_data,2),2, 1+2*(k-1))
    contourf(X',Y',S_data(:,:,k), 30)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('Sensitivity, f-stop = %.1f', f_num_data(k)));
    cb1 = colorbar;
    clim([0 quantile(S_data(:), 0.75)])
    cb1.Label.String = 'Sensitivity';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);
    %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
    % %plot(test_points(:,1),test_points(:,2), 'r*')
    % plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
    % plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
    % plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')
    
    subplot(size(f_num_data,2),2, 2+2*(k-1))
    contourf(X',Y',smin_data(:,:,k).*1e3, 30)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('smin, f-stop = %.1f', f_num_data(k)));
    cb1 = colorbar;
    clim([0 quantile(smin_data(:).*1e3,0.75)])
    cb1.Label.String = 'smin [mm]';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);
    %contour(X', Y', CoC_data(:,:,k)./px_size, [1 1], 'g', 'LineWidth', 2)
    %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);   
    %plot(test_points(:,1),test_points(:,2), 'r*')    
    % plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
    % plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
    % plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')

    % subplot(3,3,3 + 3*(k-1))
    % contourf(X',Y', LFOV_data(:,:,k), 20)
    % xlabel('l [m]');
    % ylabel('f');
    % title(sprintf('L FOV for f-stop = %d', f_num_data(k)));
    % cb1 = colorbar;
    % clim([0 1])
    % cb1.Label.String = 'L FOV [m]';
    % hold on
    % contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2);
    % contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
    % %plot(test_points(:,1),test_points(:,2), 'r*')
    % plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
    % plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
    % plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')


    figure(2)
    subplot(size(f_num_data,2),1, k)
    contourf(X',Y',speckle_size(:,:,k), 20)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('speckle size, f-stop = %.1f', f_num_data(k)));
    cb1 = colorbar;
    clim([0 quantile(speckle_size(:),0.9)])
    cb1.Label.String = 'speckle size [px]';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);
    %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
    %plot(test_points(:,1),test_points(:,2), 'r*')
    % plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
    % plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
    % plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')

    figure(3)
    subplot(size(f_num_data,2),3,1 + 3*(k-1))
    contourf(X',Y', M_data(:,:,k), 50)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('M for f-stop = %d', f_num_data(k)));
    cb1 = colorbar;
    clim([0 quantile(M_data(:),0.75)])
    cb1.Label.String = 'M [m]';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);
    %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
    % %plot(test_points(:,1),test_points(:,2), 'r*')
    % plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
    % plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
    % plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')

    subplot(size(f_num_data,2),3, 2 + 3*(k-1))
    contourf(X',Y', df_data(:,:,k), 20)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('df for f-stop = %d', f_num_data(k)));
    cb1 = colorbar;
    clim([0 quantile(df_data(:),0.75)])
    cb1.Label.String = 'df [m]';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);
    %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
    %end

    subplot(size(f_num_data,2),3, 3*(k-1) + 3)
    contourf(X',Y', m_data(:,:,k), 20)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('m for f-stop = %d', f_num_data(k)));
    cb1 = colorbar;
    clim([0 quantile(m_data(:),0.75)])
    cb1.Label.String = 'm [m]';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);

    
    % figure(4)
    % subplot(size(f_num_data,2),1,k)
    % contourf(X',Y',Z_data(:,:,k), 10)
    % xlabel('l [m]');
    % ylabel('f');
    % title(sprintf('Z, f-stop = %.1f', f_num_data(k)));
    % cb1 = colorbar;
    % clim([0 max(Z_data)])
    % cb1.Label.String = 'Z [m]';
    

    figure(5)
    subplot(size(f_num_data,2),1, k)
    contourf(X',Y',CoC_obj_data(:,:,k), 50)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('CoC obj, f-stop = %.1f', f_num_data(k)));
    cb1 = colorbar;
    clim([0 quantile(CoC_obj_data(:),0.75)])
    cb1.Label.String = 'CoC obj[m]';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);

    figure(6)
    subplot(size(f_num_data,2),1, k)
    contourf(X',Y',smin_data(:,:,k).*Px_pitch_obj, 50)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('smin, f-stop = %.1f', f_num_data(k)));
    cb1 = colorbar;
    clim([0 quantile(smin_data(:).*Px_pitch_obj,0.75)])
    cb1.Label.String = 'smin [px]';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);

    figure(8)
    subplot(size(f_num_data,2),1, k)
    contourf(X',Y',dav_data(:,:,k).*1e3, 20)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('dav, f-stop = %.1f', f_num_data(k)));
    cb1 = colorbar;
    clim([0 quantile(dav_data(:).*1e3,0.75)])
    cb1.Label.String = 'dav [mm]';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);

        
    figure(7)
    subplot(size(f_num_data,2),2,1 + 2*(k-1))
    contourf(X',Y', Px_pitch(:,:,k).*10^-3, 50)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('Px Pitch for f-stop = %d', f_num_data(k)));
    cb1 = colorbar;
    clim([0 quantile(Px_pitch(:).*10^-3,0.75)])
    cb1.Label.String = 'Px pitch [1/mm]';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);
    
    %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);


    %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
    % %plot(test_points(:,1),test_points(:,2), 'r*')
    % plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
    % plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
    % plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')
    
    subplot(size(f_num_data,2),2,2 + 2*(k-1))
    contourf(X',Y', LFOV_data(:,:,k), 20)
    xlabel('l [m]');
    ylabel('f');
    title(sprintf('LFOV for f-stop = %d', f_num_data(k)));
    cb1 = colorbar;
    clim([0 L_FOV_obj])
    cb1.Label.String = 'LFOV [m]';
    hold on
    contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    contour(X', Y', smin_data(:,:,k), [smin_max smin_max], 'r', 'LineWidth', 2);
    %contour(X', Y', smin_data(:,:,k).*Px_pitch_obj, [smin_max_px smin_max_px], 'green', 'LineWidth', 2);


    % if k == 2
    % 
    %     figure(13)
    %     subplot(1,2, 1)
    %     contourf(X',Y',S_data(:,:,k), 20)
    %     xlabel('l [m]');
    %     ylabel('f');
    %     title(sprintf('Sensitivity, f-stop = %.1f', f_num_data(k)));
    %     cb1 = colorbar;
    %     clim([0 0.08])
    %     cb1.Label.String = 'Sensitivity';
    %     hold on
    %     contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2);
    %     contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    %     %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
    %     % %plot(test_points(:,1),test_points(:,2), 'r*')
    %     % plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
    %     % plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
    %     % plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')
    % 
    %     subplot(1,2, 2)
    %     contourf(X',Y',dav_data(:,:,k), 20)
    %     xlabel('l [m]');
    %     ylabel('f');
    %     title(sprintf('dav, f-stop = %.1f', f_num_data(k)));
    %     cb1 = colorbar;
    %     clim([0 0.003])
    %     cb1.Label.String = 'dav [m]';
    %     hold on
    %     contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2);
    %     contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    %     %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);   
    %     % %plot(test_points(:,1),test_points(:,2), 'r*')    
    %     % plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
    %     % plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
    %     % plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')
    % 
    %     % subplot(1,3,3)
    %     % contourf(X',Y', LFOV_data(:,:,k), 200)
    %     % xlabel('l [m]');
    %     % ylabel('f');
    %     % title(sprintf('LFOV for f-stop = %d', f_num_data(k)));
    %     % cb1 = colorbar;
    %     % clim([0 1])
    %     % cb1.Label.String = 'L FOV [m]';
    %     % hold on
    %     % contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2);
    %     % contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    % 
    % 
    %     figure(15)
    %     subplot(1,3,1)
    %     contourf(X',Y', M_data(:,:,k), 20)
    %     xlabel('l [m]');
    %     ylabel('f');
    %     title(sprintf('M for f-stop = %d', f_num_data(k)));
    %     cb1 = colorbar;
    %     clim([0 1])
    %     cb1.Label.String = 'M [m]';
    %     hold on
    %     contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2);
    %     contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    %     %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
    %     % %plot(test_points(:,1),test_points(:,2), 'r*')
    %     % plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
    %     % plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
    %     % plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')
    % 
    %     subplot(1,3, 2)
    %     contourf(X',Y', Z_data(:,:,k), 20)
    %     xlabel('l [m]');
    %     ylabel('f');
    %     title(sprintf('Z for f-stop = %d', f_num_data(k)));
    %     cb1 = colorbar;
    %     clim([0 0.4])
    %     cb1.Label.String = 'Z [m]';
    %     hold on
    %     contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2);
    %     contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    %     %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
    %     %end
    % 
    %     subplot(1,3,  3)
    %     contourf(X',Y', df_data(:,:,k), 20)
    %     xlabel('l [m]');
    %     ylabel('f');
    %     title(sprintf('df for f-stop = %d', f_num_data(k)));
    %     cb1 = colorbar;
    %     %clim([0 20])
    %     cb1.Label.String = 'df [m]';
    %     hold on
    %     contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2)
    %     contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
    %     %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
    % 
    %     figure(16)
    %     contourf(X',Y',m_data(:,:,k), 10)
    %     xlabel('l [m]');
    %     ylabel('f');
    %     title(sprintf('m, f-stop = %.1f', f_num_data(k)));
    %     cb1 = colorbar;
    %     clim([0 0.5])
    %     cb1.Label.String = 'm [m]';
    % end

end


% 
% % lb = 0.01;   % lower bounds
% % ub = 0.05;  % upper bounds
% % 
% % P1 = lhsdesign(8, 1);
% % P1_scaled = lb + (ub - lb) .* P1;
% % 
% % test_points1 = [lb, P1_scaled', ub];
% % 
% % P2 = lhsdesign(5, 1);
% % P2_scaled = lb + (ub - lb) .* P2;
% % 
% % test_points2 = [P2_scaled'];
% % 
% % P3 = lhsdesign(8, 1);
% % P3_scaled = lb + (ub - lb) .* P3;
% % 
% % test_points3 = [lb, P3_scaled', ub];
% 
% figure(13)
% subplot(1,2, 1)
% contourf(X',Y',S_data(:,:,k), 200)
% xlabel('l [m]');
% ylabel('f');
% title(sprintf('Sensitivity, f-stop = %.1f', f_num_data(k)));
% cb1 = colorbar;
% clim([0 0.08])
% cb1.Label.String = 'Sensitivity';
% hold on
% contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2);
% contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
% %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
% %plot(test_points(:,1),test_points(:,2), 'r*')
% plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
% plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
% plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')
% 
% subplot(1,2, 2)
% contourf(X',Y',dav_data(:,:,k), 20)
% xlabel('l [m]');
% ylabel('f');
% title(sprintf('dav, f-stop = %.1f', f_num_data(k)));
% cb1 = colorbar;
% clim([0 0.003])
% cb1.Label.String = 'dav [m]';
% hold on
% contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2);
% contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
% %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);   
% %plot(test_points(:,1),test_points(:,2), 'r*')    
% plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
% plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
% plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')
% 
% 
% figure(15)
% subplot(1,3,1)
% contourf(X',Y', M_data(:,:,k), 200)
% xlabel('l [m]');
% ylabel('f');
% title(sprintf('M for f-stop = %d', f_num_data(k)));
% cb1 = colorbar;
% clim([-0.8 0])
% cb1.Label.String = 'M [m]';
% hold on
% contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2);
% contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
% %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
% %plot(test_points(:,1),test_points(:,2), 'r*')
% plot(test_points1, 0.1*ones(size(test_points1,2)), 'r*')
% plot(test_points2, 0.15*ones(size(test_points2,2)), 'r*')
% plot(test_points3, 0.2*ones(size(test_points3,2)), 'r*')
% 
% subplot(1,3, 2)
% contourf(X',Y', Z_data(:,:,k), 20)
% xlabel('l [m]');
% ylabel('f');
% title(sprintf('Z for f-stop = %d', f_num_data(k)));
% cb1 = colorbar;
% clim([0 0.4])
% cb1.Label.String = 'Z [m]';
% hold on
% contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2);
% contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
% %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
% %end
% 
% subplot(1,3,  3)
% contourf(X',Y', df_data(:,:,k), 20)
% xlabel('l [m]');
% ylabel('f');
% title(sprintf('df for f-stop = %d', f_num_data(k)));
% cb1 = colorbar;
% %clim([0 20])
% cb1.Label.String = 'df [m]';
% hold on
% contour(X', Y', dav_data(:,:,k), [0.0005 0.0005], 'r', 'LineWidth', 2)
% contour(X', Y', S_data(:,:,k), [0.02 0.02], 'white', 'LineWidth', 2);
% %contour(X', Y', speckle_size(:,:,k), [3 3], 'g', 'LineWidth', 2);
% 
% figure(16)
% hold on
% plot(f_data, m_data(1, :, k), 'DisplayName', sprintf('f-stop = %d', f_num_data(k)))
% xlabel('f')
% ylabel('m')
% legend show;