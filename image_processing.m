% Developed by Vineet Mukim, Research Fellow, University of Stavanger
clc; clear all;
close all;
warning off;
tic;
%gpuDevice; % just to initiate GPU usage

calibration_flag = 1; % 0 for OFF, 1 for ON
if calibration_flag
%     clear all;
    disp('calibrating camera');
    %--------------------------------------------------------------------------
    % Camera Calibration with option on or off, or maybe do this everytime
    % calibration images to process
    calibrationImages = imageDatastore('C:\Users\mukim\Documents\PhD\Work\Projects\Capillary Wave Surface Reconstruction\Journal Paper 2\experimental_data\fluid1\calibration');
    % Detect checkerboards in images
    [image_points, board_size, images_used] = detectCheckerboardPoints(calibrationImages.Files);
    % selecting the last image in folder list as chosen view and reading it
    chosen_camera_view = imread(calibrationImages.Files{end});
    disp('final orientation of calibeation pattern == ');
    disp(calibrationImages.Files{end});% fprintf('\n\n');    
    [mrows, ncols] = size(chosen_camera_view);
%     figure; imshow(chosen_camera_view);
    % Generate world coordinates of the corners of the squares
    square_size = 4.8;%10.0;  % in units of 'millimeters'
    world_points = generateCheckerboardPoints(board_size, square_size);
    % Calibrate the camera
    camera_params = estimateCameraParameters(image_points, world_points, 'ImageSize', [mrows, ncols]);
    % View reprojection errors
%     h1=figure; showReprojectionErrors(camera_params);
    % Visualize pattern locations
%     h2=figure; showExtrinsics(camera_params, 'PatternCentric');
    % Display parameter estimation errors
%     displayErrors(estimationErrors, cameraParams);
    % undistortion
    % calibration folder has many images, choosing the one view which is
    % used to capture the images for undistortion and calculating
    % extrinsics of the camera
    [undistorted_chosen_camera_view, new_origin] = undistortImage(chosen_camera_view, camera_params);
    %figure; imshow(undistorted_chosen_camera_view);    
    [image_points, ~] = detectCheckerboardPoints(undistorted_chosen_camera_view);
    % Adjust the imagePoints so that they are expressed in the coordinate system
    % used in the original image, before it was undistorted.  This adjustment
    % makes it compatible with the cameraParameters object computed for the original image.
    image_points = image_points + new_origin;
    % Extract camera intrinsics.
    cam_intrinsics = camera_params.IntrinsicMatrix;
    % Compute extrinsic parameters of the camera.
    [rotation_mat, translation_vec] = extrinsics(image_points, world_points, camera_params);
    % just trying to use the mapping function on the known data to verify
    % if it works, actual use of this function should be done after mean
    % finding algorithm on actual images
    mapped_world_points = pointsToWorld(camera_params, rotation_mat, translation_vec, image_points);
    mapping_error = abs(mapped_world_points - world_points);
end

%--------------------------------------------------------------------------

% % analyze all folders in a given path for various freq levels in for loop
% % and calculate the mean curve and save all files
folderPath = 'C:\Users\mukim\Documents\PhD\Work\Projects\Capillary Wave Surface Reconstruction\Journal Paper 2\experimental_data'; 
fluidList = ["fluid1"];%, "fluid2", "fluid3"];
threshold_filter = 50;%100; % intensity threshold
std_dev_gaussian_filter = 1; % standard deviation for gaussian filter
file_string_length = 126;
crop_box = [2, 640, 10, 720]; %[xmin xmax ymin ymax]
overwrite_flag = 1;

freqList = 70:100:100;%[10, 20, 50]; % cant start from 0 as zero freq treated separately
fluid_depth = 10;%13;
translate_x = -43; % distance in x axis of pattern origin and center (point and excitation)
translate_y = 26; % distance in y axis of pattern origin and center (point and excitation)

num_screen_intercept = 100; % downsizing screen intercept after interpolation fit
num_surface_intercept = 800; % num of points on surface has to be more than screen for more accurate calculation

% THESE COORIDINATES ARE WITH REFERENCE TO THE PETRY DISH CENTER AND NOT
% THE CALIBRATION PATTERN ORIGIN
z_sc = 592; % check how to measure this..wrt glass bottom or top etc
% x_laser = 43; % caluclated from zero freeq line fit to make sure this point lies on the zero freq line
y_laser = -3;
z_laser = 616; % check how to measure this..wrt glass bottom or top or liquid surface etc
% translated_light_source = [x_laser, y_laser, z_laser];

glass_th = 2.40;
mu_liquid_air = 1.33;
mu_glass_air = 1.50;

for i=1:length(fluidList)
    fluid = fluidList(i);
    
    % -------------------------------------------
    % zero frequency calculation
    folder_zero_freq = fullfile(folderPath, fluid, '00');
    imageData = imageDatastore(folder_zero_freq);
    waveImageList_zero_freq = imageData.Files;
    cd(folder_zero_freq);
    if ~isfolder('processed')
        disp('folder created in zero freq');
        mkdir('processed');
    end        
    cd(fullfile(folder_zero_freq, 'processed'));
    img_path_zero_freq = waveImageList_zero_freq{1};
   
    % zero freq line plot
    [screen_intercept_zero_freq] = processing(img_path_zero_freq, crop_box, std_dev_gaussian_filter,...
        threshold_filter, camera_params, overwrite_flag, file_string_length);
    zero_freq_flag = 1;
    smoothening_para = 0.02;
    [mapped_translated_screen_intercept_downsized, fitresult_yx_zero_freq] = coordinate_transform(camera_params, rotation_mat,...
        translation_vec, screen_intercept_zero_freq, translate_x, translate_y, num_surface_intercept, zero_freq_flag, smoothening_para);    

    x_laser = fitresult_yx_zero_freq(y_laser); % to make sure laser source exctly above the illumination, psi=0
    translated_light_source = [x_laser, y_laser, z_laser]; % this has to be laser source location correspnsing to new origin i.e. surface of center of petry dish
    scatter(x_laser, y_laser, 'o', 'DisplayName', "laser source");
    hold on;

    % illumination might not be parallel to y axis, to account for slant illumination, the cooridnates are calculated by similarity of triangles
    % calculating x and y intercept coordinates of illuminated line (most likely to be slant)
    x_intercept_list = (mapped_translated_screen_intercept_downsized(:, 1) - x_laser)*z_laser/(z_laser+z_sc) + x_laser;
    y_intercept_list = (mapped_translated_screen_intercept_downsized(:, 2) - y_laser)*z_laser/(z_laser+z_sc) + y_laser;

    scatter(x_intercept_list, y_intercept_list, 'b.', 'DisplayName', "downsized data surface intercept");
    xlim([mean(x_intercept_list)-3, mean(x_intercept_list)+3]);
    print(gcf, strcat(folder_zero_freq, '\processed\', img_path_zero_freq(157:end-4), '_mapped_translated_screen_intercept_zero_freq.png'), '-dpng');

    for j=1:length(freqList)   
        freq = num2str(freqList(j));
        folder = fullfile(folderPath, fluid, freq);
        imageData = imageDatastore(folder);
        waveImageList = imageData.Files;
        cd(folder);
        if ~isfolder('processed')
            disp('folder created');
            mkdir('processed');
        end        
        cd(fullfile(folder, 'processed'));
        disp('current working directory is == '); disp(cd);
        lambda_list = [];
        for l=5:5%1:length(waveImageList)
            img_path = waveImageList{l};
            disp(img_path(157:end-4));
            [screen_intercept] = processing(img_path, crop_box, std_dev_gaussian_filter,...
                threshold_filter, camera_params, overwrite_flag, file_string_length);
            if isempty(screen_intercept)
                % empty array returned if overwrite found
                continue;
            end
            cd('C:\Users\mukim\Documents\PhD\Work\Projects\Capillary Wave Surface Reconstruction\Code\');

% %             figure;
% %             scatter(screen_intercept(:, 1), screen_intercept(:, 2), 'b.');
% %             xlabel('x axis'); ylabel('y axis');
% %             title('refracted laser sheet on screen in image coordinates');

            %--------------------------------------------------------------            
            % inverse transform to get world coordinates from image
            % coordinates, use camera calinration and point2world functions
%             disp('image to world coordinate transfomation');            
            zero_freq_flag = 0;
            smoothening_para = 0.01;
            [mapped_translated_screen_intercept_downsized, ~] = coordinate_transform(camera_params, rotation_mat,...
                translation_vec, screen_intercept, translate_x, translate_y, num_screen_intercept, zero_freq_flag, smoothening_para);          

            % can plot actual, scaled, fitted, translated and mapped refracted screen intercept data separately
            scatter(mapped_translated_screen_intercept_downsized(:, 1), mapped_translated_screen_intercept_downsized(:, 2), 'bo', 'DisplayName', "downsized data screen intercept");
            xlim([mean(mapped_translated_screen_intercept_downsized(:, 1))-3, mean(mapped_translated_screen_intercept_downsized(:, 1))+3]);
            print(gcf, strcat(folder, '\processed\', img_path(157:end-4), '_mapped_translated_screen_intercept.png'), '-dpng');     
           
%             --------------------------------------------------------------
%             run inverse algorithm on these world cooridnates, so many input paramters needed and obtain the surface cross-section

%             disp('inverse algorithm');
            [surface_intercept_calculated, screen_vec, glass_top_intercept_inv, glass_bottom_intercept_inv,...
                unit_surface_normal_calc_final_list, surface_normal_coplanar_error, snell_error, dF_dy_list_non_smooth, dF_dy_list_smooth]...
                = inverse_solver(x_intercept_list, y_intercept_list, z_sc, fluid_depth, 0, glass_th,...
                mu_liquid_air, mu_glass_air, translated_light_source, mapped_translated_screen_intercept_downsized, 0.2);                       

%             % plot screen and calculated top and bottom glass intercept
%             % with subfigures
%             fig1 = figure;
%             ax = subplot(3,1,1, 'Parent', fig1);
%             for bb=1:length(screen_vec)
%                 scatter(ax, glass_top_intercept_inv(2, bb), glass_top_intercept_inv(1, bb));
%                 hold (ax, 'on');
%             end
%             title(ax,'Calculated Glass Top Intercept');
%             xlabel(ax, 'y-axis (mm)'); ylabel(ax, 'x-axis (mm)');
%             ax = subplot(3,1,2, 'Parent', fig1);
%             for bb=1:length(screen_vec)
%                 scatter(ax, glass_bottom_intercept_inv(2, bb), glass_bottom_intercept_inv(1, bb));
%                 hold (ax, 'on');
%             end
%             title(ax, 'Calculated Glass Bottom Intercept');
%             xlabel(ax, 'y-axis (mm)'); ylabel(ax, 'x-axis (mm)');
%             ax = subplot(3,1,3, 'Parent', fig1);
%             for bb=1:length(screen_vec)
%                 scatter(ax, screen_vec(2, bb), screen_vec(1, bb));
%                 hold (ax, 'on');
%             end
%             title(ax, 'Mapped Translated Screen Intercept');
%             xlabel(ax, 'y-axis (mm)'); ylabel(ax, 'x-axis (mm)');
%             print(fig1, strcat(folder, '\processed\', img_path(157:end-4), '_intercepts.png'), '-dpng');
            
            %-----------
      
            % without subfigure
            figure;
            for bb=1:length(screen_vec)
                scatter(glass_top_intercept_inv(2, bb), glass_top_intercept_inv(1, bb));
                hold on;
            end
%             title('Calculated Glass Top Intercept');
            xlabel('y-axis (mm)'); ylabel('x-axis (mm)');
            pbaspect([4 1 1]);
            print(gcf, strcat(folder, '\processed\', img_path(157:end-4), '_calc_glass_top_intercepts.png'), '-dpng');
            
            figure;
            for bb=1:length(screen_vec)
                scatter(glass_bottom_intercept_inv(2, bb), glass_bottom_intercept_inv(1, bb));
                hold on;
            end
%             title('Calculated Glass Bottom Intercept');
            xlabel('y-axis (mm)'); ylabel('x-axis (mm)');
            pbaspect([4 1 1]);
            print(gcf, strcat(folder, '\processed\', img_path(157:end-4), '_calc_glass_bottom_intercepts.png'), '-dpng');
            
            figure;
            for bb=1:length(screen_vec)
                scatter(screen_vec(2, bb), screen_vec(1, bb));
                hold on;
            end
%             title('Mapped Translated Screen Intercept');
            xlabel('y-axis (mm)'); ylabel('x-axis (mm)');
            pbaspect([4 1 1]);
            print(gcf, strcat(folder, '\processed\', img_path(157:end-4), '_screen_intercepts.png'), '-dpng');
            
            %---------------------------------------------------------------------------------         
            % plot screen to source ray tracing
            figure;  
            for bb=1:length(screen_vec)
                plot3([screen_vec(1, bb); glass_bottom_intercept_inv(1, bb)], [screen_vec(2, bb); glass_bottom_intercept_inv(2, bb)], [screen_vec(3, bb); glass_bottom_intercept_inv(3, bb)]);
                hold on;  
                plot3([glass_top_intercept_inv(1, bb); glass_bottom_intercept_inv(1, bb)], [glass_top_intercept_inv(2, bb); glass_bottom_intercept_inv(2, bb)], [glass_top_intercept_inv(3, bb); glass_bottom_intercept_inv(3, bb)]);
                hold on;  
                plot3([glass_top_intercept_inv(1, bb); surface_intercept_calculated(1, bb)], [glass_top_intercept_inv(2, bb); surface_intercept_calculated(2, bb)], [glass_top_intercept_inv(3, bb); fluid_depth]);
                hold on;  
                plot3([surface_intercept_calculated(1, bb); translated_light_source(1)], [surface_intercept_calculated(2, bb); translated_light_source(2)], [surface_intercept_calculated(3, bb); translated_light_source(3)]); 
                hold on;  
            end
            xlabel('x-axis (mm)'); ylabel('y-axis (mm)'); zlabel('z-axis (mm)');
%             title('Inverse Ray Tracing');
            view(45, 5);
            print(gcf, strcat(folder, '\processed\', img_path(157:end-4), '_inverse_ray_tracing.png'), '-dpng');    

            %---------------------------------------------------------------------------------
%             % plot incident, refracted vectors and unit surface normals
%             figure;
%             surface_intercept_calculated_mod = surface_intercept_calculated;
%             surface_intercept_calculated_mod(3, :) = fluid_depth;
%             final_incident_vec = surface_intercept_calculated_mod-glass_top_intercept_inv;
%             unit_final_incident_vec = final_incident_vec / norm(final_incident_vec);
%             final_refracted_vec = translated_light_source'.*ones(3, length(screen_vec)) - surface_intercept_calculated_mod;
%             unit_final_refracted_vec = final_refracted_vec / norm(final_refracted_vec);
%             translated_unit_refracted_vec = 0.1*unit_final_refracted_vec + surface_intercept_calculated_mod;    
%             %scaling only the y compoenent of unit surface normal for visualization
% %             unit_surface_normal_calc_final_list(:, i) = [0; 10*unit_surface_normal_calc_final_list(2, i); 0];
% %             disp('y coordinate non zero'); disp(unit_surface_normal_calc_final_list(:, i));
%             translated_unit_surface_normal_vec = 0.1*unit_surface_normal_calc_final_list + surface_intercept_calculated_mod;
%             translated_unit_incident_vec = -0.1*unit_final_incident_vec + surface_intercept_calculated_mod;
%             plot3([surface_intercept_calculated_mod(1, :); translated_unit_refracted_vec(1, :)], [surface_intercept_calculated_mod(2, :); translated_unit_refracted_vec(2, :)], [surface_intercept_calculated_mod(3, :); translated_unit_refracted_vec(3, :)], 'm');
%             hold on; 
%             plot3([surface_intercept_calculated_mod(1, :); translated_unit_incident_vec(1, :)], [surface_intercept_calculated_mod(2, :); translated_unit_incident_vec(2, :)], [surface_intercept_calculated_mod(3, :); translated_unit_incident_vec(3, :)], 'b');
%             hold on; 
%             plot3([surface_intercept_calculated_mod(1, :); translated_unit_surface_normal_vec(1, :)], [surface_intercept_calculated_mod(2, :); translated_unit_surface_normal_vec(2, :)], [surface_intercept_calculated_mod(3, :); translated_unit_surface_normal_vec(3, :)], 'k');
%             xlabel('x-axis (mm)'); ylabel('y-axis (mm)'); zlabel('z-axis (mm)');    

            %---------------------------------------------------------------------------------
            % plotting slope (dF/dy)
            figure;
            plot(surface_intercept_calculated(2, :), dF_dy_list_non_smooth);
            hold on;            
            plot(surface_intercept_calculated(2, :), dF_dy_list_smooth);
            xlabel('y-axis (mm)');
            ylabel('slope (-)');
            legend('calculated', 'smoothened');            
%             title('Surface Topography Shape Slope');
            pbaspect([4 1 1]);
            set(gca,'box','off'); legend('Location','northoutside', 'Orientation','horizontal');
            print(gcf, strcat(folder, '\processed\', img_path(157:end-4), '_slope.png'), '-dpng');                           

            %---------------------------------------------------------------------------------
            % plot reconstructed cross sections
            figure;
            scatter(surface_intercept_calculated(2, :), surface_intercept_calculated(3,:));
            hold on;
            ylabel('z-axis (mm)');
            xlabel('y-axis (mm)');
%             title('Surface Topography Cross-Section');
            pbaspect([4 1 1]);
            print(strcat(folder, '\processed\', img_path(157:end-4), '_cross_section.png'), '-dpng');
        %     filename = [image_file, msg, 'cross-section.png']; 
%             set(gcf, 'FontSize', 10);%,'FontName','Times');
        % %     saveas(gcf, filename, 'epsc');
        %     print(gcf, filename,'-depsc2','-r600');
        %     print(gcf, filename,'-dpng');

            %--------------------------------------------------------------
            % construct the 3 D surace and also obtain a radial cross section from this     
%             disp('3D surface reconstruction');
            id = find(surface_intercept_calculated(2, :) > 0, 1, 'first');    
            r_list_lhs = sqrt(surface_intercept_calculated(1, 1:id-1).^2+surface_intercept_calculated(2, 1:id-1).^2);
            amp_list_lhs = surface_intercept_calculated(3, 1:id-1);
            r_list_rhs = sqrt(surface_intercept_calculated(1, id:end).^2+surface_intercept_calculated(2, id:end).^2);
            amp_list_rhs = surface_intercept_calculated(3, id:end);

            % findings peaks
            skip_peak_index = 5;
            [~, peak_id_lhs1] = findpeaks(amp_list_lhs(1:end-skip_peak_index));
            [~, peak_id_lhs2] = findpeaks(-amp_list_lhs(1:end-skip_peak_index));            
            peak_id_lhs = unique([peak_id_lhs1 peak_id_lhs2]);
%             peak_id_lhs = peak_id_lhs+skip_peak_index-1;  % dont need this as indexing starts from 1           
    
            [~, peak_id_rhs1] = findpeaks(amp_list_rhs(skip_peak_index:end));
            [~, peak_id_rhs2] = findpeaks(-amp_list_rhs(skip_peak_index:end));
            peak_id_rhs = unique([peak_id_rhs1 peak_id_rhs2]);
            peak_id_rhs = peak_id_rhs+skip_peak_index-1;
            
            disp('peak r lhs = '); disp(r_list_lhs(peak_id_lhs));
            disp('peak r rhs = '); disp(r_list_rhs(peak_id_rhs));

            
                        
            % plot radial cross-section
            figure;
            plot(r_list_lhs, amp_list_lhs, 'b');
            hold on;
            plot(r_list_rhs, amp_list_rhs, 'r');
            hold on;
            scatter(r_list_lhs(peak_id_lhs), amp_list_lhs(peak_id_lhs), 'bo');
            hold on;
            scatter(r_list_rhs(peak_id_rhs), amp_list_rhs(peak_id_rhs), 'ro');
            xlabel('radius (mm)');
            ylabel('amplitude (mm)');
%             title('Reconstructed Radial Cross-Section');
            legend('left zone','right zone'); legend('Location','northoutside', 'Orientation','horizontal');
            pbaspect([4 1 1]);
            set(gca,'box','off');
            print(gcf, strcat(folder, '\processed\', img_path(157:end-4), '_radial_cross_section.png'), '-dpng');

            %--------------------------------------------------------------
            % This part first calculates wave parameters like wavelength, wavespeed and decay coeff and other 
            % parameters and then uses existing eq. to find surface tension and viscoisity

%             disp('wavelength calculation');
            peak_r_list_lhs = r_list_lhs(peak_id_lhs);
            peak_r_list_rhs = r_list_rhs(peak_id_rhs);

            lambda_lhs = abs(peak_r_list_lhs(2:end)-peak_r_list_lhs(1:end-1))*2;
            lambda_rhs = abs(peak_r_list_rhs(2:end)-peak_r_list_rhs(1:end-1))*2; 
            disp('lambda lhs = '); disp(lambda_lhs);
            disp('lambda rhs = '); disp(lambda_rhs);
            
            lambda_list = [lambda_list, lambda_lhs, lambda_rhs];
%             w = waitforbuttonpress;
% %             close all
% %           clear all;

            clear screen_intercept mapped_translated_screen_intercept...
            mapped_translated_screen_intercept_downsized surface_intercept_calculated screen_vec glass_top_intercept_inv...
            glass_bottom_intercept_inv unit_surface_normal_calc_final_list surface_normal_coplanar_error snell_error dF_dy_list_non_smooth...
            dF_dy_list_smooth id r_list_lhs r_list_rhs amp_list_lhs amp_list_rhs peak_id_lhs peak_id_rhs peak_r_list_lhs peak_r_list_rhs...
            lambda_lhs lambda_rhs; 
            close all;
        end  
        disp(lambda_list);
        mean_lambda = mean(lambda_list);
        mean_k = 2*pi/(mean_lambda*1e-3);
        omega_val = 2*pi*str2double(freq);
        rho = 1000;
        g = 9.81;
        sigma_calc = (omega_val.^2/mean_k-g)*rho/mean_k.^2;
        disp(['mean wavelength (all images for ', freq,  'Hz (lhs & rhs avg) = ', num2str(mean_lambda, '%0.2f'), ' mm']);
        disp(['interfacial tension calculated = ', num2str(sigma_calc*1000, '%0.2f'), ' mN/m']);
        disp('------------------');
    end  
%     clear screen_intercept_zero_freq mapped_translated_screen_intercept_downsized_zero_freq_line fitresult_zero_freq_line_yx...
%     x_laser translated_light_source x_intercept_list y_intercept_list;  
end


%--------------------------------------------------------------
% This part first calculates wave parameters like wavelength, wavespeed and decay coeff and other 
% parameters and then uses existing eq. to find surface tension and viscoisity and then plo0t therse 
% values for paperdo some analysis on the radial cross section profiles to obtain wave
% Some other things can be done too, like studing the diffusion of waves
% and verifying energy conservation and building a model for interface


%--------------------------------------------------------------------------
toc;
%--------------------------------------------------------------------------


















%--------------------------------------------------------------------------
% function definitions
%--------------------------------------------------------------------------

function [surface_intercept_calculated, screen_vec, glass_top_intercept_inv, glass_bottom_intercept_inv,...
                unit_surface_normal_calc_final_list, surface_normal_coplanar_error, snell_error, dF_dy_list_non_smooth, dF_dy_list_smooth]...
    = inverse_solver(x_intercept_list, y_intercept_list, z_sc, fluid_depth, bottom_level, glass_th,...
    mu_liquid_air, mu_glass_air, light_source, screen_intercept, smoothening_param)
    
    mean_surface_intercept_z = fluid_depth;
    mean_glass_top_intercept_z = bottom_level;
    mean_glass_bottom_intercept_z = bottom_level - glass_th;
    mean_screen_intercept_z = bottom_level - (z_sc + glass_th);

%     snell_check_1 = zeros(length(screen_intercept));
%     snell_check_2 = zeros(length(screen_intercept));
%     snell_check_3 = zeros(length(screen_intercept));
    
%     glass_top_intercept_inv_temp = 1e6*ones(3, length(screen_intercept));
%     glass_bottom_intercept_inv_temp = 1e6*ones(3, length(screen_intercept));
%     unit_surface_normal_calc_temp_list = 1e6*ones(3, length(screen_intercept));

    unit_surface_normal_calc_final_list = 1e6*ones(3, length(screen_intercept));
    glass_top_intercept_inv = 1e6*ones(3, length(screen_intercept));
    glass_bottom_intercept_inv = 1e6*ones(3, length(screen_intercept));
    surface_normal_coplanar_error = 1e6*ones(length(screen_intercept), 1); 
    surface_intercept_calculated = ones(3, length(screen_intercept));
    surface_intercept_calculated(3, :) = 0; 
    
%     mu_glass_liquid = mu_glass_air/mu_liquid_air;
    source = transpose(light_source);
    
    options_fermat_function = optimset('Display', 'off', 'MaxFunEvals', 5e3, 'MaxIter', 5e3, 'TolFun', 1e-12, 'TolX', 1e-12);
    fermat_guess = [0; 0; 0; 0];  % initial guess for fermat function can be interception of a straight line
      
    options_snell_function = optimoptions('fsolve', 'Display', 'off', 'FinDiffType', 'central',...
                'MaxIter', 4000, 'MaxFunEvals', 2000000, 'FunctionTolerance', 1e-20, 'StepTolerance', 1e-8);
    coeff_guess = [1; 1]; 
%     coeff_guess = [-1; -1];             
    
    for i=1:length(screen_intercept) % loop over point on screen
%         disp(strcat('screen intercept point no:', num2str(i), '/', num2str(length(screen_intercept))));        
        screen_vec = [screen_intercept(i, 1); screen_intercept(i, 2); mean_screen_intercept_z];
%         disp(strcat('screen vec: ', num2str(screen_vec')));

        parfor j=1:length(y_intercept_list)  % loop over point on surface
%             disp(j);  
            warning off;
            surface_vec = [x_intercept_list(j); y_intercept_list(j); mean_surface_intercept_z];
%             disp(strcat('surface vec: ', num2str(surface_vec')));

            fermat_function_instance = @(glass_intercept_guess) fermat_function(glass_intercept_guess,...
                mean_glass_top_intercept_z, mean_glass_bottom_intercept_z, screen_vec, surface_vec,...
                mu_liquid_air, mu_glass_air);    
    
            glass_intercept = fminsearch(fermat_function_instance, fermat_guess, options_fermat_function);

            glass_bottom_vec = [glass_intercept(1:2); mean_glass_bottom_intercept_z];
            glass_top_vec = [glass_intercept(3:4); mean_glass_top_intercept_z];
            
            glass_bottom_intercept_inv_temp(:, j) = glass_bottom_vec';
            glass_top_intercept_inv_temp(:, j) = glass_top_vec';

%             disp(strcat('glass bottom vec:', num2str(glass_bottom_vec')));
%             disp(strcat('glass top vec:', num2str(glass_top_vec')));
            
            % verifying Snell's law
%             disp('air-glass interface check');
            incident_ray = screen_vec - glass_bottom_vec; % vector direction reversed to measure angle correct
            refracted_ray = glass_top_vec - glass_bottom_vec;
            theta_i_test = acosd(dot(incident_ray,[0; 0; -1])/norm(incident_ray));
% %             disp(theta_i_test);
            theta_r_test = acosd(dot(refracted_ray,[0; 0; 1])/norm(refracted_ray));
% %             disp(theta_r_test);
            snell_check_1(i, j) = sind(theta_i_test) / sind(theta_r_test);
%             disp('snell check air glass:'); disp(snell_check_1(i, j));
            
%             disp('glass-fluid interface check');
            incident_ray = glass_bottom_vec - glass_top_vec; % vector direction reversed to measure angle correct
            refracted_ray = surface_vec - glass_top_vec;
            theta_i_test = acosd(dot(incident_ray,[0; 0; -1])/norm(incident_ray));
%             disp(theta_i_test);
            theta_r_test = acosd(dot(refracted_ray,[0; 0; 1])/norm(refracted_ray));
%             disp(theta_r_test);
            snell_check_2(i, j) = sind(theta_i_test) / sind(theta_r_test);
%             disp('snell check glass liquid:'); disp(snell_check_2(i, j));

            %------------------------------------------

            unit_incident_vec = -1*(surface_vec-glass_top_vec)/norm(surface_vec-glass_top_vec);
            unit_refracted_vec = (source-surface_vec)/norm(source-surface_vec);

%             disp('incident vec:'); disp(-1*(surface_vec-glass_top_vec)');
%             disp('unit incident vec:'); disp(unit_incident_vec');
%             disp('refracted vec:'); disp((source-surface_vec)');
%             disp('unit refracted vec:'); disp(unit_refracted_vec');

            snell_function_instance = @(coeff_guess) snell_function(coeff_guess, unit_incident_vec, unit_refracted_vec, mu_liquid_air);

            coeff = fsolve(snell_function_instance, coeff_guess, options_snell_function);
%             disp('coeff'); disp(coeff');
            surface_normal = coeff(1)*unit_incident_vec + coeff(2)*unit_refracted_vec;
%             disp('surface normal that ensures snell's law'); disp(surface_normal');
            
            % surface normal obtained above is directed downwards, towards
            % incident ray, thus we are looking for opoosite of that
            unit_surface_normal = -1*surface_normal/norm(surface_normal);
%             disp('optimum unit surface normal:'); disp(unit_surface_normal');
            unit_surface_normal_calc_temp_list(:,j) = unit_surface_normal;

            % verifying Snell's law   
    %         disp('fluid-air check');
            theta_i_test = acosd(dot(unit_incident_vec, unit_surface_normal)); % vector direction reversed to measure angle correct
%             disp(strcat('snell test 3 theta_i_test: ', num2str(theta_i_test)));
            theta_r_test = acosd(dot(unit_refracted_vec, unit_surface_normal));
%             disp(strcat('snell test 3 theta_r_test: ', num2str(theta_r_test)));
            snell_check_3(i, j) = sind(theta_i_test) / sind(theta_r_test);
%             disp(strcat('snell test 3 check liquid air: ', num2str(snell_check_3(i, j))));
            
            % verify if the surface normal is a valid surface normal for an axi-symmetric case
            % i.e. it lies in the span of 2 vectors viz. vertical and radius vector
            center_vec = [0;0;fluid_depth]; % for current surface this is zero
            vertical_vec = [0;0;1]; % always vertical
            radius_vec = surface_vec - center_vec;
%             disp('radius vec:'); disp(radius_vec');
            vert_rad_cross = cross(vertical_vec, radius_vec);
            unit_vert_rad_cross = vert_rad_cross / norm(vert_rad_cross);

            % this can be considered like autocorrelation calculations
            surface_normal_coplanar_error(i, j)= abs(dot(unit_surface_normal, unit_vert_rad_cross));
%             disp(strcat('surface normal coplanar error: ', num2str(surface_normal_coplanar_error(i, j))));
    %         scatter3(screen_vec(1), screen_vec(2), screen_vec(3));
    %         scatter3(bottom_vec(1), bottom_vec(2), bottom_vec(3));
    %         scatter3(surface_vec(1), surface_vec(2), surface_vec(3));
            
            %plot source to surface vectors
            
%            disp('glass top vec'); disp(glass_top_vec');

%             plot3([screen_vec(1);glass_bottom_vec(1)], [screen_vec(2);glass_bottom_vec(2)], [screen_vec(3);glass_bottom_vec(3)]);
%             hold on;
%             plot3([glass_top_vec(1);glass_bottom_vec(1)], [glass_top_vec(2);glass_bottom_vec(2)], [glass_top_vec(3);glass_bottom_vec(3)]);
%             hold on;
%             plot3([glass_top_vec(1);surface_vec(1)], [glass_top_vec(2);surface_vec(2)], [glass_top_vec(3);surface_vec(3)], 'r');
%             hold on;
%             plot3([surface_vec(1);source(1)], [surface_vec(2);source(2)], [surface_vec(3);source(3)]);
%             hold on;
%             translated_unit_refracted_vec = 1*unit_refracted_vec + surface_vec;
%             plot3([surface_vec(1);translated_unit_refracted_vec(1)], [surface_vec(2);translated_unit_refracted_vec(2)], [surface_vec(3);translated_unit_refracted_vec(3)], 'm');
%             hold on;            
%             translated_unit_surface_normal = 1*unit_surface_normal + surface_vec;
            % surface normal list plot
%             plot3([surface_vec(1); translated_unit_surface_normal(1)], [surface_vec(2); translated_unit_surface_normal(2)], [surface_vec(3); translated_unit_surface_normal(3)], 'b');
%             hold on;
%             disp('----------')
        end  
        snell_error = [snell_check_1 snell_check_2 snell_check_3];       
        [~, id] = min(surface_normal_coplanar_error(i, :));            
        surface_intercept_calculated(1, i) = x_intercept_list(id);        
        surface_intercept_calculated(2, i) = y_intercept_list(id);
        unit_surface_normal_calc_final_list(:, i) = unit_surface_normal_calc_temp_list(:, id);
        glass_bottom_intercept_inv(:, i) = glass_bottom_intercept_inv_temp(:, id);
        glass_top_intercept_inv(:, i) = glass_top_intercept_inv_temp(:, id);

%         disp(strcat('surface normal coplanar error min ID: ', num2str(id)));
%         disp(strcat('surface normal coplanar error min VAL: ', num2str(val)));                
%         disp(strcat('final surface intercept:', num2str(surface_intercept_calculated(i)');
%         disp('corresponding final unit surface normal:'); disp(unit_surface_normal_calc_final_list(:, i)');
%         disp(strcat('final glass intercept bottom: ', num2str(glass_bottom_intercept_inv(:, i)')));
%         disp(strcat('final glass intercept top: ', num2str(glass_top_intercept_inv(:, i)')));  
    end
    screen_vec = [screen_intercept, ones(length(screen_intercept), 1)*mean_screen_intercept_z]';   

    %--------------------------------------------------------------------------
    % calculating slope values to be used to reconstruct the surface using the surface normals calculated
    
    dF_dy_list_non_smooth = 1e6*ones(1, length(surface_intercept_calculated(2, :))-1);
    
    % numerical integration method
    for i=1:length(surface_intercept_calculated(2, :))
%         this is based on unit surface normal calculation of F(x, y, z) eq
%         5 in paper 1
        surface_normal_scale = 1 / unit_surface_normal_calc_final_list(3, i);
        dF_dy_list_non_smooth(i) = -1*unit_surface_normal_calc_final_list(2, i)*surface_normal_scale; 
    end
    %-----------------------------------------
    % smoothen the gradient in y-direction               
    [xData, yData] = prepareCurveData(surface_intercept_calculated(2, :), dF_dy_list_non_smooth);

    % Set up fittype and options
    ft = fittype('smoothingspline');
    opts = fitoptions('Method', 'SmoothingSpline');
    opts.SmoothingParam = smoothening_param;
    [fitresult, ~] = fit(xData, yData, ft, opts); 
    
    dF_dy_list_smooth = fitresult(surface_intercept_calculated(2, :));        
    %-----------------------------------------
    % numerical integration using trapezoidal method for numerical integration instead of
    % simpsons as shown above because of alternate 0s and issues in spline
    % interpolation for intermediate points
    surface_intercept_calculated(3, :) = cumtrapz(surface_intercept_calculated(2, :), dF_dy_list_smooth);   
    %     disp(surface_intercept_calculated);
end

%--------------------------------------------------------------------------

 function time = fermat_function(glass_intercept_guess, glass_top_intercept_z, glass_bottom_intercept_z, screen_vec, surface_vec,  mu_liquid_air, mu_glass_air)
    guess_vec_bottom = [glass_intercept_guess(1:2); glass_bottom_intercept_z];
    guess_vec_top = [glass_intercept_guess(3:4); glass_top_intercept_z];
    dist1 = norm(screen_vec - guess_vec_bottom); % all in air
    dist2 = norm(guess_vec_bottom - guess_vec_top); % all in glass
    dist3 = norm(guess_vec_top - surface_vec); % all in liquid
%     disp(glass_intercept_guess);
%     disp(dist1);
%     disp(dist2);
    % using relative velocities of light in different media to calculate scaled time
%     time = (dist1/mu_liquid_air)+(dist2*mu_liquid_glass)+dist3; % using velocity in liquid as reference
    time = dist1*1 + dist2*mu_glass_air + dist3*mu_liquid_air; % using velocity in air as reference    
%     disp(time);
%     disp('----');
 end
 
% %-------------------------------------------------------------------------- 

 function error = snell_function(coeff_guess, unit_incident_vec, unit_refracted_vec, mu_liquid_air)
%     disp('coeff guess'); disp(coeff_guess');
%     disp('unit incident vec'); disp(unit_incident_vec');
%     disp('unit refracted vec'); disp(unit_refracted_vec');
    surface_normal = coeff_guess(1)*unit_incident_vec + coeff_guess(2)*unit_refracted_vec;
%     disp('calculated surface normal'); disp(surface_normal');
%     disp('norm of surface normal vec'), disp(norm(surface_normal));
    unit_surface_normal = surface_normal / norm(surface_normal);
%     disp('calculated unit surface normal'); disp(unit_surface_normal');

%     theta_i = acos(dot(unit_surface_normal, unit_incident_vec));
%     theta_r = acos(dot(unit_surface_normal, unit_refracted_vec));
%     theta_i = acos(dot(-1*unit_surface_normal, unit_incident_vec));
%     theta_r = acos(dot(unit_surface_normal, unit_refracted_vec));

    theta_i = acosd(dot(unit_surface_normal, unit_incident_vec));
    theta_r = acosd(dot(unit_surface_normal, unit_refracted_vec));

    if theta_i>90
        theta_i=180-theta_i;
    end
    
    if theta_r>90
        theta_r=180-theta_r;
    end       
    error = mu_liquid_air * sind(theta_i) - sind(theta_r);
    
%     error = 1*abs(error);  % might be useful to scale the error to avoid
%     very tiny numerical value 
%     disp(error);
%     disp('---snell function iter---')
 end
 
 %--------------------------------------------------------------------------
 function [screen_intercept] = processing(img_path, crop_box, std_dev_gaussian_filter,...
     threshold_filter, camera_params, overwrite_flag, file_string_length)
    
    img_name = strcat(img_path(file_string_length:end-4));
%     disp(img_name); %disp(isfile(imageFile));
    
    % image processing and calculating mean curve
    imageFile = strcat(img_name, '_processed.png');
%     imageFile_edge = strcat(img_name, '_edge.png');
%     imageFile_mean = strcat(img_name, '_mean.png');

    if ~isfile(imageFile) || overwrite_flag
%         disp('either file not found or overwrite flag is true');           
        img = imread(img_path); % reading image
        % undistorting the expt image using the calibration data to remove effects of lens distortion
        [img, ~] = undistortImage(img, camera_params);
%         figure; imshow(img);  
    %                 w = waitforbuttonpress;

        img_crop = img;
        % area outside the crop box is darkened instead of cropping
        % since the dimension of image in camera calibration and
        % further must remain same
        [mrows, ncols] = size(img); % for matlab, x axis is col and y axis is row, maybe
    %                 disp(ncols);
    %                 disp(mrows);
        for x = 1:mrows
            for y = 1:ncols
                if ~((crop_box(1)<x) && (x<crop_box(2)) && (crop_box(3)<y) && (y<crop_box(4)))
                    img_crop(x, y) = 0;
                end                           
            end
        end
        %img_crop = img(crop_box(1):crop_box(2), crop_box(3):crop_box(4)); % cropping
        img_crop_gauss = imgaussfilt(img_crop, std_dev_gaussian_filter); % gaussian filter
        img_crop_thres = 255*(img_crop_gauss>threshold_filter); % threshold filter
%         figure; imshow(img_crop_thres);

        % saving processed images for further analysis
        %set(gcf, 'InvertHardCopy', 'off');
        %set(gcf,'color','black');
    %             print('./inverse_solution_plots/mean_curve/actual_img.png', '-dpng');
        % just filename can be extracted from string by using strfind(waveImageList{1}, '\')
        % ~ to invert the image to have white background
%         disp(cd); fprintf('\n\n');
%         imwrite(~img_crop_thres, imageFile);
    %                 saveas(f, strcat(waveImageList{i}(104:end-5), '_processed.png')); % save as .png file
    %                 imshow(img_crop_thres);
    %                 print(gcf, strcat(imageFile(1:end-4), '.eps'),'-depsc2','-r300');
    %                 close;
                       
        % need only the coordinates for image to world mapping
        [mrows, ncols] = find(img_crop_thres==255);        
        screen_intercept = [ncols mrows]; % check the image axes
        % ncols show x coordinate and mcols show Y cooridnate accoring
        % tho the axes of calibraiton pattern
    else
        disp('File already processed. Delete processed files for re-processing.')
        screen_intercept = [];
    end
 end

  %--------------------------------------------------------------------------
function [mapped_translated_screen_intercept_downsized, fitresult_yx] = coordinate_transform(camera_params,...
    rotation_mat, translation_vec, screen_intercept, translate_x, translate_y, num_intercept, zero_freq_flag,...
    smoothening_para)
  
    mapped_screen_intercept = pointsToWorld(camera_params, rotation_mat, translation_vec, screen_intercept);           

    %converting coordinate system to mtch with the inverse solver
    %system, this requires just X, Y to Y, X transformation
    mapped_screen_intercept = [mapped_screen_intercept(:, 2) mapped_screen_intercept(:, 1)];

    % process image of origin for translation un-deform images before analysis
    % then obtain world coordinates of wave surface and then translate it to to
    % actual desired location by translation of origin, hope rotation is not
    % needed as we are maintaining alignments

    % these are center coordinates in calibration pattern origin
    mapped_translated_screen_intercept = mapped_screen_intercept - [translate_x, translate_y];

    % Set up fittype and options
    if zero_freq_flag 
        ft = fittype('poly1');
        opts = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'LAR');
%         disp('line');
    else
        ft = fittype('smoothingspline');
        opts = fitoptions('Method', 'SmoothingSpline'); 
        opts.SmoothingParam = smoothening_para;  
%         disp('spline');
    end
    
    % Fit model to data, note swapped x and y locations
    [mapped_translated_screen_intercept(:,2), mapped_translated_screen_intercept(:,1)] = prepareCurveData(mapped_translated_screen_intercept(:,2), mapped_translated_screen_intercept(:,1));
    [fitresult_yx, ~] = fit(mapped_translated_screen_intercept(:,2), mapped_translated_screen_intercept(:,1), ft, opts);   
    
    % creating fine data to get more intermediate points
    dy_fine = (max(mapped_translated_screen_intercept(:, 2))-min(mapped_translated_screen_intercept(:, 2)))/num_intercept;
    mapped_translated_screen_intercept_downsized(:, 2) = mapped_translated_screen_intercept(1, 2):dy_fine:mapped_translated_screen_intercept(end, 2);
    mapped_translated_screen_intercept_downsized(:, 1) = fitresult_yx(mapped_translated_screen_intercept_downsized(:, 2));
    
    % plot fit with data
    figure;
    scatter(mapped_translated_screen_intercept(:,1), mapped_translated_screen_intercept(:,2), 'r.');
    hold on;  
    plot(mapped_translated_screen_intercept_downsized(:,1), mapped_translated_screen_intercept_downsized(:,2), 'g.');
    legend('mapped to physical world', 'downsized data screen intercept');
    ylim([min(mapped_translated_screen_intercept_downsized(:,2))-20 max(mapped_translated_screen_intercept_downsized(:,2))+60]);
    if zero_freq_flag
        title('Mapped Translated Data - Zero Freq');
    else
        title('Mapped Translated Data');
    end        
    xlabel('x-axis (mm)'); ylabel('y-axis (mm)');
    grid on; grid minor;
end
 
 