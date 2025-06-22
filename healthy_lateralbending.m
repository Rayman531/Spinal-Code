%% Load Optotrak Data (Motion Tracking)
optotrak_data = dlmread('L3_B_Healthy_LateralBending_3d.csv', ',', 5, 0);
frame_number = optotrak_data(:, 1);
sampling_rate = 100;  % Hz
duration = 30;  % seconds
num_samples = sampling_rate * duration; 

% Extract marker positions (last 27 columns)
marker_data = optotrak_data(:, end-27:end);
left = marker_data(:, 1:3); 
right = marker_data(:, 10:12); 
anterior = marker_data(:, 19:21);

% Compute coordinate system (Bottom Vertebra)
y = left - right; 
v = anterior - right; 
z = cross(v, y); 
x = cross(y, z); 
T_bottom = (left + right) / 2;  % Translation vector

% Normalize axes
i = x ./ vecnorm(x, 2, 2);
j = y ./ vecnorm(y, 2, 2);
k = z ./ vecnorm(z, 2, 2);

% Compute Rotation Matrices for Bottom Vertebra
R_bottom = zeros(3, 3, size(optotrak_data, 1));
for t = 1:size(optotrak_data, 1)
    R_bottom(:,:,t) = [i(t,:); j(t,:); k(t,:)];
end

%% Compute Augmented Transformation Matrix for Bottom Vertebra
T_bottom_hom = zeros(4, 4, size(optotrak_data, 1));
for t = 1:size(optotrak_data, 1)
    T_bottom_hom(:,:,t) = [R_bottom(:,:,t), T_bottom(t,:)'; 0, 0, 0, 1];
end

%% Extract Top Vertebra Markers
left_top = marker_data(:, 4:6);  
right_top = marker_data(:, 13:15);  
anterior_top = marker_data(:, 22:24);  

%% Convert Top Vertebra Landmarks into Bottom Coordinates
left_top_transformed = zeros(size(left_top));
right_top_transformed = zeros(size(right_top));
anterior_top_transformed = zeros(size(anterior_top));

for t = 1:size(optotrak_data, 1)
    R_b = R_bottom(:,:,t);
    T_b = T_bottom(t,:)';
    left_top_transformed(t,:) = R_b' * (left_top(t,:)' - T_b);
    right_top_transformed(t,:) = R_b' * (right_top(t,:)' - T_b);
    anterior_top_transformed(t,:) = R_b' * (anterior_top(t,:)' - T_b);
end

%% Compute Coordinate System for Top Vertebra in Bottom Frame
y_top = left_top_transformed - right_top_transformed;
v_top = anterior_top_transformed - right_top_transformed;
z_top = cross(v_top, y_top);
x_top = cross(y_top, z_top);

% Normalize axes
i_top = x_top ./ vecnorm(x_top, 2, 2);
j_top = y_top ./ vecnorm(y_top, 2, 2);
k_top = z_top ./ vecnorm(z_top, 2, 2);

% Compute Rotation Matrices for Top Vertebra in Bottom Frame
R_top_f = zeros(3, 3, size(optotrak_data, 1));
for t = 1:size(optotrak_data, 1)
    R_top_f(:,:,t) = [i_top(t,:); j_top(t,:); k_top(t,:)];
end

%% Compute Translation Vector for Top Vertebra in Bottom Frame
T_top_f = (left_top_transformed + right_top_transformed) / 2;

%% Compute Augmented Transformation Matrix for Top Relative to Bottom
M_top_f = zeros(4, 4, size(optotrak_data, 1));
for t = 1:size(optotrak_data, 1)
    M_top_f(:,:,t) = [R_top_f(:,:,t), T_top_f(t,:)'; 0, 0, 0, 1];
end

%% Compute M (Transformation from Frame 1 to Frame f)
M = zeros(4, 4, size(optotrak_data, 1));
M(:,:,1) = eye(4);  % Identity matrix for frame 1
for t = 2:size(optotrak_data, 1)
    M(:,:,t) = inv(M_top_f(:,:,1)) * M_top_f(:,:,t);
end


%% Extract Euler Angles (α, β, γ)
beta = zeros(size(optotrak_data, 1), 1);

for t = 1:size(optotrak_data, 1)
    R_xyz = M(1:3, 1:3, t);
    
    % Extract Euler angles from R_xyz (XYZ convention)
    beta(t) = asin(-R_xyz(3,1));  
end


alpha = zeros(size(optotrak_data, 1), 1);
for t = 1:size(optotrak_data, 1)
    alpha(t) = asin((R_xyz(3,2))/cos(beta(t)));   % α (Axial Rotation)
end

% Convert β to degrees
alpha_deg = rad2deg(alpha);


%% Load Spine Flexor Data (Force Measurement)
lvm_data = dlmread('L3_B_Healthy_LateralBending.lvm');
time_lvm = lvm_data(:, 8);  
My_moment = lvm_data(:, 6);

% Interpolate My to match Optotrak sampling rate
time_optotrak = linspace(time_lvm(1), time_lvm(end), num_samples);
My_interp = interp1(time_lvm, My_moment, time_optotrak, 'linear');

%% Compute ROM and Neutral Zone for α, β, γ
ROM = max(alpha_deg) - min(alpha_deg);
NZ = max(alpha_deg(abs(My_interp) < 0.1));

%% Plot Hysteresis Curve for β
figure;
hold on;
plot(My_interp(1,1500:2000), alpha_deg(1500:2000,1), 'r', 'LineWidth', 1.5);
xlabel('Moment [Nm]');
ylabel('\alpha Angle [°]');
title('Spinal LateralBending');

yline(NZ, '--k', 'NZ', 'LabelHorizontalAlignment', 'left');
yline(max(alpha_deg), '-.b', '+ROM', 'LabelHorizontalAlignment', 'left');
yline(min(alpha_deg), '-.b', '-ROM', 'LabelHorizontalAlignment', 'left');

grid on;
hold off;


%% Bar Plot of NZ and ROM for α, β, γ
figure;
bar([NZ, ROM]);
xticklabels({'Lateral Bending (alpha)'});
ylabel('Angle (°)');
title('Comparison of NZ and ROM for alpha');
grid on;
