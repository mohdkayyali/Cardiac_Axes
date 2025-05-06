
clearvars
clear all

vs =readtable('/home/mka23/Desktop/Anatomical_electrical_axes/Anatomical_electrical_analysis/vector_data_healthy.csv');
Electrical_vectors_byMethod = vs(:,["EM1","EM2","EM3","EM4","EM5"]);
Anatomical_vectors_byMethod = vs(:,["AM1","AM2","AM3","AM4","AM5"]);
EVs = Electrical_vectors_byMethod; 
AVs = Anatomical_vectors_byMethod;

%% 
A_vs_E_cosine = [];
A_vs_E_angles = [];
angle_from_mean_E_allM = [];

for EM = 1:5
        strVectors = Electrical_vectors_byMethod{:, EM};
        numVectors = cellfun(@(x) str2num(x), strVectors, 'UniformOutput', false);
        electrical_axis_M = cell2mat(numVectors)';
        method_mean_vector = mean(electrical_axis_M,2);
        for p = 1:size(Electrical_vectors_byMethod,1)
            electrical_curr_vector = electrical_axis_M(:,p);
            angle_from_mean_E_allM(p,EM) = rad2deg(acos(dot(method_mean_vector, electrical_curr_vector) / (norm(method_mean_vector) * norm(electrical_curr_vector))));
        end
    for AM = 1:5
        % anatomical_axis_M = Anatomical_vectors_byMethod.(['A_M',num2str(AM)]);
        strVectors = Anatomical_vectors_byMethod{:, AM};
        numVectors = cellfun(@(x) str2num(x), strVectors, 'UniformOutput', false);
        anatomical_axis_M = cell2mat(numVectors)';
        for p = 1:size(Electrical_vectors_byMethod,1)
            anatomical_curr_vector = -1.*anatomical_axis_M(:,p);
            electrical_curr_vector = electrical_axis_M(:,p);
            A_vs_E_cosine(p,(EM - 1) * 5 + AM) = dot(anatomical_curr_vector, electrical_curr_vector) / (norm(anatomical_curr_vector) * norm(electrical_curr_vector)); 
            A_vs_E_angles(p,(EM - 1) * 5 + AM) = rad2deg(acos(dot(anatomical_curr_vector, electrical_curr_vector) / (norm(anatomical_curr_vector) * norm(electrical_curr_vector))));
        end
    end
end


angle_filter = abs(angle_from_mean_E_allM - mean(angle_from_mean_E_allM)) <= 3* std(angle_from_mean_E_allM);
sum(angle_filter)

%% 

geodesic_dist_newE_vs_oldE = cell(5);
geodesic_dist_newE_vs_anat = cell(5);
geodesic_dist_E_vs_anat = cell(5);
for EM = 1:5
    for AM = 1:5
        % Extract and normalize E_vectors
        strVectors = EVs{:, EM};
        numVectors = cellfun(@(x) str2num(x), strVectors, 'UniformOutput', false);
        E_vectors = cell2mat(numVectors)';
        strVectors = AVs{:, AM};
        numVectors = cellfun(@(x) str2num(x), strVectors, 'UniformOutput', false);
        A_vectors = -1*cell2mat(numVectors)';
        E_vectors = E_vectors ./ vecnorm(E_vectors, 2, 1);
        A_vectors = A_vectors ./ vecnorm(A_vectors, 2, 1);
        E_vectors = E_vectors(:,logical(angle_filter(:, EM)));
        A_vectors = A_vectors(:,logical(angle_filter(:, EM)));



        % Number of samples
        numSamples = size(E_vectors, 2);
        
        % Define split ratio
        trainRatio = 0.8; % 80% training data
        splitIndex = floor(trainRatio * numSamples);
        
        % Randomly shuffle indices
        rng(1); % For reproducibility
        indices = randperm(numSamples);
        
        % Split data into training and testing sets
        trainIndices = indices(1:splitIndex);
        testIndices = indices(splitIndex+1:end);
        
        % Training data
        E_train = E_vectors(:, trainIndices);
        A_train = A_vectors(:, trainIndices);
        
        % Testing data
        E_test = E_vectors(:, testIndices);
        A_test = A_vectors(:, testIndices);



        % Compute optimal rotation matrix using SVD
        H = E_train * A_train';
        [U, ~, V] = svd(H);
        T_rotation = U * V';  % This is guaranteed to be a pure rotation matrix

        % Transform A_vectors to obtain E_new
        E_new_rotation = T_rotation * A_test;
        % Normalize the transformed vectors (should be very close to unit vectors already, but just to be safe)
        E_new_rotation = E_new_rotation ./ vecnorm(E_new_rotation, 2, 1);

        % Compute geodesic distance between corresponding vectors in E_vectors and E_new_rotation
        geo_dist_error_rotation = zeros(size(E_test, 2), 1);
        for i = 1:size(E_test, 2)
            geo_dist_error_rotation(i) = acos(dot(E_test(:,i), E_new_rotation(:,i)));
        end
        geodesic_dist_newE_vs_oldE{EM,AM} = geo_dist_error_rotation;

        % Compute geodesic distance between corresponding vectors in E_new_rotation and A_vectors
        geo_dist_error_anatelectnew = zeros(size(E_test, 2), 1);
        for i = 1:size(E_test, 2)
            geo_dist_error_anatelectnew(i) = acos(dot(A_test(:,i), E_new_rotation(:,i)));
        end
        geodesic_dist_newE_vs_anat{EM,AM} = geo_dist_error_anatelectnew;

        % Compute geodesic distance between corresponding vectors in A_vectors and E_vectors
        geo_dist_error_anatelect = zeros(size(E_test, 2), 1);
        for i = 1:size(E_test, 2)
            geo_dist_error_anatelect(i) = acos(dot(A_test(:,i), E_test(:,i)));
        end
        geodesic_dist_E_vs_anat{EM,AM} = geo_dist_error_anatelect;
    end
end


%% 
% Prepare data matrices for heatmaps
methods_to_use = [1 2 3 4]; % Indices for maxQRS, meanQRS, v-avgQRS, eig1QRS
n_methods = length(methods_to_use);

% Initialize matrices for all comparisons
MGD_E_vs_A = zeros(n_methods, 5);
STD_E_vs_A = zeros(n_methods, 5);
MGD_newE_vs_E = zeros(n_methods, 5);
STD_newE_vs_E = zeros(n_methods, 5);
MGD_newE_vs_A = zeros(n_methods, 5);
STD_newE_vs_A = zeros(n_methods, 5);

% Fill matrices
for i = 1:n_methods
    for j = 1:5
        EM = methods_to_use(i);
        
        % Original E vs A distances
        orig_distances = geodesic_dist_E_vs_anat{EM,j};
        MGD_E_vs_A(i,j) = mean(orig_distances);
        STD_E_vs_A(i,j) = std(orig_distances);
        
        % Transformed E vs true E distances
        trans_distances = geodesic_dist_newE_vs_oldE{EM,j};
        MGD_newE_vs_E(i,j) = mean(trans_distances);
        STD_newE_vs_E(i,j) = std(trans_distances);
        
        % Transformed E vs A distances
        transA_distances = geodesic_dist_newE_vs_anat{EM,j};
        MGD_newE_vs_A(i,j) = mean(transA_distances);
        STD_newE_vs_A(i,j) = std(transA_distances);
    end
end



%% 
figure();


heatmap(MGD_newE_vs_E, 'ColorLimits', [0.3, 0.4], 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 13, 'GridVisible', 'off');
title('MGD Between Predicted and True Electrical Vectors');
ylabel('Electrical Axis Methods');
xlabel('Anatomical Axis Methods');
colorbar;
colormap('gray');
ax = gca;
ax.XDisplayLabels = {'$PC1_{LV}$', '$PC1_{LRV}$', 'MVA', 'MAVA', 'VPA'}; 
ax.YDisplayLabels = {'maxQRS', 'meanQRS','v-avgQRS', 'eig1QRS'}; 



%% 

% Create figure with subplots
figure('Position', [100 100 1800 800]);

% Labels for axes
EM_labels = {'maxQRS', 'meanQRS', 'v-avgQRS', 'eig1QRS'};
AM_labels = {'PC1_{LV}', 'PC1_{LRV}', 'MVA', 'MAVA', 'VPA'};

% Plot 1: Original E vs A MGD
subplot(2,2,1);
heatmap(MGD_E_vs_A);
title('MGD Between Electrical Vectors and Anatomical Vectors');
ylabel('Electrical Axis Methods');
xlabel('Anatomical Axis Methods');
colorbar;
colormap('gray');
ax = gca;
ax.XDisplayLabels = {'PC1_{LV}', 'PC1_{LRV}', 'MVA', 'MAVA', 'VPA'}; 
ax.YDisplayLabels = {'maxQRS', 'meanQRS','v-avgQRS', 'eig1QRS'}; 

% Plot 2: Original E vs A STD
subplot(2,2,3);
heatmap(STD_E_vs_A);
title('Standard Deviation of GD');
ylabel('Electrical Axis Methods');
xlabel('Anatomical Axis Methods');
colorbar;
colormap('gray');
ax = gca;
ax.XDisplayLabels = {'PC1_{LV}', 'PC1_{LRV}', 'MVA', 'MAVA', 'VPA'}; 
ax.YDisplayLabels = {'maxQRS', 'meanQRS','v-avgQRS', 'eig1QRS'}; 

% Plot 3: Transformed E vs true E MGD
subplot(2,2,2);
heatmap(MGD_newE_vs_E);
title('MGD Between Predicted and True Electrical Vectors');
ylabel('Electrical Axis Methods');
xlabel('Anatomical Axis Methods');
colorbar;
colormap('gray');
ax = gca;
ax.XDisplayLabels = {'PC1_{LV}', 'PC1_{LRV}', 'MVA', 'MAVA', 'VPA'}; 
ax.YDisplayLabels = {'maxQRS', 'meanQRS','v-avgQRS', 'eig1QRS'}; 


% Plot 4: Transformed E vs true E STD
subplot(2,2,4);
heatmap(STD_newE_vs_E);
title('Standard Deviation of GD');
ylabel('Electrical Axis Methods');
xlabel('Anatomical Axis Methods');
colorbar;
colormap('gray');
ax = gca;
ax.XDisplayLabels = {'PC1_{LV}', 'PC1_{LRV}', 'MVA', 'MAVA', 'VPA'}; 
ax.YDisplayLabels = {'maxQRS', 'meanQRS','v-avgQRS', 'eig1QRS'}; 



