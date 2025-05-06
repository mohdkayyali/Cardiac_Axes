clearvars
clear all

vs =readtable('/home/mka23/Desktop/Anatomical_electrical_axes/Anatomical_electrical_analysis/vector_data_healthy.csv');
Electrical_vectors_byMethod = vs(:,["EM1","EM2","EM3","EM4","EM5"]);
Anatomical_vectors_byMethod = vs(:,["AM1","AM2","AM3","AM4","AM5"]);

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
computeAngle = @(vector) atan2d(vector(2), vector(1));
EVs = Electrical_vectors_byMethod; 
AVs = Anatomical_vectors_byMethod;
alpha_E = zeros(height(EVs), 5); 
beta_E = zeros(height(EVs), 5); 
gamma_E = zeros(height(EVs), 5); 
alpha_A = zeros(height(EVs), 5); 
beta_A = zeros(height(EVs), 5); 
gamma_A = zeros(height(EVs), 5); 
for M = 1:5
    strVectors = EVs{:, M};
    numVectors = cellfun(@(x) str2num(x), strVectors, 'UniformOutput', false);
    E_vectors = cell2mat(numVectors)';

    strVectors = AVs{:, M};
    numVectors = cellfun(@(x) str2num(x), strVectors, 'UniformOutput', false);
    A_vectors = -1.*cell2mat(numVectors)';

    % A_vectors = table2array(AVs(:,M));
    % A_vectors = -1 * cell2mat(A_vectors');

    % Calculate angles between vectors and respective horizontal axes for each plane
    for i = 1:length(A_vectors)
        for j = 1:3
            switch j
                case 1 % ZX
                    projected_vector_E = [E_vectors(1,i); E_vectors(3,i)];
                    alpha_E(i,M) = computeAngle(projected_vector_E);
                    projected_vector_A = [A_vectors(1,i); A_vectors(3,i)];
                    alpha_A(i,M) = computeAngle(projected_vector_A);
                case 2 % ZY
                    projected_vector_E = [E_vectors(2,i); E_vectors(3,i)];
                    beta_E(i,M) = computeAngle(projected_vector_E);
                    projected_vector_A = [A_vectors(2,i); A_vectors(3,i)];
                    beta_A(i,M) = computeAngle(projected_vector_A);
                case 3 % YX
                    projected_vector_E = [E_vectors(1,i); E_vectors(2,i)];
                    gamma_E(i,M) = computeAngle(projected_vector_E);
                    projected_vector_A = [A_vectors(1,i); A_vectors(2,i)];
                    gamma_A(i,M) = computeAngle(projected_vector_A);
                end        
        end
    end
end

delta_alpha = zeros(size(alpha_E, 1), 25);
delta_beta = zeros(size(beta_E, 1), 25);
delta_gamma = zeros(size(gamma_E, 1), 25);

% Calculate delta values for alpha, beta, and gamma
for EM = 1:5
    for AM = 1:5
        delta_alpha(:, (EM-1)*5 + AM) =  alpha_A(:, AM) - alpha_E(:, EM);

        delta_beta(:, (EM-1)*5 + AM) =  beta_A(:, AM) - beta_E(:, EM);

        delta_gamma(:, (EM-1)*5 + AM) =  gamma_A(:, AM) - gamma_E(:, EM) ;

        % delta_alpha(:, (EM-1)*5 + AM) =  alpha_E(:, EM) - alpha_A(:, AM);
        % 
        % delta_beta(:, (EM-1)*5 + AM) =  beta_E(:, EM) - beta_A(:, AM);
        % 
        % delta_gamma(:, (EM-1)*5 + AM) =  gamma_E(:, EM) - gamma_A(:, AM);
    end
end

delta_alpha(delta_alpha < -180) = 360 + delta_alpha(delta_alpha < -180);
delta_beta(delta_beta < -180) = 360 + delta_beta(delta_beta < -180);
delta_gamma(delta_gamma < -180) = 360 + delta_gamma(delta_gamma < -180);

delta_alpha(delta_alpha > 180) = delta_alpha(delta_alpha > 180)-360;
delta_beta(delta_beta > 180) = delta_beta(delta_beta > 180)-360;
delta_gamma(delta_gamma > 180) = delta_gamma(delta_gamma > 180)-360;

angle_sets = {delta_alpha, delta_beta, delta_gamma};
data_e = {alpha_E,beta_E,gamma_E};
data_a = {alpha_A,beta_A,gamma_A};
%% SC figure planar
figure;

for w = 1:3
    subplot(2, 2, w+1);   
    hold on;
    min_std_index = 0;
    min_std = min(std(angle_sets{w})); % Initialize minimum std value
    stds = zeros(25,1);
    for i = 1:25
        data_m = angle_sets{w}(logical(angle_filter(:, ceil(i/5))), i);
        stds(i,1) = std(data_m);
        pd = fitdist(data_m, 'Normal');
        x_values = linspace(min(data_m), max(data_m), 150);
        y = pdf(pd, x_values);
        % disp(pd)
        plot(x_values, y, 'Color', [.7 .7 .7], 'LineWidth', 1.5,'HandleVisibility','off');
    end
    plot(x_values, y, 'Color', [.7 .7 .7], 'LineWidth', 1.5);

    % [~,o] = min(stds);
    o=5;
    data_m = angle_sets{w}(logical(angle_filter(:, ceil(o/5))), o);
    pd = fitdist(data_m, 'Normal');
    x_values = linspace(min(data_m), max(data_m), 150);
    y = pdf(pd, x_values);
    plot(x_values, y, 'black', 'LineWidth', 1.5);
    % o=25;
    % data_m = angle_sets{w}(logical(angle_filter(:, ceil(o/5))), o);
    % pd = fitdist(data_m, 'Normal');
    % x_values = linspace(min(data_m), max(data_m), 150);
    % y = pdf(pd, x_values);
    % plot(x_values, y, 'red', 'LineWidth', 1.5);

    disp(pd)
    disp(o)
    grid on;ylabel('Probability Density','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    if w==1
        title('Frontal Plane','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14)
    elseif w==2
        title('Sagittal Plane','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14)
    else
        title('Transverse Plane','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14)
    end
    xlim([-200 200]);
    xticks([-180:60:180])
    ylim([0 0.035])
    xlabel('Angle between anatomical-electrical pair','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

end
% legend('Other Methods','Anatomical Axis: VPA, Electrical Axis: maxQRS',...
%     'Anatomical Axis: VPA, Electrical Axis: eig1QRS', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12)


% title(['Relationship between A-M',num2str(A_M), '& E-M',num2str(E_M)]);
%% SC figure 3d

figure;
stdevs = [];
% subplot(2,2,1);
hold on;
title('Absolute 3D Angle','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14)
for i = [1:5,11:25,1,5]
    % i=20
    % data_m = A_vs_E_cosine(logical(angle_filter(:,ceil(i/5))),i);
    data_m = A_vs_E_angles(logical(angle_filter(:,ceil(i/5))),i);

    pd = fitdist(data_m, 'Normal');
    x_values = linspace(min(data_m), max(data_m), 500);
    y = pdf(pd, x_values);
    disp(pd)
    % plot(x_values, y, 'red','LineWidth', 0.5);

    stdevs(i) = pd.sigma;
    if i == 5
        plot(x_values, y, 'black','LineWidth', 2);
    % elseif i==25 
    %     plot(x_values, y, 'red','LineWidth', 2);
    elseif i ==2
        plot(x_values, y, 'Color', [.7 .7 .7],'LineWidth', 1.5);
    else
        plot(x_values, y, 'Color', [.7 .7 .7],'LineWidth', 1.5,'HandleVisibility','off');
    end

    
end
hold off;grid on;
xlabel('Angle between anatomical-electrical pair','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Other Methods','Anatomical Axis: VPA, Electrical Axis: maxQRS',...
     'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14)

% Add labels and title if needed
% xlabel('cosine(angle between anatomical-electrical pair)','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

ylabel('Probability Density','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

%% std heatmap
stdevs = stdevs([1:5,11:25])';
figure();
P = reshape(round(stdevs,2), [5,4]) %/max(mse_geo_lse)

xx = {'$PC1_{LV}$', '$PC1_{LRV}$', 'MVA', 'MAVA', 'VPA'};
% yy = {'$v_{maxQRS}$', '$v_{maxXYZ}$', '$v_{time-avg}$', '$v_{v-avg}$', '$v_{SVD}$'};
yy = {'maxQRS', 'meanQRS', 'v-avgQRS', 'eig1QRS'};
L = heatmap(xx, yy, P', 'ColorLimits', [16, 30], 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 13, 'GridVisible', 'off');

% L = heatmap(xx, yy, P', 'ColorLimits', [0.8, 2], 'FontSize', 14, 'GridVisible', 'off');
colormap(gray)
xlabel('Anatomical Axis Methods')
ylabel('Electrical Axis Methods');
title('Standard Deviation of angular differences in 3D')

%% 2D Gaussian

AM = 5; 
EM = 1;
EM_AM= (EM-1)*5 + AM;
filter = logical(angle_filter(:, ceil(EM_AM/5)));
% filter = Group_2;

% Get the 2-d samples for the "floor"
electrical_F_T = [alpha_E(filter,EM) , gamma_E(filter,EM)];
anatomical_F_T = [alpha_A(filter,AM) , gamma_A(filter,AM)];
delta_F_T = [delta_alpha(filter,EM_AM), delta_gamma(filter,EM_AM)];
figure; hold on; 

% Define grid for electrical data
x_electrical = linspace(min(electrical_F_T(:,1)), max(electrical_F_T(:,1)), 100);
y_electrical = linspace(min(electrical_F_T(:,2)), max(electrical_F_T(:,2)), 100);
[X_electrical, Y_electrical] = meshgrid(x_electrical, y_electrical);
% Compute kernel density estimation (KDE) for electrical data
Z_electrical = ksdensity(electrical_F_T, [X_electrical(:), Y_electrical(:)]);

% Define grid for anatomical data
x_anatomical = linspace(min(anatomical_F_T(:,1)), max(anatomical_F_T(:,1)), 100);
y_anatomical = linspace(min(anatomical_F_T(:,2)), max(anatomical_F_T(:,2)), 100);
[X_anatomical, Y_anatomical] = meshgrid(x_anatomical, y_anatomical);
% Compute kernel density estimation (KDE) for anatomical data
Z_anatomical = ksdensity(anatomical_F_T, [X_anatomical(:), Y_anatomical(:)]);

% Define grid for delta data
x_delta = linspace(min(delta_F_T(:,1)), max(delta_F_T(:,1)), 100);
y_delta = linspace(min(delta_F_T(:,2)), max(delta_F_T(:,2)), 100);
[X_delta, Y_delta] = meshgrid(x_delta, y_delta);
% Compute kernel density estimation (KDE) for anatomical data
Z_delta = ksdensity(delta_F_T, [X_delta(:), Y_delta(:)]);

% Plot contours of density distributions
contour(X_electrical, Y_electrical, reshape(Z_electrical, size(X_electrical)),6, 'LineColor', 'k','LineWidth',1.5);
hold on;
contour(X_anatomical, Y_anatomical, reshape(Z_anatomical, size(X_anatomical)),6, 'LineColor', 'r','LineWidth',1.5);
contour(X_delta, Y_delta, reshape(Z_delta, size(X_delta)),6,'LineColor','b', 'LineWidth',1.5)
% Make the figure look nice
grid on; 

xlabel('Frontal Plane','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12)
ylabel('Trasnverse Plane','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12)
xticks(-90:30:90)
xlim([-90 90])
yticks(-90:30:90)
ylim([-90 90])
zlim([-.001 0.05])
zlabel('Probabilty Density','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12)

box on
electrical = data_e{1}(filter, EM);
pd = fitdist(electrical, 'Normal');
x_values = linspace(min(electrical), max(electrical), 150);
y = pdf(pd, x_values);
plot3(x_values, ones(size(x_values)) * 90, y, '.black', 'LineWidth', 2);

electrical = data_e{3}(filter, EM);
pd = fitdist(electrical, 'Normal');
x_values = linspace(min(electrical), max(electrical), 150);
y = pdf(pd, x_values);
plot3(ones(size(x_values)) * 90,x_values, y, '.black', 'LineWidth', 2);

anatomical = data_a{1}(filter, AM);
pd = fitdist(anatomical, 'Normal');
x_values = linspace(min(anatomical), max(anatomical), 150);
y = pdf(pd, x_values);
plot3(x_values, ones(size(x_values)) * 90, y, '-red', 'LineWidth', 2);
anatomical = data_a{3}(filter, AM);
pd = fitdist(anatomical, 'Normal');
x_values = linspace(min(anatomical), max(anatomical), 150);
y = pdf(pd, x_values);
plot3(ones(size(x_values)) * 90,x_values, y, '-red', 'LineWidth', 2);

delta = angle_sets{1}(filter, EM_AM);
pd = fitdist(delta, 'Normal');
x_values = linspace(min(delta), max(delta), 150);
y = pdf(pd, x_values);
plot3(x_values, ones(size(x_values)) * 90, y, '--blue', 'LineWidth', 2);
delta = angle_sets{3}(filter, EM_AM);
pd = fitdist(delta, 'Normal');
x_values = linspace(min(delta), max(delta), 150);
y = pdf(pd, x_values);
plot3(ones(size(x_values)) * 90,x_values, y, '--blue', 'LineWidth', 2);


% scatter3(electrical_F_T(:, 1), electrical_F_T(:, 2), zeros(size(electrical_F_T, 1), 1), 'k', 'filled','MarkerFaceAlpha',0.12)
% scatter3(anatomical_F_T(:, 1), anatomical_F_T(:, 2), zeros(size(anatomical_F_T, 1), 1), 'r', 'filled','MarkerFaceAlpha',0.12)
% scatter3(delta_F_T(:, 1), delta_F_T(:, 2), zeros(size(delta_F_T, 1), 1), 'b', 'filled','MarkerFaceAlpha',0.12,'Marker','.')

% Legend
legend('Electrical angles density', 'Anatomical angle density', 'Angular difference density', ...
    'Electrical angles', '', 'Anatomical angles', '', 'Angular difference', ...
    'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 18)

view(3)
title('Anatomical axis: VPA, Electrical axis: maxQRS', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14)

%% 
AM = 5; 
EM = 1;
EM_AM= (EM-1)*5 + AM;
filter = logical(angle_filter(:, ceil(EM_AM/5)));
polar_coordinates_healthy = readtable('/home/mka23/Desktop/Anatomical_electrical_axes/Anatomical_electrical_analysis/polar_coordinates_results/polar_coordinates_healthy_inv.csv');


anatomical_theta_phi =  [polar_coordinates_healthy.theta_A(filter), polar_coordinates_healthy.phi_A(filter)];
electrical_theta_phi =  [polar_coordinates_healthy.theta_E(filter), polar_coordinates_healthy.phi_E(filter)];
delta_theta_phi =  [polar_coordinates_healthy.dtheta(filter), polar_coordinates_healthy.dphi(filter)];

figure; hold on; 

% Define grid for electrical data
x_electrical = linspace(min(electrical_theta_phi(:,1)), max(electrical_theta_phi(:,1)), 100);
y_electrical = linspace(min(electrical_theta_phi(:,2)), max(electrical_theta_phi(:,2)), 100);
[X_electrical, Y_electrical] = meshgrid(x_electrical, y_electrical);

% Compute Gaussian estimation for electrical data
Z_electrical = mvnpdf([X_electrical(:) Y_electrical(:)], mean(electrical_theta_phi), cov(electrical_theta_phi));

% Define grid for anatomical data
x_anatomical = linspace(min(anatomical_theta_phi(:,1)), max(anatomical_theta_phi(:,1)), 100);
y_anatomical = linspace(min(anatomical_theta_phi(:,2)), max(anatomical_theta_phi(:,2)), 100);
[X_anatomical, Y_anatomical] = meshgrid(x_anatomical, y_anatomical);

% Compute Gaussian estimation for anatomical data
Z_anatomical = mvnpdf([X_anatomical(:) Y_anatomical(:)], mean(anatomical_theta_phi), cov(anatomical_theta_phi));

% Define grid for delta data
x_delta = linspace(min(delta_theta_phi(:,1)), max(delta_theta_phi(:,1)), 100);
y_delta = linspace(min(delta_theta_phi(:,2)), max(delta_theta_phi(:,2)), 100);
[X_delta, Y_delta] = meshgrid(x_delta, y_delta);

% Compute Gaussian estimation for delta data
Z_delta = mvnpdf([X_delta(:) Y_delta(:)], mean(delta_theta_phi), cov(delta_theta_phi));

% Plot contours of density distributions
contour(X_electrical, Y_electrical, reshape(Z_electrical, size(X_electrical)),6, 'LineColor', 'k','LineWidth',1.5);
contour(X_anatomical, Y_anatomical, reshape(Z_anatomical, size(X_anatomical)),6, 'LineColor', 'r','LineWidth',1.5);
contour(X_delta, Y_delta, reshape(Z_delta, size(X_delta)),6,'LineColor','b', 'LineWidth',1.5)

% Make the figure look nice
grid on; view(45, 55);
view(3)
xlabel('$\theta$ (Frontal Plane)','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12)
ylabel('$\phi$ Transverse | Sagittal Plane','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12)
xticks(-120:30:60)
xlim([-120 60])
yticks(-30:30:150)
ylim([-30 150])
zlim([-.001 0.06])
zlabel('Probability Density','Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12)

box on

% Plot PDFs along respective axes
electrical = electrical_theta_phi(:,1);
pd = fitdist(electrical, 'Normal');
x_values = linspace(min(electrical), max(electrical), 150);
y = pdf(pd, x_values);
plot3(x_values, ones(size(x_values)) * 150, y, '-.black', 'LineWidth', 2);

electrical = electrical_theta_phi(:,2);
pd = fitdist(electrical, 'Normal');
x_values = linspace(min(electrical), max(electrical), 100);
y = pdf(pd, x_values);
plot3(ones(size(x_values)) * 60,x_values, y, '-.black', 'LineWidth', 2);

anatomical = anatomical_theta_phi(:,1);
pd = fitdist(anatomical, 'Normal');
x_values = linspace(min(anatomical), max(anatomical), 150);
y = pdf(pd, x_values);
plot3(x_values, ones(size(x_values)) * 150, y, '-red', 'LineWidth', 2);

anatomical = anatomical_theta_phi(:,2);
pd = fitdist(anatomical, 'Normal');
x_values = linspace(min(anatomical), max(anatomical), 150);
y = pdf(pd, x_values);
plot3(ones(size(x_values)) * 60,x_values, y, '-red', 'LineWidth', 2);

delta = delta_theta_phi(:,1);
pd = fitdist(delta, 'Normal');
x_values = linspace(min(delta), max(delta), 150);
y = pdf(pd, x_values);
plot3(x_values, ones(size(x_values)) * 150, y, '--blue', 'LineWidth', 2);

delta = delta_theta_phi(:,2);
pd = fitdist(delta, 'Normal');
x_values = linspace(min(delta), max(delta), 150);
y = pdf(pd, x_values);
plot3(ones(size(x_values)) * 60,x_values, y, '--blue', 'LineWidth', 2);

scatter3(electrical_theta_phi(:, 1), electrical_theta_phi(:, 2), zeros(size(electrical_theta_phi, 1), 1), 'xk','MarkerEdgeAlpha',0.05)
scatter3(anatomical_theta_phi(:, 1), anatomical_theta_phi(:, 2), zeros(size(anatomical_theta_phi, 1), 1), 'xr','MarkerEdgeAlpha',0.05)
scatter3(delta_theta_phi(:, 1), delta_theta_phi(:, 2), zeros(size(delta_theta_phi, 1), 1), 'xb', 'MarkerEdgeAlpha',0.05)

% Legend
legend('Electrical angles density', 'Anatomical angle density', 'Angular difference density', ...
    'Electrical angles', '', 'Anatomical angles', '', 'Angular difference', ...
    'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 18)

view(3)
title('Anatomical axis: VPA, Electrical axis: maxQRS', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14)


%% 
AM = 5; EM=1;
polar_coordinates_healthy = readtable('/home/mka23/Desktop/Anatomical_electrical_axes/Anatomical_electrical_analysis/polar_coordinates_results/polar_coordinates_healthy.csv');

     
x =  [polar_coordinates_healthy.theta_A(filter);polar_coordinates_healthy.theta_E(filter);polar_coordinates_healthy.dtheta(filter)];
y =  [polar_coordinates_healthy.phi_A(filter);polar_coordinates_healthy.phi_E(filter);polar_coordinates_healthy.dphi(filter)];
classifier = [repmat({'Anatomical'}, length(x)/3, 1); repmat({'Electrical'}, length(x)/3, 1);repmat({'delta'}, length(x)/3, 1)];
figure
h=scatterhist(x,y,'Group',classifier,'Kernel','on','Location','NorthWest',...
    'Direction','out','Color','rkb',...
    'LineWidth',[2,2,2],'Marker','...','MarkerSize',[2,2,2]);



%% 

x =  [data_a{1}(filter, AM);data_e{1}(filter, EM)];
y =  [data_a{2}(filter, AM);data_e{2}(filter, EM)];
classifier = [repmat({'Anatomical'}, length(x)/2, 1); repmat({'Electrical'}, length(x)/2, 1)];
figure
scatterhist(x,y,'Group',classifier,'Kernel','on','Location','NorthWest',...
    'Direction','out','Color','kbr',...
    'LineWidth',[2,2,2],'Marker','..','MarkerSize',[2,2,2]);



