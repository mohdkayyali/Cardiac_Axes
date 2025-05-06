function [Loop_area] = Calc_area(xyz_data, start_time,end_time,method)
%can add another variable for different ways of calculating the area
    
    %% Using sum of triangle areas
    QRS_loop = xyz_data(:,start_time:end_time);
    N = size(QRS_loop, 2);
    Loop_area = 0;    
    
    if strcmp(method,'Qstart_O') 
        P0 = xyz_data(:,start_time+1);  % Loop Origin

        % Loop through the triangles made from consecutive vectors
        for i = 1:N
            % Define vertices of the triangle 
            P1 = QRS_loop(:, i);
            if i < N
                P2 = QRS_loop(:, i + 1);
            else
                P2 = QRS_loop(:, 1);  % Connect the last point to the first point
            end
            
            Area_triangle = 0.5* norm(cross(P1 - P0, P2 - P0));
            Loop_area = Loop_area + Area_triangle;
    
        end
    elseif strcmp(method,'mean_O')
        P0 = mean(QRS_loop,2);
        for i = 1:N
            % Define vertices of the triangle 
            P1 = QRS_loop(:, i);
            if i < N
                P2 = QRS_loop(:, i + 1);
            end
            
            Area_triangle = 0.5* norm(cross(P1 - P0, P2 - P0));
            Loop_area = Loop_area + Area_triangle;
    
        end

    end
end


%% plot
% QRS_start = xyz_data(:,t_Qonset+1);
% figure(), plot3(xyz_data(1,t_Qonset:t_Qoffset),xyz_data(2,t_Qonset:t_Qoffset),xyz_data(3,t_Qonset:t_Qoffset), linewidth = 2), xlabel('X Amplitude (mV)'), ylabel('Y Amplitude (mV)'), zlabel('Z Amplitude (mV)'), hold on
% scatter3(QRS_start(1), QRS_start(2), QRS_start(3), 'red', 'filled');  hold on
% for i = 1:N
%     if i < N
%         plot3([QRS_start(1) xyz_data(1,t_Qonset+i)],[QRS_start(2) xyz_data(2,t_Qonset+i)],[QRS_start(3) xyz_data(3,t_Qonset+i)],'black'),hold on
%     else
%         plot3([QRS_start(1) xyz_data(1,t_Qonset)],[QRS_start(2) xyz_data(2,t_Qonset)],[QRS_start(3) xyz_data(3,t_Qonset)],'black'),hold on
%     end
% end
% plot3([QRS_start(1) axis_vector_M1(1)],[QRS_start(2) axis_vector_M1(2)],[QRS_start(3) axis_vector_M1(3)], 'red', linewidth = 2)
% 
% grid on

% QRS_start = xyz_data(:,t_Qonset+1);end_time=t_Qoffset;start_time = t_Qonset;
% figure(), plot3(xyz_data(1,t_Qonset:t_Qoffset),xyz_data(2,t_Qonset:t_Qoffset),xyz_data(3,t_Qonset:t_Qoffset), linewidth = 2), xlabel('X Amplitude (mV)'), ylabel('Y Amplitude (mV)'), zlabel('Z Amplitude (mV)'), hold on
% scatter3(P0(1), P0(2), P0(3), 'red', 'filled');  hold on
% for i = 1:N
%     if i < N
%         plot3([P0(1) xyz_data(1,t_Qonset+i)],[P0(2) xyz_data(2,t_Qonset+i)],[P0(3) xyz_data(3,t_Qonset+i)],'black'),hold on
%     end
% end
% plot3([P0(1) axis_vector_M1(1)],[P0(2) axis_vector_M1(2)],[P0(3) axis_vector_M1(3)], 'red', linewidth = 2)
% 
% grid on