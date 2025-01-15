function [axis_vector] = VCG_axis_maxNorm(xyz_data, t_Qonset, t_Qoffset, PLOT)

    % Axis -- largest distance between end of QRS loop and any point on the loop 
    QRS_start = xyz_data(:,t_Qonset + 1); % +1 for true isolelectric start
    
    Loop_data = xyz_data(:,t_Qonset:t_Qoffset);

    dist = [];
    for i = 2:length(Loop_data)
        dist(i - 1) = norm(QRS_start - Loop_data(:,i));
    end
    
    [dist_max, i_max] = max(dist);
    t_max = i_max+1; %time stamp at which max occurs
    axis_vector = Loop_data(:,t_max) - QRS_start;
    % 
    if PLOT == true 
    figure(), plot3(Loop_data(1,:),Loop_data(2,:),Loop_data(3,:)), xlabel('X Amplitude (mV)'), ylabel('Y Amplitude (mV)'), zlabel('Z Amplitude (mV)'), hold on
    scatter3(QRS_start(1), QRS_start(2), QRS_start(3), 'red', 'filled'); legend('QRS Loop'), hold on
    scatter3(Loop_data(1,t_max), Loop_data(2,t_max), Loop_data(3,t_max), 'green', 'filled'), hold on
    plot3([QRS_start(1) Loop_data(1,t_max)],[QRS_start(2) Loop_data(2,t_max)],[QRS_start(3) Loop_data(3,t_max)])
    legend('QRS Loop','QRS start','Point of max normal distance','Electrical Axis')
    grid on
    
    figure(), subplot(311),plot(xyz_data(1,:)), hold on, scatter(t_max,axis_vector(1),'red'), ylabel('X Amplitude (mV'), grid on,...
                subplot(312),plot(xyz_data(2,:)), hold on, scatter(t_max,axis_vector(2)),ylabel('Y Amplitude (mV'),grid on, ...
                subplot(313),plot(xyz_data(3,:)), hold on, scatter(t_max,axis_vector(3)),ylabel('Z Amplitude (mV'),grid on;
    end
end


