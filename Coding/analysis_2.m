clc,clear;

% 自定义开始时间、结束时间和步长
start_time =0;      % 开始时间 (秒)
end_time = 414;      % 结束时间 (秒)
dt = 1;              % 时间步长 (秒)

% 设定总时长和步长，生成时间数组
time_steps = start_time:dt:end_time;
n_sections = 224;    % 总板凳数量
theta_initial = 2 * pi * 16;  % 初始角度 (第16圈)，即 theta = 32*pi
pitch = 0.55;  % 螺距

% 初始化存储位置的矩阵
positions = zeros(n_sections * 2, length(time_steps));

% 定义一个函数，输入时间t，返回t时刻的所有把手位置
function positions = calculate_positions_at_time(t)
    n_sections = 224;  % 把手数量
    pitch = 0.55;  % 螺距 55cm
    v = 1; %龙头的速度
    distance_1 = 2.86;  % 龙头两个把手的距离
    distance_2 = 1.65;  % 龙身两个把手的距离
    theta_initial = 2 * pi * 16;  % 初始角度 (第16圈)
    positions = zeros(n_sections * 2, 1);  % 存储位置数据
    
    theta = sqrt(theta_initial^2 - 4 * pi * v / pitch * t);  % 当前时刻的角度变化
    radius = pitch * theta / (2 * pi);  % 螺旋半径
    
    % 计算龙头位置
    x_head = radius * cos(theta); 
    y_head = radius * sin(theta);
    positions(1) = x_head;
    positions(2) = y_head;
    
    % 计算每节龙身的位置
    for i = 2:n_sections
        if i == 2
            distance = distance_1;
        else
            distance = distance_2;
        end
        delta_theta = distance / radius;
        theta = theta + delta_theta;
        radius = 0.55 * theta / (2 * pi);
        
        % 计算每节板凳的位置
        x_i = radius * cos(theta);
        y_i = radius * sin(theta);
        positions(2*i-1) = x_i;
        positions(2*i) = y_i;
    end
end

% 生成每秒的位置信息
for t_idx = 1:length(time_steps)
    t = time_steps(t_idx);
    positions(:, t_idx) = calculate_positions_at_time(t);
end

% 计算完整的螺旋轨迹：从 theta = 32*pi 到 theta = 0
theta_values = linspace(32*pi, 0, 1000);  % 从 32*pi 到 0 的角度变化
spiral_x = zeros(size(theta_values));
spiral_y = zeros(size(theta_values));

for i = 1:length(theta_values)
    theta = theta_values(i);
    radius = pitch * theta / (2 * pi);  % 螺旋半径
    spiral_x(i) = radius * cos(theta);
    spiral_y(i) = radius * sin(theta);
end

% 位置数据生成完成后，直接用于分析和绘图
% ---------------- 动态仿真 ------------------

figure_handle = figure;  % 获取图像句柄
hold on;
grid on;
axis equal;
xlim([-12, 12]);  % 调整显示区域
ylim([-12, 12]);  % 调整显示区域
title('板凳龙位置仿真');
xlabel('X Position (m)');
ylabel('Y Position (m)');

% 在图上绘制螺旋线轨迹 (虚线)
plot(spiral_x, spiral_y, '--k', 'LineWidth', 1);  % 用虚线显示完整的螺旋轨迹

% 动态显示板凳龙的运动
for t_idx = 1:length(time_steps)
    % 如果图像窗口已关闭，则中止动画
    if ~isvalid(figure_handle)
        disp('窗口已关闭，动画中止');
        break;
    end
    
    % 清除之前的绘图（不清除螺旋轨迹）
    cla;
    plot(spiral_x, spiral_y, '--k', 'LineWidth', 1);  % 重新绘制螺旋轨迹
    
    % 获取当前时间步的所有板凳位置信息
    x_positions = positions(1:2:end, t_idx);  % x坐标
    y_positions = positions(2:2:end, t_idx);  % y坐标
    
    % 绘制每节板凳的位置
    for i = 2:n_sections
        % 获取当前节的前后把手位置
        x_back = x_positions(i);
        y_back = y_positions(i);
        x_front = x_positions(i - 1);
        y_front = y_positions(i - 1);
        
        % 计算方向向量
        dx = x_front - x_back;
        dy = y_front - y_back;
        length_vector = sqrt(dx^2 + dy^2);
        if length_vector == 0
            continue;  % 如果前后把手重合，则跳过当前节
        end
        direction_x = dx / length_vector;
        direction_y = dy / length_vector;
        
        % 垂直方向向量
        perpendicular_x = -direction_y;
        perpendicular_y = direction_x;
        
        % 计算矩形四个角的位置
        board_width = 0.3;
        handle_offset = 0.275;
        
        x_back_left = x_back + (board_width / 2) * perpendicular_x - handle_offset * direction_x;
        y_back_left = y_back + (board_width / 2) * perpendicular_y - handle_offset * direction_y;
        x_back_right = x_back - (board_width / 2) * perpendicular_x - handle_offset * direction_x;
        y_back_right = y_back - (board_width / 2) * perpendicular_y - handle_offset * direction_y;
        x_front_left = x_front + (board_width / 2) * perpendicular_x + handle_offset * direction_x;
        y_front_left = y_front + (board_width / 2) * perpendicular_y + handle_offset * direction_y;
        x_front_right = x_front - (board_width / 2) * perpendicular_x + handle_offset * direction_x;
        y_front_right = y_front - (board_width / 2) * perpendicular_y + handle_offset * direction_y;
        
        % 绘制矩形 (表示板凳)
        fill([x_back_left, x_back_right, x_front_right, x_front_left], ...
             [y_back_left, y_back_right, y_front_right, y_front_left], 'r');
    end
    
    % 显示当前时间
    current_time = time_steps(t_idx);  % 获取当前时间
    time_text = sprintf('Time: %.1f s', current_time);  % 格式化时间文本
    text(-11, 11, time_text, 'FontSize', 12, 'FontWeight', 'bold');  % 在图像左上角显示时间
    % 暂停一段时间来创建动画效果
    pause(0.1);
end

% Hold the final plot
hold off;
