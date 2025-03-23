clc, clear;
time_total = 423; % 总时间 300秒
dt = 1; % 时间步长 1秒
dt_small = 0.0001; % 用于计算速度的小时间步长
n_sections = 224; % 把手数量
v = 1; %龙头的速度

% 时间数组
time_steps = 0:dt:time_total;

% 初始化存储位置的矩阵
positions = zeros(n_sections * 2, length(time_steps)); 
velocity = zeros(n_sections, length(time_steps)); % 因为速度在相邻的整数秒之间计算

% 定义一个函数，输入时间t，返回t时刻的所有把手位置
function positions = calculate_positions_at_time(t)
    % 设定参数
    n_sections = 224; % 把手数量
    pitch = 0.55; % 螺距55cm
    v = 1; %龙头的速度
    distance_1 = 2.86; %龙头两个把手的距离
    distance_2 = 1.65; %龙身两个把手的距离
    theta_initial = 2 * pi * 16; % 第16圈起始角度
    % 初始化存储位置的矩阵
    positions = zeros(n_sections * 2, 1);
    
    % 计算当前时刻的龙头角度和半径
    theta = sqrt(theta_initial^2 - 4 * pi *v / pitch * t); % 当前时刻的角度变化
    radius = pitch * theta / (2 * pi); % 螺线半径
    
    % 计算龙头的位置
    x_head = radius * cos(theta); % x 位置
    y_head = radius * sin(theta); % y 位置
    
    % 将龙头位置保存到矩阵（第1行和第2行分别为x和y）
    positions(1) = x_head;
    positions(2) = y_head;
    
    % 计算每一节龙身的位置
    for i = 2:n_sections
        %两个把手之间的距离
        if (i == 2)
            distance = distance_1;
        else
            distance = distance_2;
        end
        delta_theta = distance / radius;
        theta = theta + delta_theta;
        radius = 0.55 * theta / (2 * pi);
        % 计算该节板凳的 x 和 y 坐标
        x_i = radius * cos(theta);
        y_i = radius * sin(theta);
        
        % 保存该节板凳的位置 (注意行号是 2i-1 和 2i)
        positions(2*i-1) = x_i;
        positions(2*i) = y_i;
    end
end

% 计算每秒钟的位置信息
for t_idx = 1:length(time_steps)
    t = time_steps(t_idx);
    positions(:, t_idx) = calculate_positions_at_time(t);
end

% 计算整数秒之间的速度，使用 t 和 t+0.0001 之间的位移来计算
for t_idx = 1:length(time_steps)
    % 当前时刻 t_now 和 t+0.001 秒时刻
    t_now = time_steps(t_idx);
    t_small = t_now + dt_small;
    
    % 获取所有把手在 t_now 时刻的位置
    pos_now = calculate_positions_at_time(t_now);  % 返回一个 2*n_sections 大小的向量
    
    % 获取所有把手在 t+0.001s 时刻的位置
    pos_next = calculate_positions_at_time(t_small);  % 返回一个 2*n_sections 大小的向量
    
    % 将 pos_now 和 pos_next 重新分为 x, y 坐标
    x_now = pos_now(1:2:end);   % 当前时刻所有把手的 x 坐标
    y_now = pos_now(2:2:end);   % 当前时刻所有把手的 y 坐标
    x_next = pos_next(1:2:end); % t+0.001s 时刻所有把手的 x 坐标
    y_next = pos_next(2:2:end); % t+0.001s 时刻所有把手的 y 坐标
    
    % 计算速度 (vx, vy)，然后计算总速度
    v_x = (x_next - x_now) / dt_small;
    v_y = (y_next - y_now) / dt_small;
    
    % 计算每个把手的速度并存储
    velocity(:, t_idx) = sqrt(v_x.^2 + v_y.^2);
end


%将结果保留六位小数
positions = round(positions, 6);
velocity = round(velocity, 6);

% 构建第一页行头，表示每节板凳的 x(m) 和 y(m) 坐标
headers_rows_1 = {'', '龙头x(m)', '龙头y(m)'};
for i = 1:n_sections-3
    headers_rows_1 = [headers_rows_1, {['第' num2str(i) '节龙身x(m)'], ['第' num2str(i) '节龙身y(m)']}];
end
headers_rows_1 = [headers_rows_1, {'龙尾x(m)', '龙尾y(m)'}];
headers_rows_1 = [headers_rows_1, {'龙尾（后）x(m)', '龙尾（后）y(m)'}];

% 构建第一页列头，表示每个时间点 (0s, 1s, 2s, ...)
headers_columns_1 = arrayfun(@(x) [num2str(x) 's'], time_steps, 'UniformOutput', false);

% 将 headers_rows_1 和 headers_columns_1 拼接到数据中
final_data_1 = [headers_rows_1', [headers_columns_1; num2cell(positions)]];

% 构建第二页行头
headers_rows_2 = {'', '龙头（m/s）'};
for i = 1:n_sections-3
    headers_rows_2 = [headers_rows_2, {['第' num2str(i) '节龙身（m/s）']}];
end
headers_rows_2 = [headers_rows_2, {'龙尾（m/s）'}];
headers_rows_2 = [headers_rows_2, {'龙尾（后）（m/s）'}];

% 构建第二页列头，表示每个时间点 (1s, 2s, 3s, ...)
headers_columns_2 = arrayfun(@(x) [num2str(x) 's'], time_steps(1:end), 'UniformOutput', false);

% 将 headers_rows_2 和 headers_columns_2 拼接到数据中
final_data_2 = [headers_rows_2', [headers_columns_2; num2cell(velocity)]];

% 保存到Excel文件，位置数据到第一页，速度数据到第二页
writecell(final_data_1, 'result1.xlsx', 'Sheet', '位置');
writecell(final_data_2, 'result1.xlsx', 'Sheet', '速度');
