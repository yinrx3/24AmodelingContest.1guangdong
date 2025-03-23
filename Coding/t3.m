clc; clear;
start_time = 0; %开始时间
dt = 0.01; % 默认时间步长
end_time = 440;%结束时间
n_sections = 224; % 把手数量
width = 0.3; % 板凳宽度30cm
offset = 0.275; % 把手到边界的距离为 0.275m
collision_time = -1; % 初始化碰撞时间
initial_pitch = 0.4; % 螺距初始值
pitch_step = 0.01; % 螺距增加的步长
radius_threshold = 4.5; % 半径阈值
max_pitch = 1.0; % 螺距最大值，防止无限循环

% 时间数组
time_steps = start_time:dt:end_time;

% 定义函数计算每时刻所有把手位置
function positions = calculate_positions_at_time(t, pitch)
    n_sections = 224;
    v = 1; %龙头的速度
    distance_1 = 2.86; %龙头两个把手的距离
    distance_2 = 1.65; %龙身两个把手的距离
    theta_initial = 2 * pi * 16; % 初始角度
    positions = zeros(n_sections * 2, 1);
    
    theta = sqrt(theta_initial^2 - 4 * pi * v/ pitch * t); % 当前角度
    radius = pitch * theta / (2 * pi); % 螺旋半径
    
    x_head = radius * cos(theta); 
    y_head = radius * sin(theta);
    positions(1) = x_head;
    positions(2) = y_head;
    
    for i = 2:n_sections
        if (i == 2)
            distance = distance_1;
        else
            distance = distance_2;
        end
        delta_theta = distance / radius;
        theta = theta + delta_theta;
        radius = pitch * theta / (2 * pi);
        x_i = radius * cos(theta);
        y_i = radius * sin(theta);
        positions(2*i-1) = x_i;
        positions(2*i) = y_i;
    end
end

% 定义辅助函数，检查点是否在矩形内
function inside = is_point_in_rectangle(px, py, rect)
    x1 = rect(1,1); y1 = rect(1,2);
    x2 = rect(2,1); y2 = rect(2,2);
    x3 = rect(3,1); y3 = rect(3,2);
    x4 = rect(4,1); y4 = rect(4,2);
    
    % 使用向量叉积判断点是否在矩形内
    v1 = [x2 - x1, y2 - y1]; % 向量 v1
    v2 = [x3 - x2, y3 - y2]; % 向量 v2
    v3 = [x4 - x3, y4 - y3]; % 向量 v3
    v4 = [x1 - x4, y1 - y4]; % 向量 v4
    
    % 点相对于每条边的叉积
    cp1 = (px - x1)*(y2 - y1) - (py - y1)*(x2 - x1); % 点在边1左侧
    cp2 = (px - x2)*(y3 - y2) - (py - y2)*(x3 - x2); % 点在边2左侧
    cp3 = (px - x3)*(y4 - y3) - (py - y3)*(x4 - x3); % 点在边3左侧
    cp4 = (px - x4)*(y1 - y4) - (py - y4)*(x1 - x4); % 点在边4左侧
    
    % 如果所有叉积都为正或都为负，则点在矩形内
    inside = (cp1 >= 0 && cp2 >= 0 && cp3 >= 0 && cp4 >= 0) || ...
             (cp1 <= 0 && cp2 <= 0 && cp3 <= 0 && cp4 <= 0);
end

% 初始化最小螺距
found_pitch = -1;

% 开始调整 pitch 直到找到最小满足条件的 pitch
pitch = initial_pitch;
while pitch <= max_pitch
    collision_time = -1; % 重置碰撞时间

    % 主循环，检测碰撞
    for t_idx = 1:length(time_steps)
        t = time_steps(t_idx);
        positions_now = calculate_positions_at_time(t, pitch);
        % 检查螺旋中心，当 theta 接近 0 时，停止计算
        theta_2 = (2 * pi * 16)^2 - 4 * pi / 0.55 * t;
        if theta_2 < 0
            disp(['龙头在 t = ', num2str(t), ' 秒到达螺旋中心，无法继续旋入']);
            collision_time = t;
            break;
        end

        % 计算龙头的两个前角，基于前后把手的位置
        x_head_front = positions_now(1);   % 龙头前把手的 x 坐标
        y_head_front = positions_now(2);   % 龙头前把手的 y 坐标
        x_head_back = positions_now(3);    % 龙头后把手的 x 坐标
        y_head_back = positions_now(4);    % 龙头后把手的 y 坐标

        % 计算龙头的方向向量
        dx_head = x_head_front - x_head_back;
        dy_head = y_head_front - y_head_back;
        length_vector_head = sqrt(dx_head^2 + dy_head^2);

        % 归一化方向向量 (使其长度为1)
        direction_x_head = dx_head / length_vector_head;
        direction_y_head = dy_head / length_vector_head;

        % 垂直方向向量 (用于计算龙头的左右两侧)
        perpendicular_x_head = -direction_y_head;
        perpendicular_y_head = direction_x_head;

        % 计算龙头的前左和后左角坐标，考虑宽度和偏移量
        x_head_front_left = x_head_front + (width / 2) * perpendicular_x_head + offset * direction_x_head;
        y_head_front_left = y_head_front + (width / 2) * perpendicular_y_head + offset * direction_y_head;
    
        x_head_back_left = x_head_front + (width / 2) * perpendicular_x_head - offset * direction_x_head;
        y_head_back_left = y_head_front + (width / 2) * perpendicular_y_head - offset * direction_y_head;
        % 遍历之后的所有板凳，检查是否发生碰撞
        for i = 3:n_sections
            % 当前节的后把手和前把手的位置
            x_back = positions_now(2*i-1); % 后把手
            y_back = positions_now(2*i);   % 后把手
            x_front = positions_now(2*(i-1)-1); % 前把手
            y_front = positions_now(2*(i-1));   % 前把手

            % 计算板凳的方向向量 (dx, dy)
            dx = x_back - x_front; % 计算前把手到后把手的向量
            dy = y_back - y_front;
            length_vector = sqrt(dx^2 + dy^2);

            % 归一化方向向量 (使其长度为1)
            direction_x = dx / length_vector;
            direction_y = dy / length_vector;

            % 垂直方向向量 (用于计算板凳的左右两侧)
            perpendicular_x = -direction_y;
            perpendicular_y = direction_x;

            % 计算板凳的四个角，考虑偏移量 offset 和宽度 width
            % 后右角
            x_back_left = x_back + (width / 2) * perpendicular_x + offset * direction_x;
            y_back_left = y_back + (width / 2) * perpendicular_y + offset * direction_y;

            % 后左角
            x_back_right = x_back - (width / 2) * perpendicular_x + offset * direction_x;
            y_back_right = y_back - (width / 2) * perpendicular_y + offset * direction_y;

            % 前右角
            x_front_left = x_front + (width / 2) * perpendicular_x - offset * direction_x;
            y_front_left = y_front + (width / 2) * perpendicular_y - offset * direction_y;

            % 前左角
            x_front_right = x_front - (width / 2) * perpendicular_x - offset * direction_x;
            y_front_right = y_front - (width / 2) * perpendicular_y - offset * direction_y;

            % 将这四个角点存储为矩形
            rect = [x_back_left, y_back_left;
                    x_back_right, y_back_right;
                    x_front_right, y_front_right;
                    x_front_left, y_front_left];

            % 检查龙头前左、后左是否进入该矩形内
            if is_point_in_rectangle(x_head_front_left, y_head_front_left, rect) || ...
               is_point_in_rectangle(x_head_back_left, y_head_back_left, rect)
                collision_time = t;
                disp(['和第',num2str(i-2),'节龙身碰撞'])
                disp(['碰撞发生在 t = ', num2str(collision_time), ' 秒']);
                break;
            end
        end
        if collision_time > 0
            break;
        end
    end
    
    % 如果发生碰撞，检查龙头是否位于目标圆内
    if collision_time > 0
        x_head = positions_now(1);   % 龙头前把手的 x 坐标
        y_head = positions_now(2);   % 龙头前把手的 y 坐标
        
        if sqrt(x_head^2 + y_head^2) <= radius_threshold
            found_pitch = pitch;  % 找到满足条件的最小螺距
            disp(['找到了满足条件的最小螺距：pitch = ', num2str(found_pitch)]);
            break;
        end
    end
    
    % 增加 pitch
    pitch = pitch + pitch_step;
end

% 如果没有找到满足条件的螺距，输出提示
if found_pitch == -1
    disp('未找到使得龙头位于圆内的最小螺距');
end
