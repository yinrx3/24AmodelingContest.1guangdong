clc, clear;
start_time = 0; % 从进入调头空间开始计算
dt = 1; % 时间步长 1秒
end_time = 100;%结束时间足够大，确保不会遗漏
dt_small = 1; % 用于计算速度的小时间步长
n_sections = 224; % 把手数量

% 时间数组
time_steps = start_time:dt:end_time;

%定义一个函数，用于计算任意时刻所有把手的坐标
function positions = calculate_positions_at_time(t,v)
    % 参数设定
    pitch = 1.7; %螺距
    n_sections = 224;      % 总共224节
    distance_1 = 2.86;     % 龙头前两个把手的距离
    distance_2 = 1.65;     % 龙身和龙尾每两个把手的距离
    R_space = 4.5;         % 调头空间半径
    positions = zeros(n_sections * 2, 1);  % 初始化存储位置的矩阵
    
    % 计算龙头前把手的位置
    [x_head, y_head] = calculate_head_position(t,v);

    % 从龙头开始，依次计算每个把手的位置
    current_x = x_head;
    current_y = y_head;
    previous_x = x_head;
    previous_y = y_head;

    for i = 1:n_sections
        % 根据当前是第一个把手还是后续把手确定距离
        if i == 1
            distance = distance_1;  % 如果是龙头前把手，使用龙头前两个把手的距离
        else
            distance = distance_2;  % 龙身和龙尾两个把手的距离
        end
        % 判断当前把手是否位于调头空间内
        if sqrt(current_x^2 + current_y^2) <= R_space
            % 如果当前把手位于调头空间内
            if i>=2 && sqrt(previous_x^2 + previous_y^2) >R_space
                % 如果是龙身把手且上一把手在调头空间外，代表当前位置不准确，需要更新当前位置
                if i==2
                    [current_x, current_y,~] = search_on_arc(previous_x, previous_y, distance_1);
                else
                    [current_x, current_y,~] = search_on_arc(previous_x, previous_y, distance_2);
                end
            end
            %如果是最后一个节点，直接存储
            if (i==n_sections)
                break;
            end
            %在圆弧内计算下一把手
            [next_x, next_y, found] = search_on_arc(current_x, current_y, distance);

            %如果在圆弧内没有找到
            if(~found)
                [next_x, next_y] = search_on_spiral(current_x, current_y, distance);
            end

        else
            % 如果当前把手不在调头空间内，按照近似方式计算下一把手
            current_r = sqrt(current_x^2 + current_y^2);
            current_theta = 2 * pi * current_r / pitch;
            delta_theta = distance / current_r;

            %如果在盘入螺旋线
            if(cos(current_theta) * current_x >= 0 && sin(current_theta) * current_y >= 0)
                next_theta = current_theta + delta_theta;
                next_r = pitch * next_theta / (2 * pi);
                next_x = next_r * cos(next_theta);
                next_y = next_r * sin(next_theta);

            %如果在盘出螺旋线
            else
                next_theta = current_theta - delta_theta;
                next_r = pitch * next_theta / (2 * pi);
                next_x = -next_r * cos(next_theta);
                next_y = -next_r * sin(next_theta);
            end
        end

        %记录当前位置并更新坐标
        positions(2*i-1) = current_x;
        positions(2*i) = current_y;
        previous_x = current_x;
        previous_y = current_y;
        current_x = next_x;
        current_y = next_y;

    end
end

% 计算龙头位置的函数
function [x_head, y_head] = calculate_head_position(t,v)
    % 参数设定
    pitch = 1.7;  % 螺距
    r1 = 2.1118;  % 调头弧线1的半径
    r2 = 2.3963;  % 调头弧线2的半径
    theta1 = 3.0215;  % 第一段圆弧旋转的角度
    theta2 = 3.0215;  % 第二段圆弧旋转的角度
    T_turn = (r1 * theta1 + r2 * theta2) / v;  % 调头阶段总时间
    
    % 1. 螺旋线阶段，适用于 -100s 到 0s
    if t <= 0
        theta_initial_in = 16.6320;  % 螺旋线初始角度
        theta = sqrt(theta_initial_in^2 - 4 * pi * v / pitch * t); % 角度
        radius = pitch * theta / (2 * pi); % 半径
        x_head = radius * cos(theta); % x坐标
        y_head = radius * sin(theta); % y坐标
    
    % 2. 调头阶段，适用于 0s 到 T_turn
    elseif t > 0 && t <= T_turn
        % 第一段圆弧运动：0 到 r1 * theta1 时间内
        if t <= (r1 * theta1)/v
            % 圆心坐标
            x_in_circle_center = -1.3404;
            y_in_circle_center = -1.9852;

            % 圆弧和盘入螺旋线的交点
            x_in_intersection = -2.7119;
            y_in_intersection = -3.5911;

            % 当前时间对应的角度
            theta_arc1 = t / r1;  
            % 从圆心到交点的向量方向
            theta_start_arc1 = atan2(y_in_intersection - y_in_circle_center, ...
                                     x_in_intersection - x_in_circle_center);
                                 
            % 计算龙头在第一段圆弧上的位置
            x_head = x_in_circle_center + r1 * cos(theta_start_arc1 - theta_arc1);  % 顺时针旋转
            y_head = y_in_circle_center + r1 * sin(theta_start_arc1 - theta_arc1);  % 顺时针旋转
            
        % 第二段圆弧运动：从 r1 * theta1 到 T_turn
        else
            t_remaining = t - r1 * theta1 / v;  % 剩余时间
            x_out_circle_center = 1.1556;
            y_out_circle_center = 1.7986;
            
            % 两段圆弧交点坐标
            x_intersection = -0.1711;
            y_intersection = -0.2267;
            
            theta_arc2 = t_remaining * v / r2;  % 当前时间对应的圆弧角度
            % 从交点到圆心的向量方向
            theta_start_arc2 = atan2(y_intersection - y_out_circle_center, ...
                                     x_intersection - x_out_circle_center);
                                 
            % 计算龙头在第二段圆弧上的位置 (逆时针旋转 theta_arc2)
            x_head = x_out_circle_center + r2 * cos(theta_start_arc2 + theta_arc2);  
            y_head = y_out_circle_center + r2 * sin(theta_start_arc2 + theta_arc2);  
        end
    
    % 3. 盘出螺旋线阶段，适用于 T_turn 到 100s
    else
        t_out = t - T_turn;  % 调头结束后的时间
        [x_head, y_head] = calculate_head_position(-t_out, v);
        x_head = -x_head;
        y_head = -y_head;
        
    end
end

%搜索在螺旋线上的下一个点
function [next_x, next_y] = search_on_spiral(x, y, distance)
    % 初始theta值为16.6320
    theta = 16.6320;
    
    % 设定一个小的角度步长
    delta_theta = 0.00001;  % 可以根据需求调整步长的大小
    while true
        % 根据螺旋线方程 r = 1.7 * theta / (2 * pi) 计算当前的半径
        current_r = 1.7 * theta / (2 * pi);
        % 根据当前半径和theta计算螺旋线上对应的x, y坐标
        current_x = current_r * cos(theta);
        current_y = current_r * sin(theta);
        
        % 计算螺旋线上的点与给定点(x, y)的距离
        d = sqrt((current_x - x)^2 + (current_y - y)^2);
        
        % 如果距离达到了目标距离，则返回当前螺旋线上的点
        if abs(d - distance) < 0.5
            next_x = current_x;
            next_y = current_y;
            break;
        end
        
        % 增加theta，继续沿着螺旋线移动
        theta = theta + delta_theta;

    end
end


%搜索在圆弧上的下一个把手
function [next_x, next_y, found] = search_on_arc(current_x, current_y, distance)
    % 圆弧的参数设定
    r1 = 2.1118;  % 第一段圆弧的半径
    r2 = 2.3963;  % 第二段圆弧的半径
    x_in_circle_center = -1.3404;
    y_in_circle_center = -1.9852;
    x_out_circle_center = 1.1556;
    y_out_circle_center = 1.7986;

    % 圆弧交点和圆弧与调头空间交点坐标
    x_intersection = -0.1711;
    y_intersection = -0.2267;
    x_start_arc2 = 2.7119;
    y_start_arc2 = 3.5911;
    x_start_arc1 = -2.7119;
    y_start_arc1 = -3.5911;

    % 初始角度步长
    angle_step = 0.00001;  % 使用小步长搜索
    
    % 1. 检查当前点是否位于第二段圆弧上
    dist_to_out_center = sqrt((current_x - x_out_circle_center)^2 + (current_y - y_out_circle_center)^2);
    if abs(dist_to_out_center - r2) < 0.5
        % 计算当前点的极角
        theta_current = atan2(current_y - y_out_circle_center, current_x - x_out_circle_center);
        % 计算第二段圆弧的角度范围
        theta_start_arc2 = atan2(y_start_arc2 - y_out_circle_center, x_start_arc2 - x_out_circle_center);
        theta_end_arc2 = atan2(y_intersection - y_out_circle_center, x_intersection - x_out_circle_center);

        % 判断当前点的角度是否在第二段圆弧的角度范围内
        if theta_current >= theta_end_arc2 && theta_current <= theta_start_arc2
            %先在第二段圆弧上搜索
            [next_x, next_y, found] = findNextArcPoint( ...
                                        current_x, current_y, ...                % 当前点的坐标
                                        x_out_circle_center, y_out_circle_center, ... % 外圆的圆心坐标
                                        r2, ...                                  % 外圆的半径
                                        current_x, current_y, ...          % 第二段圆弧的起点坐标
                                        x_intersection, y_intersection, ...      % 圆弧终点交点坐标
                                        distance, ...                            % 目标距离
                                        angle_step);                              % 角度步长
            if (found)
                return;
            end
            %再在第二段圆弧上搜索
            [next_x, next_y, found] = findNextArcPoint_r( ...
                                        current_x, current_y, ...                % 当前点的坐标
                                        x_in_circle_center, y_in_circle_center, ... % 内圆的圆心坐标
                                        r1, ...                                  % 内圆的半径
                                        x_start_arc1, y_start_arc1, ...          % 第一段圆弧的起点坐标
                                        x_intersection, y_intersection, ...      % 圆弧终点交点坐标
                                        distance, ...                            % 目标距离
                                        angle_step);
            return;
        end
    end
    
    % 2. 检查当前点是否位于第一段圆弧上
    dist_to_in_center = sqrt((current_x - x_in_circle_center)^2 + (current_y - y_in_circle_center)^2);
    if abs(dist_to_in_center - r1) < 0.5
        % 计算当前点的极角
        theta_current = atan2(current_y - y_in_circle_center, current_x - x_in_circle_center);
        % 计算第一段圆弧的角度范围
        theta_start_arc1 = atan2(y_start_arc1 - y_in_circle_center, x_start_arc1 - x_in_circle_center) + 2*pi;
        theta_end_arc1 = atan2(y_intersection - y_in_circle_center, x_intersection - x_in_circle_center);

        % 判断当前点的角度是否在第一段圆弧的角度范围内
        if theta_current >= theta_end_arc1 && theta_current <= theta_start_arc1
            %在第一段圆弧上搜索
            [next_x, next_y, found] = findNextArcPoint_r( ...
                                        current_x, current_y, ...                % 当前点的坐标
                                        x_in_circle_center, y_in_circle_center, ... % 内圆的圆心坐标
                                        r1, ...                                  % 内圆的半径
                                        x_start_arc1, y_start_arc1, ...          % 第一段圆弧的起点坐标
                                        current_x, current_y, ...      % 圆弧终点交点坐标
                                        distance, ...                            % 目标距离
                                        angle_step);                              % 角度步长
            if (found)
                return
            end
            found = false;
            next_x = NaN;
            next_y = NaN;
            return;
        end
    end
    [next_x, next_y, found] = findNextArcPoint( ...
        current_x, current_y, ...                % 当前点的坐标
        x_out_circle_center, y_out_circle_center, ... % 外圆的圆心坐标
        r2, ...                                  % 外圆的半径
        x_start_arc2, y_start_arc2, ...          % 第二段圆弧的起点坐标
        x_intersection, y_intersection, ...      % 圆弧终点交点坐标
        distance, ...                            % 目标距离
        angle_step);                              % 角度步长

end

%顺时针搜索圆弧上的下一个点
function [next_x, next_y, found] = findNextArcPoint( ...
    current_x, current_y, ...                % 当前点的坐标
    x_out_circle_center, y_out_circle_center, ... % 外圆的圆心坐标
    r2, ...                                  % 外圆的半径
    x_start_arc2, y_start_arc2, ...          % 第二段圆弧的起点坐标
    x_intersection, y_intersection, ...      % 圆弧终点交点坐标
    distance, ...                            % 目标距离
    angle_step)                              % 角度步长

    % 初始化搜索起始角度
    theta_current = atan2(y_start_arc2 - y_out_circle_center, x_start_arc2 - x_out_circle_center);
    % 初始化返回值
    found = false;
    next_x = NaN;
    next_y = NaN;

    % 搜索圆弧上的点
    while theta_current > atan2(y_intersection - y_out_circle_center, x_intersection - x_out_circle_center)
        % 计算当前圆弧点的位置
        next_x = x_out_circle_center + r2 * cos(theta_current);
        next_y = y_out_circle_center + r2 * sin(theta_current);
        
        % 计算该点与 current_x, current_y 的距离
        dist = sqrt((next_x - current_x)^2 + (next_y - current_y)^2);
        if abs(dist - distance) < 0.5
            % 找到符合条件的点，返回结果
            found = true;
            return;
        end
        
        % 没找到则继续沿着圆弧向内移动，减小角度
        theta_current = theta_current - angle_step;
    end
end


%逆时针搜索
function [next_x, next_y, found] = findNextArcPoint_r( ...
    current_x, current_y, ...                % 当前点的坐标
    x_in_circle_center, y_in_circle_center, ... % 内圆的圆心坐标
    r1, ...                                  % 内圆的半径
    x_start_arc1, y_start_arc1, ...          % 第一段圆弧的起点坐标
    x_intersection, y_intersection, ...      % 圆弧终点交点坐标
    distance, ...                            % 目标距离
    angle_step)                              % 角度步长

    % 初始化搜索起始角度
    theta_current = atan2(y_intersection - y_in_circle_center, x_intersection - x_in_circle_center);

    % 搜索圆弧上的点（逆时针）
    while theta_current < atan2(y_start_arc1 - y_in_circle_center, x_start_arc1 - x_in_circle_center) + 2*pi
        % 计算当前圆弧点的位置
        next_x = x_in_circle_center + r1 * cos(theta_current);
        next_y = y_in_circle_center + r1 * sin(theta_current);
        
        % 计算该点与 current_x, current_y 的距离
        dist = sqrt((next_x - current_x)^2 + (next_y - current_y)^2);
        if abs(dist - distance) < 0.5
            % 找到符合条件的点，返回结果
            found = true;
            return;
        end
        
        % 没找到则继续沿着圆弧向外移动，增加角度
        theta_current = theta_current + angle_step;
    end

    % 如果没有找到符合条件的点，返回 false
    found = false;
    next_x = NaN;
    next_y = NaN;
end


% 定义初始龙头速度和步长
v_head = 1.00;  % 龙头速度初始为1m/s
v_step = 0.01;  % 每次增加的速度步长
max_v = 2;  % 把手速度的最大值
R_space = 4.5;  % 调头空间半径

% 标记是否找到最大速度
found_max_speed = false;

% 开始迭代速度，从1m/s逐步增加
while ~found_max_speed
    % 初始化时间数组
    time_steps = start_time:dt:end_time;
    fprintf("当前速度:%0.2f m/s\n", v_head);

    % 循环时间步长
    for t_idx = 1:length(time_steps)
        % 当前时刻 t_now 和 t+0.001 秒时刻
        t_now = time_steps(t_idx);
        fprintf("当前时刻:%0.2f s\n", t_now);
        delta_t = dt_small;
        if t_now == 0
            delta_t = 1;
        end
        t_small = t_now + delta_t;

        % 获取所有把手在 t_now 时刻的位置
        pos_now = calculate_positions_at_time(t_now, v_head);  % 返回一个 2*n_sections 大小的向量

        % 获取所有把手在 t+0.001s 时刻的位置
        pos_next = calculate_positions_at_time(t_small, v_head);  % 返回一个 2*n_sections 大小的向量

        % 将 pos_now 和 pos_next 重新分为 x, y 坐标
        x_now = pos_now(1:2:end);   % 当前时刻所有把手的 x 坐标
        y_now = pos_now(2:2:end);   % 当前时刻所有把手的 y 坐标
        x_next = pos_next(1:2:end); % t+0.001s 时刻所有把手的 x 坐标
        y_next = pos_next(2:2:end); % t+0.001s 时刻所有把手的 y 坐标

        % 计算速度 (vx, vy)，然后计算总速度
        v_x = (x_next - x_now) / delta_t;
        v_y = (y_next - y_now) / delta_t;

        % 计算每个把手的速度并存储
        velocity = sqrt(v_x.^2 + v_y.^2);

        % 检查是否有把手速度超过了2m/s
        exceeding_indices = find(velocity > max_v);  % 超过2m/s的把手索引
        if ~isempty(exceeding_indices)
            % 记录第一个超出2m/s的把手信息
            exceeding_speed = velocity(exceeding_indices(1));
            exceeding_section = exceeding_indices(1);
            exceeding_time = t_now;
            found_max_speed = true;
            break;
        end

        % 检查所有把手是否都在调头空间外
        in_space = sqrt(x_now.^2 + y_now.^2) - R_space <= 0.5;
        if all(~in_space)  % 如果所有把手都在调头空间外
            break;  % 结束当前速度的检测
        end
    end

    % 如果发现速度超过2m/s，跳出循环
    if found_max_speed
        v_head = v_head - v_step;  % 回退到上一个速度作为最大速度
        break;
    else
        v_head = v_head + v_step;  % 增加速度
    end
end

% 输出最终的最大速度
fprintf('龙头的最大行进速度为: %.2f m/s\n', v_head);
% 输出超出2m/s的具体速度、把手编号和时刻
if found_max_speed
    fprintf('超出2m/s的速度为: %.6f m/s\n', exceeding_speed);
    fprintf('超速的把手编号为: %d\n', exceeding_section);
    fprintf('超速发生的时间为: %.2f 秒\n', exceeding_time);
end


