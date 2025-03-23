clc; clear;

% 第一段代码中的螺旋线方程参数
a = 1.7;  % 螺距常数
R_space = 4.5;  % 调头空间半径

% 1. 盘入螺旋线的角度范围 (从 theta_in_start 到 theta_in_end)
theta_in_start = 16.6320;
theta_in_end = 32 * pi;  % 盘入螺旋线的结束角度
theta_in = linspace(theta_in_start, theta_in_end, 1000);  % 盘入螺旋线的角度范围
r_in = a * theta_in / (2 * pi);  % 盘入螺旋线的半径

% 转换为笛卡尔坐标
x_in = r_in .* cos(theta_in);
y_in = r_in .* sin(theta_in);

% 2. 盘出螺旋线 (与盘入螺旋线中心对称)
x_out = -x_in;
y_out = -y_in;

% 3. 调头空间的圆
theta_circle = linspace(0, 2 * pi, 500);  % 角度范围
x_circle = R_space * cos(theta_circle);  % 圆的 x 坐标
y_circle = R_space * sin(theta_circle);  % 圆的 y 坐标

% 优化后的 r1, r2, theta1, theta2 (假设的数值)
r1 = 2.1118;
r2 = 2.3963;
theta1 = 3.0215;
theta2 = 3.0215;

% 计算圆心坐标
% 盘入螺旋线的切线向量
theta_in_intersection = theta_in_start;  % 交点发生在盘入螺旋线的起始角度
r_in_intersection = a * theta_in_intersection / (2 * pi);
x_in_intersection = r_in_intersection * cos(theta_in_intersection);
y_in_intersection = r_in_intersection * sin(theta_in_intersection);
dr_dtheta_in = a / (2 * pi);
dx_dtheta_in = -(dr_dtheta_in * cos(theta_in_intersection) - r_in_intersection * sin(theta_in_intersection));
dy_dtheta_in = -(dr_dtheta_in * sin(theta_in_intersection) + r_in_intersection * cos(theta_in_intersection));
tangent_in_unit_vector = [dx_dtheta_in, dy_dtheta_in] / norm([dx_dtheta_in, dy_dtheta_in]);

% 顺时针旋转90度的切线向量
tangent_in_rotated = [tangent_in_unit_vector(2), -tangent_in_unit_vector(1)];

% 盘入螺旋线的圆心坐标
x_in_circle_center = x_in_intersection + r1 * tangent_in_rotated(1);
y_in_circle_center = y_in_intersection + r1 * tangent_in_rotated(2);

% 盘出螺旋线的切线向量 (与盘入螺旋线对称)
tangent_out_rotated = -tangent_in_rotated;

% 盘出螺旋线的圆心坐标
x_out_intersection = -x_in_intersection;
y_out_intersection = -y_in_intersection;
x_out_circle_center = x_out_intersection + r2 * tangent_out_rotated(1);
y_out_circle_center = y_out_intersection + r2 * tangent_out_rotated(2);

% 计算调头弧线
theta_arc1 = linspace(0, -theta1, 100);  % 第一段圆弧的角度范围
arc1_x = x_in_circle_center + r1 * cos(theta_arc1 + atan2(y_in_intersection - y_in_circle_center, x_in_intersection - x_in_circle_center));
arc1_y = y_in_circle_center + r1 * sin(theta_arc1 + atan2(y_in_intersection - y_in_circle_center, x_in_intersection - x_in_circle_center));

theta_arc2 = linspace(0, -theta2, 100);  % 第二段圆弧的角度范围
arc2_x = x_out_circle_center + r2 * cos(theta_arc2 + atan2(y_out_intersection - y_out_circle_center, x_out_intersection - x_out_circle_center));
arc2_y = y_out_circle_center + r2 * sin(theta_arc2 + atan2(y_out_intersection - y_out_circle_center, x_out_intersection - x_out_circle_center));

% ---------------- 动态仿真设置 ------------------

% 读取 Excel 文件中的位置数据
data = readmatrix('result4.xlsx', 'Sheet', '位置');

% 数据的第一行是列头，第一列是行头，所以需要去掉
positions_data = data(1:end, 2:end);

% 将数据转换为矩阵
positions = reshape(positions_data', [], size(positions_data, 1))';

% 设置时间段参数
start_time = -100;   % 开始时间 (秒)
end_time = 100;     % 结束时间 (秒)
time_steps = start_time:1:end_time;  % 时间数组，对应201个时刻
n_sections = 224;  % 总板凳数量

% ---------------- 创建图像和动画 ------------------

figure_handle = figure;  % 获取图像句柄
hold on;
grid on;
axis equal;
xlim([-12, 12]);  % 调整显示区域
ylim([-12, 12]);  % 调整显示区域
title('板凳龙位置仿真与螺旋线图示');
xlabel('X Position (m)');
ylabel('Y Position (m)');

% 绘制静态的螺旋线轨迹和调头空间
plot(x_in, y_in, 'b', 'LineWidth', 1.5);  % 盘入螺旋线
plot(x_out, y_out, 'r', 'LineWidth', 1.5);  % 盘出螺旋线
plot(x_circle, y_circle, 'k--', 'LineWidth', 1.5);  % 调头空间的圆

% 绘制调头弧线
plot(arc1_x, arc1_y, 'g', 'LineWidth', 2);  % 第一段调头弧线
plot(arc2_x, arc2_y, 'g', 'LineWidth', 2);  % 第二段调头弧线

% 动态绘制板凳龙的位置
for t_idx = 100:201  % 时间范围对应的位置索引，从第101个时间步到第201个
    % 如果图像窗口已关闭，则中止动画
    if ~isvalid(figure_handle)
        disp('窗口已关闭，动画中止');
        break;
    end
    
    % 清除之前的板凳绘图（不清除螺旋线、调头空间和调头弧线）
    cla;
    plot(x_in, y_in, 'b', 'LineWidth', 1.5);  % 重新绘制螺旋轨迹
    plot(x_out, y_out, 'r', 'LineWidth', 1.5);  % 重新绘制盘出螺旋线
    plot(x_circle, y_circle, 'k--', 'LineWidth', 1.5);  % 重新绘制调头圆
    plot(arc1_x, arc1_y, 'g', 'LineWidth', 2);  % 重新绘制第一段调头弧线
    plot(arc2_x, arc2_y, 'g', 'LineWidth', 2);  % 重新绘制第二段调头弧线

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
    pause(0.1);  % 调整 pause 值可以改变动画速度
end

% Hold the final plot
hold off;
