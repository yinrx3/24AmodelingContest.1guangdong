clc; clear;

% 螺旋线方程参数
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

% 优化后的 r1, r2, theta1, theta2 (替换为优化结果)
r1 = 2.1118;
r2 = 2.3963;
theta1 = 3.0215;
theta2 = 3.0215;

% 4. 计算圆心坐标
% 盘入螺旋线切线向量
dr_dtheta_in = a / (2 * pi);
theta_in_intersection = theta_in_start;  % 交点发生在盘入螺旋线的起始角度
r_in_intersection = a * theta_in_intersection / (2 * pi);
x_in_intersection = r_in_intersection * cos(theta_in_intersection);
y_in_intersection = r_in_intersection * sin(theta_in_intersection);
dx_dtheta_in = -(dr_dtheta_in * cos(theta_in_intersection) - r_in_intersection * sin(theta_in_intersection));
dy_dtheta_in = -(dr_dtheta_in * sin(theta_in_intersection) + r_in_intersection * cos(theta_in_intersection));
tangent_in_unit_vector = [dx_dtheta_in, dy_dtheta_in] / norm([dx_dtheta_in, dy_dtheta_in]);

% 顺时针旋转90度的切线向量
tangent_in_rotated = [tangent_in_unit_vector(2), -tangent_in_unit_vector(1)];

% 盘入螺旋线的圆心坐标
x_in_circle_center = x_in_intersection + r1 * tangent_in_rotated(1)
y_in_circle_center = y_in_intersection + r1 * tangent_in_rotated(2)

% 盘出螺旋线的切线向量 (与盘入螺旋线对称)
tangent_out_rotated = -tangent_in_rotated;

% 盘出螺旋线的圆心坐标
x_out_intersection = -x_in_intersection;
y_out_intersection = -y_in_intersection;
x_out_circle_center = x_out_intersection + r2 * tangent_out_rotated(1)
y_out_circle_center = y_out_intersection + r2 * tangent_out_rotated(2)

% 5. 计算两段圆弧
% 第一段圆弧从盘入螺旋线的圆心出发，沿盘入方向顺时针旋转 theta1
% 第二段圆弧从盘出螺旋线的圆心出发，沿盘出方向顺时针旋转 theta2

% 调头弧线1: 从圆心到螺旋线交点的向量，顺时针旋转 theta1
theta_arc1 = linspace(0, -theta1, 100);  % 顺时针旋转 theta1
arc1_x = x_in_circle_center + r1 * cos(theta_arc1 + atan2(y_in_intersection - y_in_circle_center, x_in_intersection - x_in_circle_center));
arc1_y = y_in_circle_center + r1 * sin(theta_arc1 + atan2(y_in_intersection - y_in_circle_center, x_in_intersection - x_in_circle_center));

% 调头弧线2: 从盘出螺旋线的圆心出发，顺时针旋转 theta2
theta_arc2 = linspace(0, -theta2, 100);  % 顺时针旋转 theta2
arc2_x = x_out_circle_center + r2 * cos(theta_arc2 + atan2(y_out_intersection - y_out_circle_center, x_out_intersection - x_out_circle_center));
arc2_y = y_out_circle_center + r2 * sin(theta_arc2 + atan2(y_out_intersection - y_out_circle_center, x_out_intersection - x_out_circle_center));

% 6. 画图
figure;
hold on;
axis equal;
grid on;

% 盘入螺旋线
plot(x_in, y_in, 'b', 'LineWidth', 1.5);

% 盘出螺旋线
plot(x_out, y_out, 'r', 'LineWidth', 1.5);

% 调头空间的圆
plot(x_circle, y_circle, 'k--', 'LineWidth', 1.5);

% 圆心
plot(x_in_circle_center, y_in_circle_center, 'bo', 'MarkerFaceColor', 'b');
plot(x_out_circle_center, y_out_circle_center, 'ro', 'MarkerFaceColor', 'r');

% 调头弧线第一段
plot(arc1_x, arc1_y, 'g', 'LineWidth', 2);

% 调头弧线第二段
plot(arc2_x, arc2_y, 'g', 'LineWidth', 2);

title('螺旋线与调头空间的图示');
legend('盘入螺旋线', '盘出螺旋线', '调头空间', '盘入圆心', '盘出圆心', '调头弧线');
hold off;
