clc; clear;

% 螺旋线方程参数
a = 1.7;  % 螺距常数
R_space = 4.5;  % 调头空间半径

% 1. 盘入螺旋线交点 (极坐标 -> 笛卡尔坐标)
theta_in_intersection = (R_space * 2 * pi) / a;  % 盘入交点角度
r_in_intersection = a * theta_in_intersection / (2 * pi);  % 盘入交点极径
x_in_intersection = r_in_intersection * cos(theta_in_intersection); % 转换为 x 坐标
y_in_intersection = r_in_intersection * sin(theta_in_intersection); % 转换为 y 坐标

% 2. 盘出螺旋线交点 (极坐标 -> 笛卡尔坐标)
x_out_intersection = -x_in_intersection;
y_out_intersection = -y_in_intersection;

% 3. 计算盘入螺旋线的切线向量并旋转 90 度 (顺时针旋转)
dr_dtheta_in = a / (2 * pi);
dx_dtheta_in = -(dr_dtheta_in * cos(theta_in_intersection) - r_in_intersection * sin(theta_in_intersection));
dy_dtheta_in = -(dr_dtheta_in * sin(theta_in_intersection) + r_in_intersection * cos(theta_in_intersection));
tangent_in_unit_vector = [dx_dtheta_in, dy_dtheta_in] / norm([dx_dtheta_in, dy_dtheta_in]);  % 单位向量

% 顺时针旋转90度
tangent_in_rotated = [tangent_in_unit_vector(2), -tangent_in_unit_vector(1)];

% 4. 盘出螺旋线的切线向量 (与盘入螺旋线对称)
tangent_out_rotated = -tangent_in_rotated;

% 优化后的 r1, r2, theta1, theta2 (替换为优化结果)
r1 = 2.1118;
r2 = 2.3963;
theta1 = 3.0215;
theta2 = 3.0215;
s = theta1 * r1 + theta2 * r2;
fprintf("最短调头路径是为%0.6f \n",s);

% 5. 盘入螺旋线的圆心坐标
x_in_circle_center = x_in_intersection + r1 * tangent_in_rotated(1);
y_in_circle_center = y_in_intersection + r1 * tangent_in_rotated(2);

% 6. 盘出螺旋线的圆心坐标
x_out_circle_center = x_out_intersection + r2 * tangent_out_rotated(1);
y_out_circle_center = y_out_intersection + r2 * tangent_out_rotated(2);

% 7. 从圆心到交点的向量计算
vector_in = [x_in_intersection - x_in_circle_center, y_in_intersection - y_in_circle_center];
vector_out = [x_out_intersection - x_out_circle_center, y_out_intersection - y_out_circle_center];

% 8. 顺时针旋转后的盘入向量
vector_in_rotated = [vector_in(1) * cos(theta1) + vector_in(2) * sin(theta1), ...
                     -vector_in(1) * sin(theta1) + vector_in(2) * cos(theta1)];

% 9. 顺时针旋转后的盘出向量
vector_out_rotated = [vector_out(1) * cos(theta2) + vector_out(2) * sin(theta2), ...
                      -vector_out(1) * sin(theta2) + vector_out(2) * cos(theta2)];

% 10. 单位化两个旋转后的向量
vector_in_rotated = vector_in_rotated / norm(vector_in_rotated);
vector_out_rotated = vector_out_rotated / norm(vector_out_rotated);

% 11. 检查方向相反的约束
dot_product = dot(vector_in_rotated, vector_out_rotated);
fprintf('方向约束：内积 = %.6f (应接近 -1)\n', dot_product);

% 12. 检查终点坐标相等的约束
x_in_rotated_end = x_in_circle_center + r1 * vector_in_rotated(1);
y_in_rotated_end = y_in_circle_center + r1 * vector_in_rotated(2);

x_out_rotated_end = x_out_circle_center + r2 * vector_out_rotated(1);
y_out_rotated_end = y_out_circle_center + r2 * vector_out_rotated(2);

% 13. 终点差距
x_diff = abs(x_in_rotated_end - x_out_rotated_end);
y_diff = abs(y_in_rotated_end - y_out_rotated_end);
fprintf('终点约束：x 坐标差 = %.6f, y 坐标差 = %.6f (应接近 0)\n', x_diff, y_diff);
