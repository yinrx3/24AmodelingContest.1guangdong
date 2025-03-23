clc; clear;

% 螺旋线方程参数
a = 1.7;  % 螺距常数
R_space = 4.5;  % 调头空间半径

% 1. 盘入螺旋线交点 (极坐标 -> 笛卡尔坐标)
theta_in_intersection = (R_space * 2 * pi) / a;  % 盘入交点角度
r_in_intersection = a * theta_in_intersection / (2 * pi);  % 盘入交点极径
x_in_intersection = r_in_intersection * cos(theta_in_intersection);  % 转换为 x 坐标
y_in_intersection = r_in_intersection * sin(theta_in_intersection);  % 转换为 y 坐标

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

% 优化目标函数
objective = @(x) x(1) * x(3) + x(2) * x(4);  % r1*theta1 + r2*theta2

% 约束函数定义
function [c, ceq] = constraints(x, x_in_intersection, y_in_intersection, x_out_intersection, y_out_intersection, ...
                                tangent_in_rotated, tangent_out_rotated)

    r1 = x(1); r2 = x(2); theta1 = x(3); theta2 = x(4);

    % 盘入螺旋线的圆心坐标
    x_in_circle_center = x_in_intersection + r1 * tangent_in_rotated(1);
    y_in_circle_center = y_in_intersection + r1 * tangent_in_rotated(2);

    % 盘出螺旋线的圆心坐标
    x_out_circle_center = x_out_intersection + r2 * tangent_out_rotated(1);
    y_out_circle_center = y_out_intersection + r2 * tangent_out_rotated(2);

    % 从圆心到交点的向量计算
    vector_in = [x_in_intersection - x_in_circle_center, y_in_intersection - y_in_circle_center];
    vector_out = [x_out_intersection - x_out_circle_center, y_out_intersection - y_out_circle_center];

    % 顺时针旋转后的盘入向量
    vector_in_rotated = [vector_in(1) * cos(theta1) + vector_in(2) * sin(theta1), ...
                         -vector_in(1) * sin(theta1) + vector_in(2) * cos(theta1)];

    % 顺时针旋转后的盘出向量
    vector_out_rotated = [vector_out(1) * cos(theta2) + vector_out(2) * sin(theta2), ...
                          -vector_out(1) * sin(theta2) + vector_out(2) * cos(theta2)];

    % 单位化两个旋转后的向量
    vector_in_rotated = vector_in_rotated / norm(vector_in_rotated);
    vector_out_rotated = vector_out_rotated / norm(vector_out_rotated);

    % 约束1：两个向量方向相反
    c = dot(vector_in_rotated, vector_out_rotated) + 1;

    % 约束2：两个终点坐标相等
    x_in_rotated_end = x_in_circle_center + r1 * vector_in_rotated(1);
    y_in_rotated_end = y_in_circle_center + r1 * vector_in_rotated(2);

    x_out_rotated_end = x_out_circle_center + r2 * vector_out_rotated(1);
    y_out_rotated_end = y_out_circle_center + r2 * vector_out_rotated(2);

    ceq = [x_in_rotated_end - x_out_rotated_end, y_in_rotated_end - y_out_rotated_end];
end

% 初始猜测值 [r1, r2, theta1, theta2]
x0 = [1, 1, pi/4, pi/4];

% 优化求解选项
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1000);

% 定义边界，确保r1和r2为正数
lb = [0, 0, -Inf, -Inf];  % r1 和 r2 的下界
ub = [];  % 不设置上界

% 定义优化问题的约束函数
constr_fun = @(x) constraints(x, x_in_intersection, y_in_intersection, x_out_intersection, y_out_intersection, ...
                              tangent_in_rotated, tangent_out_rotated);

% 使用fmincon进行优化求解
x_opt = fmincon(objective, x0, [], [], [], [], lb, ub, constr_fun, options);

% 输出结果
fprintf('优化后的 r1: %.4f, r2: %.4f, theta1: %.4f, theta2: %.4f\n', x_opt(1), x_opt(2), x_opt(3), x_opt(4));
