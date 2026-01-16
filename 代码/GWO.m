function [best_score, convergence_curve, best_pos, f1_hist, f2_hist, f3_hist, f4_hist] = GWO(Size, dimension, Tmax, ub, lb, pos_tag, all_param)


r_num = dimension/2;
Tma = Tmax-1;


xma = ub(:)'; % 强制为 1 x dim 行向量
xmi = lb(:)'; % 强制为 1 x dim 行向量

% 初始化种群 (狼群)
x = xmi + rand(Size, dimension) .* (xma - xmi);

% 计算初始适应度
fit = zeros(1, Size);
f1_x = zeros(1, Size);
f2_x = zeros(1, Size);
f3_x = zeros(1, Size);
f4_x = zeros(1, Size);

for i = 1:Size
    PP = reshape(x(i,:), 2, r_num)';
    [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
    fit(i) = W;
    f1_x(i) = f1;
    f2_x(i) = f2;
    f3_x(i) = f3;
    f4_x(i) = f4;
end

% 初始化 Alpha, Beta, Delta 狼
[~, sorted_indices] = sort(fit);
Alpha_pos = x(sorted_indices(1), :);
Alpha_score = fit(sorted_indices(1));

Beta_pos = x(sorted_indices(2), :);
Beta_score = fit(sorted_indices(2));

Delta_pos = x(sorted_indices(3), :);
Delta_score = fit(sorted_indices(3));

% 初始化历史记录
convergence_curve = zeros(1, Tma + 1);
f1_hist = zeros(1, Tma + 1);
f2_hist = zeros(1, Tma + 1);
f3_hist = zeros(1, Tma + 1);
f4_hist = zeros(1, Tma + 1);

convergence_curve(1) = Alpha_score;
f1_hist(1) = f1_x(sorted_indices(1));
f2_hist(1) = f2_x(sorted_indices(1));
f3_hist(1) = f3_x(sorted_indices(1));
f4_hist(1) = f4_x(sorted_indices(1));

% --- 主循环 ---
for t = 1 : Tma
    disp(['GWO 迭代: ',num2str(t+1)]);
    
    a = 2 - t * (2 / Tma); % a 从 2 线性递减到 0
    
    for i = 1:Size
        for j = 1:dimension
            r1 = rand(); r2 = rand();
            A1 = 2*a*r1 - a;
            C1 = 2*r2;
            D_alpha = abs(C1*Alpha_pos(j) - x(i,j));
            X1 = Alpha_pos(j) - A1*D_alpha;
            
            r1 = rand(); r2 = rand();
            A2 = 2*a*r1 - a;
            C2 = 2*r2;
            D_beta = abs(C2*Beta_pos(j) - x(i,j));
            X2 = Beta_pos(j) - A2*D_beta;
            
            r1 = rand(); r2 = rand();
            A3 = 2*a*r1 - a;
            C3 = 2*r2;
            D_delta = abs(C3*Delta_pos(j) - x(i,j));
            X3 = Delta_pos(j) - A3*D_delta;
            
            x(i,j) = (X1 + X2 + X3)/3;
        end
        % 边界检查
        x(i, :) = Bounds(x(i, :), xmi, xma);
    end
    
    % 更新所有狼的适应度
    for j = 1:Size
        PP = reshape(x(j,:), 2, r_num)';
        [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
        
        % 如果新的解更好，则更新
        if W < fit(j)
            fit(j) = W;
            f1_x(j) = f1;
            f2_x(j) = f2;
            f3_x(j) = f3;
            f4_x(j) = f4;
        end
    end
    
    % 更新 Alpha, Beta, Delta
    [~, sorted_indices] = sort(fit);
    
    % 如果有更好的解，则更新全局最优
    if fit(sorted_indices(1)) < Alpha_score
        Alpha_score = fit(sorted_indices(1));
        Alpha_pos = x(sorted_indices(1),:);
    end
    
    Beta_score = fit(sorted_indices(2));
    Beta_pos = x(sorted_indices(2),:);
    
    Delta_score = fit(sorted_indices(3));
    Delta_pos = x(sorted_indices(3),:);

    % 记录历史
    convergence_curve(t+1) = Alpha_score;
    best_idx = find(fit == Alpha_score, 1, 'first');
    f1_hist(t+1) = f1_x(best_idx);
    f2_hist(t+1) = f2_x(best_idx);
    f3_hist(t+1) = f3_x(best_idx);
    f4_hist(t+1) = f4_x(best_idx);
end

best_score = Alpha_score;
best_pos = Alpha_pos;
end


function s = Bounds(s, Lb, Ub)
    temp = s;
    I = temp < Lb;
    temp(I) = Lb(I) + rand(1,sum(I)).*(Ub(I)-Lb(I)); % 越界后随机生成
    
    J = temp > Ub;
    temp(J) = Lb(J) + rand(1,sum(J)).*(Ub(J)-Lb(J)); % 越界后随机生成
    
    s = temp;
end


