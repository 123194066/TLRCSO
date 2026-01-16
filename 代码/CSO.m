function [Best_Score, Convergence_curve, BestFit, f1_zong, f2_zong, f3_zong, f4_zong] = CSO(SearchAgents_no, dim, Max_iteration, ub, lb, pos_tag, all_param)
% 鸡群优化算法（修正版）：修复Link_Table参数不足问题

Tma = Max_iteration - 1;  % 修正变量名，与输入参数一致
G = 2;
rPercent = 0.1;
hPercent = 0.8;
mPercent = 0.5;
rNum = round(SearchAgents_no * rPercent);    % 公鸡数量
hNum = round(SearchAgents_no * hPercent);    % 母鸡数量
cNum = SearchAgents_no - rNum - hNum;        % 小鸡数量
mNum = round(hNum * mPercent);               % 母母鸡数量

% 转换上下界为行向量
ub_row = ub';
lb_row = lb';

% 初始化种群
TOP = initialization(SearchAgents_no, dim);
x = TOP';  % 一行一个个体
 f1_x = zeros(1, SearchAgents_no);
    f2_x = zeros(1, SearchAgents_no);
    f3_x = zeros(1, SearchAgents_no);
    f4_x = zeros(1, SearchAgents_no);
    fit = zeros(1, SearchAgents_no);
% 初始化适应度（修正Link_Table调用，添加all_param参数）
for i = 1 : SearchAgents_no
    PP = reshape(x(i,:), 2, dim/2)';
    % 关键修复：传入第三个参数all_param，并用~接收多余的返回值
    [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
    fit(i) = W;
    f1_x(i) = f1;
    f2_x(i) = f2;
    f3_x(i) = f3;
    f4_x(i) = f4;
end

pFit = fit;  % 个体最优适应度
f1Fit = f1_x;
f2Fit = f2_x;
f3Fit = f3_x;
f4Fit = f4_x;
pX = x;      % 个体最优位置
[ fMin, bestIndex ] = min( fit );  % 全局最优
bestX = x(bestIndex, :);

% 初始化收敛曲线
Convergence_curve = zeros(1, Max_iteration);
f1_zong = zeros(1, Max_iteration);
f2_zong = zeros(1, Max_iteration);
f3_zong = zeros(1, Max_iteration);
f4_zong = zeros(1, Max_iteration);
Convergence_curve(1) = fMin;
f1_zong(1) = f1Fit(bestIndex);
f2_zong(1) = f2Fit(bestIndex);
f3_zong(1) = f3Fit(bestIndex);
f4_zong(1) = f4Fit(bestIndex);

% 主循环
for t = 1 : Tma
    disp(['CSO迭代次数 = ', num2str(t)]);
    FL = rand(SearchAgents_no, 1) .* 0.4 + 0.5;  % 小鸡跟随系数
    
    % 每G代更新群体结构
    if mod(t, G) == 1
        [~, sortIndex] = sort(pFit);
        motherLib = randpermD(hNum, mNum) + rNum;
        mate = randpermF(rNum, hNum);
        mother = motherLib(randi(mNum, cNum, 1));
    end
    
    % 公鸡位置更新（修正Link_Table调用）
    for i = 1 : rNum
        anotherRooster = randiTabu(1, rNum, i, 1);
        if pFit(sortIndex(i)) <= pFit(sortIndex(anotherRooster))
            tempSigma = 1;
        else
            tempSigma = exp((pFit(sortIndex(anotherRooster)) - pFit(sortIndex(i))) / ...
                           (abs(pFit(sortIndex(i))) + realmin));
        end
        
        x(sortIndex(i), :) = pX(sortIndex(i), :) .* (1 + tempSigma .* randn(1, dim));
        x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), lb_row, ub_row);
        
        % 关键修复：传入all_param参数
        PP = reshape(x(sortIndex(i), :), 2, dim/2)';
        [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
        fit(sortIndex(i)) = W;
        f1_x(sortIndex(i)) = f1;
        f2_x(sortIndex(i)) = f2;
        f3_x(sortIndex(i)) = f3;
        f4_x(sortIndex(i)) = f4;
    end
    
    % 母鸡位置更新（修正Link_Table调用）
    for i = (rNum + 1) : (rNum + hNum)
        other = randiTabu(1, i, mate(i - rNum), 1);
        c1 = exp((pFit(sortIndex(i)) - pFit(sortIndex(mate(i - rNum)))) / ...
                (abs(pFit(sortIndex(i))) + realmin));
        c2 = exp(-pFit(sortIndex(i)) + pFit(sortIndex(other)));
        c2 = min(c2, 2);
        
        x(sortIndex(i), :) = pX(sortIndex(i), :) + ...
            (pX(sortIndex(mate(i - rNum)), :) - pX(sortIndex(i), :)) .* c1 .* rand(1, dim) + ...
            (pX(sortIndex(other), :) - pX(sortIndex(i), :)) .* c2 .* rand(1, dim);
        x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), lb_row, ub_row);
        
        % 关键修复：传入all_param参数
        PP = reshape(x(sortIndex(i), :), 2, dim/2)';
        [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
        fit(sortIndex(i)) = W;
        f1_x(sortIndex(i)) = f1;
        f2_x(sortIndex(i)) = f2;
        f3_x(sortIndex(i)) = f3;
        f4_x(sortIndex(i)) = f4;
    end
    
    % 小鸡位置更新（修正Link_Table调用）
    for i = (rNum + hNum + 1) : SearchAgents_no
        x(sortIndex(i), :) = pX(sortIndex(i), :) + ...
            (pX(sortIndex(mother(i - rNum - hNum)), :) - pX(sortIndex(i), :)) .* FL(i);
        x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), lb_row, ub_row);
        
        % 关键修复：传入all_param参数
        PP = reshape(x(sortIndex(i), :), 2, dim/2)';
        [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
        fit(sortIndex(i)) = W;
        f1_x(sortIndex(i)) = f1;
        f2_x(sortIndex(i)) = f2;
        f3_x(sortIndex(i)) = f3;
        f4_x(sortIndex(i)) = f4;
    end
    
    % 更新最优解
    for i = 1 : SearchAgents_no
        if fit(i) < pFit(i)
            pFit(i) = fit(i);
            pX(i, :) = x(i, :);
            f1Fit(i) = f1_x(i);
            f2Fit(i) = f2_x(i);
            f3Fit(i) = f3_x(i);
            f4Fit(i) = f4_x(i);
        end
        
        if pFit(i) < fMin
            fMin = pFit(i);
            bestX = pX(i, :);
        end
    end
    
    % 记录收敛曲线
    Convergence_curve(t + 1) = fMin;
    [~, bestIndex] = find(pFit == fMin);
    f1_zong(t + 1) = f1Fit(min(bestIndex));
    f2_zong(t + 1) = f2Fit(min(bestIndex));
    f3_zong(t + 1) = f3Fit(min(bestIndex));
    f4_zong(t + 1) = f4Fit(min(bestIndex));
end

Best_Score = fMin;
Convergence_curve = Convergence_curve;
BestFit = bestX;
end

