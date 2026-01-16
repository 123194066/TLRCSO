function [Best_Score, Convergence_curve, BestFit, f1_zong, f2_zong, f3_zong, f4_zong, BestIndex] = EMPCSO23(SearchAgents_no, dim, Max_iteration, ub, lb, pos_tag, all_param)

%levy
r_num = dim / 2;
Tma = Max_iteration - 1;
G = 2;
hPercent = 0.8;
rPercent = 0.1;
NWmax = 2;
NWmin = 0.5;
mPercent = 0.5;
% pd_normal = makedist('Normal', 'mu', 0, 'sigma', 1); % 标准高斯对象

rNum = round(SearchAgents_no * rPercent);
hNum = round(SearchAgents_no * hPercent);
cNum = SearchAgents_no - rNum - hNum;
mNum = round(hNum * mPercent);

xma = ub';
xmi = lb';

TOP = initialization(SearchAgents_no, dim);
x = TOP';

f1_x = zeros(1, SearchAgents_no);
f2_x = zeros(1, SearchAgents_no);
f3_x = zeros(1, SearchAgents_no);
f4_x = zeros(1, SearchAgents_no);
fit = zeros(1, SearchAgents_no);
contribution = cell(1, SearchAgents_no);  % 存储每个个体的贡献度向量

f1Fit = zeros(1, SearchAgents_no);
f2Fit = zeros(1, SearchAgents_no);
f3Fit = zeros(1, SearchAgents_no);
f4Fit = zeros(1, SearchAgents_no);

for i = 1:SearchAgents_no
    PP = reshape(x(i,:), 2, r_num)';
    [W, f1, f2, f3, f4, temp_alpha_m] = Link_Table(PP, pos_tag, all_param);
    fit(i) = W;
    f1_x(i) = f1;
    f2_x(i) = f2;
    f3_x(i) = f3;
    f4_x(i) = f4;
    contribution{i} = temp_alpha_m;  % 保存贡献度向量
    f1Fit(i) = f1;
    f2Fit(i) = f2;
    f3Fit(i) = f3;
    f4Fit(i) = f4;
end

pFit = fit;  % 个体历史最优适应度
pX = x;      % 个体历史最优位置
  pContribution = contribution;  % 个体历史最优对应的贡献度

[fMin, bestIndex] = min(pFit);
bestX = x(bestIndex, :);
Convergence_curve(1) = fMin;


f1_zong(1) = f1Fit(bestIndex);
f2_zong(1) = f2Fit(bestIndex);
f3_zong(1) = f3Fit(bestIndex);
f4_zong(1) = f4Fit(bestIndex);
window_size = 20;  
conv_threshold = Convergence_curve(1) * 0.001; 
last_scores = zeros(1, window_size);  
last_scores(1) = Convergence_curve(1);  
is_fine_mode = false;  
for t = 1:Tma
      idx = mod(t-1, window_size) + 1;  
    last_scores(idx) = Convergence_curve(t);  

    if t >= window_size
        conv_speed = (max(last_scores) - min(last_scores)) / window_size;
   
        if conv_speed < conv_threshold && ~is_fine_mode
            is_fine_mode = true;
            disp(['Iteration ', num2str(t), ': 收敛速度缓慢，切换搜索模式']);
        end
    end
       disp(['EMPCSO23 iteration = ', num2str(t), '   ', 'BEST Fitness = ', num2str(fMin)]);
    FL = rand(SearchAgents_no, 1) .* 0.4 + 0.5;

    if mod(t, G) == 1
        [~, sortIndex] = sort(pFit);
        x = pX;

        motherLib = randpermD(hNum, mNum) + rNum;

        for i = 1:SearchAgents_no
            Xx(i,:) = x(sortIndex(i), :);
        end

        oudist = pdist2(Xx, Xx);
        hdist = oudist(rNum+1:end, 1:rNum);
        [~, h_sortIndex] = sort(hdist');
        mate = h_sortIndex(1:3, 1:hNum)';
%         mother = motherLib(randi(mNum, cNum, 1));
        mother= randi(rNum, cNum, 1);  % 公鸡索引
    end

    % 公鸡更新
    for i = 1:rNum
        NW = NWmax - (NWmax - NWmin) * t / Tma;
        anotherRooster = randiTabu(1, rNum, i, 1);

        if pFit(sortIndex(i)) <= pFit(sortIndex(anotherRooster))
            tempSigma = NW .* [0.5 0.5];
        else
            tempSigma = NW .* exp( (pFit(sortIndex(anotherRooster)) - ...
                pFit(sortIndex(i))) / (abs(pFit(sortIndex(i))) + realmin) ) .* [0.5 0.5];
        end

        PP = reshape(pX(sortIndex(i), :), 2, r_num)';
        [W, f1, f2, f3, f4, temp_alpha_m] = Link_Table(PP, pos_tag, all_param);
        fit(sortIndex(i)) = W;
        f1_x(sortIndex(i)) = f1;
        f2_x(sortIndex(i)) = f2;
        f3_x(sortIndex(i)) = f3;
        f4_x(sortIndex(i)) = f4;
 pContribution{sortIndex(i)} = temp_alpha_m;  % 更新贡献度
       zero_contribution_idx = find(temp_alpha_m == 0);
            if ~isempty(zero_contribution_idx)
                % 若存在贡献度为0的单体，优先从这些单体中选择扰动对象
                r_index = zero_contribution_idx(randi(length(zero_contribution_idx)));
            else
                % 若不存在，保持原随机选择逻辑
                r_index = randpermD(r_num, 1);
            end
        xx = zeros(5, dim);
        refit = zeros(1, 5);
        f1_buf = zeros(1, 5); f2_buf = zeros(1, 5); f3_buf = zeros(1, 5); f4_buf = zeros(1, 5);

        for iii = 1:5
            x(sortIndex(i), (r_index-1)*2+1:r_index*2) = ...
                pX(sortIndex(i), (r_index-1)*2+1:r_index*2) .* (1 + tempSigma .* randn(1, 2));
            x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), xmi, xma);
            xx(iii, :) = x(sortIndex(i), :);

            PP = reshape(x(sortIndex(i), :), 2, r_num)';
            [W, f1, f2, f3, f4, ~,] = Link_Table(PP, pos_tag, all_param);

            refit(iii) = W;
            f1_buf(iii) = f1; f2_buf(iii) = f2;
            f3_buf(iii) = f3; f4_buf(iii) = f4;
        end

        [min_re, sort_re] = min(refit);
        x(sortIndex(i), :) = xx(sort_re, :);
        fit(sortIndex(i)) = min_re;
        f1_x(sortIndex(i)) = f1_buf(sort_re);
        f2_x(sortIndex(i)) = f2_buf(sort_re);
        f3_x(sortIndex(i)) = f3_buf(sort_re);
        f4_x(sortIndex(i)) = f4_buf(sort_re);
    end
    % 母鸡更新

aa = 2 * exp(-t / 120) - 0.003 * t+ 1.4;
beta = 1.0 + (t / Tma) * 0.8; 
    sigma_u = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

    for i = (rNum + 1):(rNum + hNum)
         if ~is_fine_mode  % % 正常灰狼更新
            BB = mate(i - rNum, 1:3);
            [CC, ~] = sort(BB);
            
            f11=pFit(sortIndex(CC(1)));
            f22=pFit(sortIndex(CC(2)));
            f33=pFit(sortIndex(CC(3)));       
            
            w1=1/(1+exp(-f11*2));% 放大最优个体的权重
            w2=1/(1+exp(-f22));
            w3=1/(1+exp(-f33));
            w=w1+w2+w3;
            w1=w1/w;  w2=w2/w;  w3=w3/w;
            
            u1 = sigma_u * randn(1, dim);
            v1 = randn(1, dim);
            L_steps1 = (u1 ./ abs(v1).^(1/beta));  
            
            u2 = sigma_u * randn(1, dim);
            v2 = randn(1, dim);
            L_steps2 =  (u2 ./ abs(v2).^(1/beta));  
            
            u3 = sigma_u * randn(1, dim);
            v3 = randn(1, dim);
            L_steps3 =  (u3 ./ abs(v3).^(1/beta));  
            for j = 1:dim
               L_step1 = L_steps1(j);       r11 = 2*rand();
               D1 = pX(sortIndex(CC(1)), j) - x(sortIndex(i), j);
               X1 = pX(sortIndex(CC(1)), j) - aa * L_step1 * D1* r11;
       
               L_step2 = L_steps2(j); r22 = 2*rand();
               D2 = pX(sortIndex(CC(2)), j) - x(sortIndex(i), j);
               X2 = pX(sortIndex(CC(2)), j) - aa * L_step2 * D2*r22;
       
               L_step3 = L_steps3(j); r33 = 2*rand();
               D3 = pX(sortIndex(CC(3)), j) - x(sortIndex(i), j);
               X3 = pX(sortIndex(CC(3)), j) - aa * L_step3 * D3*r33;
       
               x(sortIndex(i), j) = w1*X1 + w2*X2 +w3* X3;
            end

        else % 只更新2到5个阅读器，且只围绕距离最近的公鸡
            BB = mate(i - rNum, 1:3);
            [CC, ~] = sort(BB);
            
            % f11 = pFit(sortIndex(CC(1)));
            u1 = sigma_u * randn(1, dim);
            v1 = randn(1, dim);
            L_steps1 = (u1 ./ abs(v1).^(1/beta)); 

    update_count = randi([2,5]);  % 每次更新2-5个阅读器
    random_index = randi([1, dim/2], [1, update_count]);  % 随机选阅读器索引
    hen_update_index = zeros(1, dim);
    hen_update_index([random_index.*2-1; random_index.*2]) = 1;  % 标记更新维度

    if rand < 0.08
        random_disturb_index = randi(dim);  % 随机选1个维度扰动
        hen_update_index(random_disturb_index) = 1;
    end

    for j = 1:dim
        if hen_update_index(j) == 0; continue; end
        L_step1 = L_steps1(j);      
        D1 = pX(sortIndex(CC(1)), j) - x(sortIndex(i), j);
        X1 = 0.8 * (pX(sortIndex(CC(1)), j) - aa * L_step1 * D1) + 0.2 * bestX(j);  
        x(sortIndex(i), j) = X1;
    end            

        end

        x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), xmi, xma); % 检查边界
        PP = reshape(x(sortIndex(i), :), 2, r_num)';
        [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
        fit(sortIndex(i)) = W;
        f1_x(sortIndex(i)) = f1;
        f2_x(sortIndex(i)) = f2;
        f3_x(sortIndex(i)) = f3;
        f4_x(sortIndex(i)) = f4;
    end

    % 小鸡更新
    for i = (rNum + hNum + 1):SearchAgents_no
        x(sortIndex(i), :) = pX(sortIndex(i), :) + ...
            (pX(sortIndex(mother(i - rNum - hNum)), :) - pX(sortIndex(i), :)) .* FL(i);
        x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), xmi, xma);
        PP = reshape(x(sortIndex(i), :), 2, r_num)';
        [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
        fit(sortIndex(i)) = W;
        f1_x(sortIndex(i)) = f1;
        f2_x(sortIndex(i)) = f2;
        f3_x(sortIndex(i)) = f3;
        f4_x(sortIndex(i)) = f4;
    end

  
    for i = 1:SearchAgents_no
   
        if fit(i) < pFit(i)
            pFit(i) = fit(i);
            pX(i, :) = x(i, :);
          
            f1Fit(i) = f1_x(i);
            f2Fit(i) = f2_x(i);
            f3Fit(i) = f3_x(i);
            f4Fit(i) = f4_x(i);
             % 更新历史最优对应的贡献度
            PP = reshape(x(i, :), 2, r_num)';
            [~, ~, ~, ~, ~,  temp_alpha_m] = Link_Table(PP, pos_tag, all_param);
            pContribution{i} = temp_alpha_m;
        end
        
        if pFit(i) < fMin
            fMin = pFit(i);
            bestX = pX(i, :);
        end
    end

    [~, bestIndex] = min(pFit);
    Convergence_curve(1 + t) = fMin;
    f1_zong(t + 1) = f1Fit(bestIndex);  
    f2_zong(t + 1) = f2Fit(bestIndex);  
    f3_zong(t + 1) = f3Fit(bestIndex);  
    f4_zong(t + 1) = f4Fit(bestIndex);  
end

Best_Score = fMin;
BestFit = bestX';
BestIndex = bestIndex;

end
