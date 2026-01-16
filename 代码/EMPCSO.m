function [Best_Score, Convergence_curve, BestFit, f1_zong, f2_zong, f3_zong, f4_zong, BestIndex] = EMPCSO(SearchAgents_no, dim, Max_iteration, ub, lb, pos_tag, all_param)
%本文贡献度
r_num = dim / 2;
Tma = Max_iteration - 1;
G = 2;
hPercent = 0.8;
rPercent = 0.1;
NWmax = 2;
NWmin = 0.5;
mPercent = 0.5;

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
    [W, f1, f2, f3, f4,temp_alpha_m] = Link_Table(PP, pos_tag, all_param);
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

for t = 1:Tma
    disp(['EMPCSO 代数 = ', num2str(t)]);
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
        mother = motherLib(randi(mNum, cNum, 1));
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
            [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);

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
    aa = 2 * exp((1 - t) / 60);
    for i = (rNum + 1):(rNum + hNum)
        BB = mate(i - rNum, 1:3);
        [CC, ~] = sort(BB);

        for j = 1:dim
            r11 = rand(); r22 = rand();
            A1 = 2 * aa * r11 - aa; C1 = 2 * r22;
            D1 = abs(C1 * pX(sortIndex(CC(1)), j) - x(sortIndex(i), j));
            X1 = pX(sortIndex(CC(1)), j) - A1 * D1;

            r11 = rand(); r22 = rand();
            A2 = 2 * aa * r11 - aa; C2 = 2 * r22;
            D2 = abs(C2 * pX(sortIndex(CC(2)), j) - x(sortIndex(i), j));
            X2 = pX(sortIndex(CC(2)), j) - A2 * D2;

            r11 = rand(); r22 = rand();
            A3 = 2 * aa * r11 - aa; C3 = 2 * r22;
            D3 = abs(C3 * pX(sortIndex(CC(3)), j) - x(sortIndex(i), j));
            X3 = pX(sortIndex(CC(3)), j) - A3 * D3;

            x(sortIndex(i), j) = (X1 + X2 + X3) / 3;
        end

        x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), xmi, xma);
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
