function [Best_Score, Convergence_curve, BestFit, f1_zong, f2_zong, f3_zong, f4_zong] = MPCSO(SearchAgents_no, dim, Max_iteration, ub, lb, pos_tag, all_param)
    
    r_num = dim / 2;
    Tma = Max_iteration - 1;
    G = 2;       
    hPercent = 0.8; 
    rPercent = 0.1; 
    NWmax = 2;
    NWmin = 0.5;
    mPercent = 0.5;                  
    w_max = 0.9;                       
    w_min = 0.4;

    rNum = round(SearchAgents_no * rPercent);    % 公鸡数量
    hNum = round(SearchAgents_no * hPercent);    % 母鸡数量
    cNum = SearchAgents_no - rNum - hNum;        % 小鸡数量
    mNum = round(hNum * mPercent);               % 母 hens数量
    xma = ub';
    xmi = lb';  

    % 初始化速度矩阵
    Vmax = [4 4];  Vmin = [-4 -4];  
    Vmax = repmat(Vmax, [1, r_num])';
    Vmin = repmat(Vmin, [1, r_num])';
    for i = 1:SearchAgents_no
        V(:, i) = Vmin + (Vmax - Vmin) .* rand(dim, 1);
    end  
    V = V';

    % 初始化种群
    TOP = initialization(SearchAgents_no, dim);
    x = TOP';    
    f1_x = zeros(1, SearchAgents_no);
    f2_x = zeros(1, SearchAgents_no);
    f3_x = zeros(1, SearchAgents_no);
    f4_x = zeros(1, SearchAgents_no);
    fit = zeros(1, SearchAgents_no);
    contribution = cell(1, SearchAgents_no);  % 存储每个个体的覆盖贡献度向量
    for i = 1:SearchAgents_no
        PP = reshape(x(i, :), 2, r_num)';
        [W, f1, f2, f3, f4,  temp_alpha_m] = Link_Table(PP, pos_tag, all_param);
        fit(i) = W;
        f1_x(i) = f1;
        f2_x(i) = f2;
        f3_x(i) = f3;
        f4_x(i) = f4;
        contribution{i} = temp_alpha_m;  % 保存贡献度
    end


    pFit = fit;   
    pX = x;                      
    f1Fit = f1_x;                
    f2Fit = f2_x;                
    f3Fit = f3_x;               
    f4Fit = f4_x;    
      pContribution = contribution;  % 个体历史最优对应的贡献度


    [fMin, bestIndex] = min(fit);  
    bestX = x(bestIndex, :); 

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

    % 主迭代循环
    for t = 1:Tma
        disp(['MPCSO代数=', num2str(t)]);
        FL = rand(SearchAgents_no, 1) .* 0.4 + 0.5;  % 跟随因子

        % 每G代更新种群关系
        if mod(t, G) == 1   
            [~, sortIndex] = sort(pFit);   
            x = pX;

            motherLib = randpermD(hNum, mNum) + rNum;   

            for i = 1:SearchAgents_no
                Xx(i, :) = x(sortIndex(i), :); 
            end
            oudist = pdist2(Xx, Xx);
            hdist = oudist(rNum+1:rNum+hNum+cNum, 1:rNum);
            [~, h_sortIndex] = sort(hdist');
            mate = h_sortIndex(1:2, 1:hNum)';
            mother = motherLib(randi(mNum, cNum, 1));  
        end

        % 公鸡位置更新
        for i = 1:rNum    
            NW = NWmax - (NWmax - NWmin) * t / Tma;

            % 随机选择另一只公鸡
            anotherRooster = randiTabu(1, rNum, i, 1);  
            if pFit(sortIndex(i)) <= pFit(sortIndex(anotherRooster))
                tempSigma = NW .* [0.5 0.5];
            else
                tempSigma = NW .* exp((pFit(sortIndex(anotherRooster)) - ...
                    pFit(sortIndex(i))) / (abs(pFit(sortIndex(i))) + realmin)) .* [0.5 0.5];
            end

            % 计算当前位置适应度
            PP = reshape(pX(sortIndex(i), :), 2, r_num)';
            [W, f1, f2, f3, f4,  temp_alpha_m] = Link_Table(PP, pos_tag, all_param);
            fit(sortIndex(i)) = W;
            f1_x(sortIndex(i)) = f1;
            f2_x(sortIndex(i)) = f2;
            f3_x(sortIndex(i)) = f3;
            f4_x(sortIndex(i)) = f4;
              pContribution{sortIndex(i)} = temp_alpha_m;  % 更新贡献度

            % 随机选择子维度进行扰动
            zero_contribution_idx = find(temp_alpha_m == 0);
            if ~isempty(zero_contribution_idx)
                % 存在覆盖标签数为0的单体，优先扰动
                r_index = zero_contribution_idx(randi(length(zero_contribution_idx)));
            else
                % 无覆盖为0的单体，随机选择
                r_index = randpermD(r_num, 1);
            end
            for n = 1:length(r_index)
                x(sortIndex(i), (r_index(n)-1)*2+1:r_index(n)*2) = ...
                    pX(sortIndex(i), (r_index(n)-1)*2+1:r_index(n)*2) .* ...
                    (1 + tempSigma .* randn(1, 2));
            end

            % 边界处理与适应度更新
            x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), xmi, xma);
            PP = reshape(x(sortIndex(i), :), 2, r_num)';
            [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
            fit(sortIndex(i)) = W;
            f1_x(sortIndex(i)) = f1;
            f2_x(sortIndex(i)) = f2;
            f3_x(sortIndex(i)) = f3;
            f4_x(sortIndex(i)) = f4;
        end

        % 母鸡位置更新
        for i = (rNum + 1):(rNum + hNum)  
            c1 = 2; c2 = 2;
            r1 = rand(1); r2 = rand(1);

            % 速度更新
            V(sortIndex(i), :) = (w_max - (w_max - w_min) * t / Max_iteration) * V(sortIndex(i), :) + ...
                c1 * r1 * (pX(sortIndex(mate(i - rNum, 1)), :) - x(sortIndex(i), :)) + ...
                c2 * r2 * (pX(sortIndex(i), :) - x(sortIndex(i), :));

            % 速度边界处理
            for j = 1:SearchAgents_no
                for ii = 1:dim
                    if V(j, ii) > Vmax(ii) 
                        V(j, ii) = Vmax(ii);
                    elseif V(j, ii) < Vmin(ii)
                        V(j, ii) = Vmin(ii);
                    end
                end
            end

            % 位置更新
            x(sortIndex(i), :) = x(sortIndex(i), :) + V(sortIndex(i), :);
            x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), xmi, xma);

            PP = reshape(x(sortIndex(i), :), 2, r_num)';
            [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
            fit(sortIndex(i)) = W;
            f1_x(sortIndex(i)) = f1;
            f2_x(sortIndex(i)) = f2;
            f3_x(sortIndex(i)) = f3;
            f4_x(sortIndex(i)) = f4;
        end

        % 小鸡位置更新
        for i = (rNum + hNum + 1):SearchAgents_no    
            x(sortIndex(i), :) = pX(sortIndex(i), :) + ...
                (pX(sortIndex(mother(i - rNum - hNum)), :) - pX(sortIndex(i), :)) .* FL(i);
            x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), xmi, xma);

            % 适应度计算
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
                [~, ~, ~, ~, ~, temp_alpha_m] = Link_Table(PP, pos_tag, all_param);
                pContribution{i} = temp_alpha_m;
            end

            if pFit(i) < fMin
                fMin = pFit(i);
                bestX = pX(i, :);
            end
        end

        Convergence_curve(1 + t) = fMin;
        [~, bestIndex] = min(pFit); 
        f1_zong(t + 1) = f1Fit(bestIndex);
        f2_zong(t + 1) = f2Fit(bestIndex);
        f3_zong(t + 1) = f3Fit(bestIndex);
        f4_zong(t + 1) = f4Fit(bestIndex);
    end

    Best_Score = fMin;
    Convergence_curve = Convergence_curve;
    BestFit = bestX;
end

