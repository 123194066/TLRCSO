function [best_score,convergence_curve,best_pos, f1_hist, f2_hist, f3_hist, f4_hist]=GA(SearchAgents_no,dimension,Max_iteration,ub,lb,pos_tag,all_param)
% dim=10;
% SearchAgents_no = 100;                 %染色体数目
% Xs=20;
% Xx=-20;
% Max_iteration = 100;                  %最大遗传代数

r_num = dimension/2;
nf = zeros(dimension,SearchAgents_no);         %子种群赋空间
Pc = 0.8;                 %交叉概率
Pm = 0.1;                 %变异概率
x = initialization(SearchAgents_no,dimension);%随机获得初始种群
x = x';
% 计算初始适应度
fit = zeros(1, SearchAgents_no);
f1_x = zeros(1, SearchAgents_no);
f2_x = zeros(1, SearchAgents_no);
f3_x = zeros(1, SearchAgents_no);
f4_x = zeros(1, SearchAgents_no);
%%%%%%%%%%%按适应度升序排列%%%%%%%%%%%%%%%%%
for np=1:SearchAgents_no
    PP = reshape(x(np,:), 2, r_num)';
    [W, f1, f2, f3, f4, ~]=Link_Table(PP,pos_tag,all_param);
    fit(np) = W;
    f1_x(np) = f1;
    f2_x(np) = f2;
    f3_x(np) = f3;
    f4_x(np) = f4;
end
[SortFIT,Index]=sort(fit);
Sortf1 = f1_x(Index);
Sortf2 = f2_x(Index);
Sortf3 = f3_x(Index);
Sortf4 = f4_x(Index);

Sortf=x(Index,:);
Sortf = Sortf';




%%%%%%%%%%遗传算法循环%%%%%%%%%%%%%%%%%%%%%
for gen=1:Max_iteration
    %%%%%%%%%%%%%采用君主方案进行交叉操作%%%%%%%%%%%%%%
    Emper=Sortf(:,1);
    NoPoint=round(dimension*Pc);
    PoPoint=randi([1 dimension],NoPoint,SearchAgents_no/2);
    nf=Sortf;
    for i=1:SearchAgents_no/2
        nf(:,2*i-1)=Emper;
        nf(:,2*i)=Sortf(:,2*i);
        for k=1:NoPoint
            nf(PoPoint(k,i),2*i-1)=nf(PoPoint(k,i),2*i);
            nf(PoPoint(k,i),2*i)=Emper(PoPoint(k,i));
        end
    end
    %%%%%%%%%%%%变异操作%%%%%%%%%%%%%%%%%%%%%%
    for m=1:SearchAgents_no
        for n=1:dimension
            r=rand(1,1);
            if r<Pm
                nf(n,m)=rand(1,1)*(ub-lb)+lb;
            end
        end
    end
    %%%%%%%%%%子种群按适应度升序排列%%%%%%%%%%%%%%%%%
    for np=1:SearchAgents_no
        PP = reshape(nf(:,np), 2, r_num)';
        [W, f1, f2, f3, f4, ~]=Link_Table(PP,pos_tag,all_param);
        NFIT(np) = W;
        Nf1_x(np) = f1;
        Nf2_x(np) = f2;
        Nf3_x(np) = f3;
        Nf4_x(np) = f4;
    end

    [NSortFIT,Index]=sort(NFIT);
    NSortf=nf(:,Index);
    NSortf1 = Nf1_x(Index);
    NSortf2 = Nf2_x(Index);
    NSortf3 = Nf3_x(Index);
    NSortf4 = Nf4_x(Index);

    %%%%%%%%%%%产生新种群%%%%%%%%%%%%
    f_all=[Sortf,NSortf]; % 合并两次迭代
    f1_all = [Sortf1,NSortf1];
    f2_all = [Sortf2,NSortf2];
    f3_all = [Sortf3,NSortf3];
    f4_all = [Sortf4,NSortf4];

    FIT1=[SortFIT,NSortFIT];
    [SortFIT1,Index]=sort(FIT1);
    Sortf_all=f_all(:,Index);
    Sortf1_all = f1_all(Index);
    Sortf2_all = f2_all(Index);
    Sortf3_all = f3_all(Index);
    Sortf4_all = f4_all(Index);


    SortFIT=SortFIT1(1:SearchAgents_no); % 取前"子代数目"个子代
    Sortf=Sortf_all(:,1:SearchAgents_no);
    Sortf1 = Sortf1(1:SearchAgents_no);
    Sortf2 = Sortf2(1:SearchAgents_no);
    Sortf3 = Sortf3(1:SearchAgents_no);
    Sortf4 = Sortf4(1:SearchAgents_no);


    convergence_curve(gen)=SortFIT(1); % 记录曲线
    f1_hist(gen) = Sortf1(1);
    f2_hist(gen) = Sortf2(1);
    f3_hist(gen) = Sortf3(1);
    f4_hist(gen) = Sortf4(1);
    
end
best_pos=Sortf(:,1)';
best_score = SortFIT(1);
% trace(end);
% figure
% plot(trace)
% xlabel('迭代次数')
% ylabel('目标适应值')
% title('适应度进化曲线')
end
%%%%%%%%%%适应度函数%%%%%%%%%%%%%%%

% function result=fobj(x)
% summ=sum(x.^2);
% result=summ;
% end



            
            
            
            
            
            
            
            
            
            