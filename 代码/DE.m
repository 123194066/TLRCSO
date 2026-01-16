%以函数f(x,y)=-20e^((-0.2（x^2+y^2）/2))^(1/2))-e^((cos2πx+cos2πy)/2)+e为
function [Best_Score, convergence_curve, Best_pos, f1_hist, f2_hist, f3_hist, f4_hist] = DE(SearchAgents_no,dimension,Max_iteration,MaxX,MinX,pos_tag,all_param)
% SearchAgents_no=50;%群体个数
% Codel=2;%所求的变量个数
%  MinX(1)=-5;%未知量范围
%  MinX(2)=-5;
%  MaxX(1)=5;
%  MaxX(2)=5;
%  G=1;%迭代次数

F=1.2;%变异因子[0 2]
cr=0.8;%交叉因子[0.6 0.9]
r_num = dimension/2;

%初始化种群
P = initialization(SearchAgents_no,dimension);
P = P';


for i=1:SearchAgents_no
    PP = reshape(P(i,:), 2, r_num)';
    [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
    fit(i) = W;
    f1_x(i) = f1;
    f2_x(i) = f2;
    f3_x(i) = f3;
    f4_x(i) = f4;
end
[Sortfit, Index] = sort(fit);
SortP = P(Index,:);
Best_pos = SortP(1,:);%全局最优个体 之后不断更新 


Best_Score = Sortfit(1); % 不是C语言 一定要记得给初始变量否则程序跑飞

%%进入循环直到满足精度要求或者迭代次数达到
for Kg=1:1:Max_iteration
    time(Kg)=Kg;
    %第二步 变异
    for i=1:SearchAgents_no
        r1=1;r2=1;r3=1;r4=1;%使得个体满足变异条件
        while(r1==r2||r1==r3||r1==r4||r2==r3||r2==r4||r3==r4||r1==i||r2==i||r3==i||r4==i)
            r1=ceil(SearchAgents_no*rand(1));%大小匹配 
            r2=ceil(SearchAgents_no*rand(1));
            r3=ceil(SearchAgents_no*rand(1));
            r4=ceil(SearchAgents_no*rand(1));
        end

        h(i,:)=P(r1,:)+F*(P(r2,:)-P(r3,:));
        %h(i,:)=Best+F*(P(r2,:)-P(r3,:));
        for j=1:dimension %检查是否越界
            if(h(i,j)<MinX)
                h(i,j)=MinX;
            elseif(h(i,j)>MaxX) 
                h(i,j)=MaxX;
            end
        end
        %交叉
        for j=1:dimension
            temper=rand(1);
            if(temper<cr)
                v(i,j)=h(i,j);
            else
                v(i,j)=P(i,j);
            end
        end
        
        %选择
        PP = reshape(v(i,:), 2, r_num)';
        [W, f1, f2, f3, f4, ~]=Link_Table(PP,pos_tag,all_param);
        if W<fit(i)
            P(i,:)=v(i,:);
            fit(i) = W;
            f1_x(i) = f1;
            f2_x(i) = f2;
            f3_x(i) = f3;
            f4_x(i) = f4;
        end
          
    end

    [Sortfit, Index] = sort(fit);
    SortP = P(Index,:);
    Sortf1 = f1_x(Index);
    Sortf2 = f2_x(Index);
    Sortf3 = f3_x(Index);
    Sortf4 = f4_x(Index);

    Best_pos = SortP(1,:);%全局最优个体 之后不断更新 
    Best_Score = Sortfit(1);
    
    convergence_curve(Kg)=Best_Score;    
    f1_hist(Kg) = Sortf1(1);
    f2_hist(Kg) = Sortf2(1);
    f3_hist(Kg) = Sortf3(1);
    f4_hist(Kg) = Sortf4(1);

end

%   fprintf('最优解结果为%f,%f',Best(1),Best(2));
%    fprintf('最大函数值为%f',Best_f(Kg));
%    plot(time,Best_f(time));
