function [YX,YY,Pg,f1_zong,f2_zong,f3_zong,f4_zong]=PSO(Size,dimension,Tmax,xmax,xmin,pos_tag,all_param)
r_num=dimension/2;
c1=2.0;c2=2.0; %学习因子
w = 0.6;
%%初始化位置和速度
Vmax=[4 4];  Vmin=[-4 -4];    %sphere函数的位置和速度的限制
Vmax=repmat(Vmax,[1,dimension/2])';
Vmin=repmat(Vmin,[1,dimension/2])';
%Xmax=600; -600;
%f_asd=zeros(Tmax,5);

TOP=[];
TOP=initialization(Size,dimension);
X=TOP;
for i=1:Size
    V(:,i) = Vmin+(Vmax-Vmin).*rand(dimension,1);
end   %初始化速度
%计算适应值，找出个体最优Pi和全局最优Pg
Pi=X; %将当前值设为个体历史最优
% Pg=zeros(dimension,1);%找出全局最优，先设为一列
for j=1:Size
    XX=reshape(X(:,j),2,dimension/2)';
    [W,f1,f2,f3,f4]=Link_Table(XX,pos_tag,all_param);
    f_x(j)=W;%计算每个粒子的适应值
    f1_x(j)=f1;
    f2_x(j)=f2;
    f3_x(j)=f3;
    f4_x(i)=f4;
    f_xbest(j)=f_x(j);              %个体最优适应值
end
[f_g,I]=min(f_xbest);     %f_g为全局最优的适应值，I为最小粒子的具体位置
Pg=X(:,I); %最优粒子,为第I行
f1_zong(1)=f1_x(I);
f2_zong(1)=f2_x(I);
f3_zong(1)=f3_x(I);
f4_zong(1)=f4_x(I);
f_gg(1) =f_g;
XX1=reshape(X(:,I),2,dimension/2)';
% f_asd(1,:)=FOUN(AA,XX1);



for t=2:Tmax
    disp(['PSO代数','=',num2str(t)]);
    time(t)=t;
    r1=rand(1);r2=rand(1);%0-1之间的随机数
    for j=1:Size
        V(:,j) =w*V(:,j)+c1*r1*(Pi(:,j)-X(:,j))+c2*r2*(Pg-X(:,j));  %速度更新
    end
    %速度限制

    for j=1:Size
        for i=1:dimension
            if V(i,j)>Vmax(i)
                V(i,j)=Vmax(i);
            elseif V(i,j)<Vmin(i)
                V(i,j)=Vmin(i);
            else
            end
        end
    end
    X=X+V;  %位置更新
    %位置限制
    for j=1:Size
        for i=1:dimension
            if X(i,j)>xmax(i)
                X(i,j)=xmax(i);
            elseif X(i,j)<xmin(i)
                X(i,j)=xmin(i);
            else
            end
        end
    end

    %重新计算适应值
    kkk=0;
    for j=1:Size
        XX=reshape(X(:,j),2,dimension/2)';
        [W,f1,f2,f3,f4]=Link_Table(XX,pos_tag,all_param);
        f_x(j)=W;%计算每个粒子的适应值
        if f_x(j)<f_xbest(j)
            f_xbest(j)=f_x(j);
            Pi(:,j)=X(:,j);
            f1_x(j)=f1;
            f2_x(j)=f2;
            f3_x(j)=f3;
            f4_x(j)=f4;

        end

    end
    [f_g,I]=min(f_xbest);     %f_g为全局最优的适应值，I为最小粒子的具体位置
    Pg=Pi(:,I); %最优粒子,为第I行
    f_gg(t)=f_g;%f_gg 每一代全局最优的适应值矩阵
    f1_zong(t)=f1_x(I);
    f2_zong(t)=f2_x(I);
    f3_zong(t)=f3_x(I);
    f4_zong(t)=f4_x(I);
end
YX=f_g;
YY=f_gg;
end




