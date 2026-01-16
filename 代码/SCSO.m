function [Best_Score,Convergence_curve,BestFit_hang,f1_zong,f2_zong,f3_zong,f4_zong]=SCSO(SearchAgents_no,dim,Max_iteration,ub,lb,pos_tag,all_param)
BestFit=zeros(2,dim/2); % 定义种群初始位置
Best_Score=Inf; % 定义最优初始值，可根据求极大值和极小值进行更换，此处为求极小值，若求极大值则换成-inf
f1_best = Inf;
f2_best = Inf;
f3_best = Inf;
f4_best = Inf;

Positions = initialization(SearchAgents_no,dim); % 生成初始种群
Positions = Positions';
Convergence_curve=zeros(1,Max_iteration);
t=0;
p=[1:360];

while t<Max_iteration
    for i=1:size(Positions,1)
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        XX=reshape(Positions(i,:),2,dim/2)';
        [Fitness,f1,f2,f3,f4]=Link_Table(XX,pos_tag,all_param);%计算每个种群个体的适应度值
        fitness_zong(i)= Fitness;
        if Fitness<Best_Score
            Best_Score=Fitness;
            f1_best=f1;
            f2_best=f2;
            f3_best=f3;
            f4_best=f4;
            BestFit=reshape(Positions(i,:),2,dim/2);
        end
    end
    S=2;                                    %S鏄渶澶ф晱鎰熻寖鍥?
    rg=S-((S)*t/(Max_iteration));                %涓?鑸伒鏁忓害鑼冨洿
   for i=1:size(Positions,1)
        r=rand*rg;                          %姣忓彧鐚笉鍚岀殑鐏垫晱搴﹁寖鍥?
        R=((2*rg)*rand)-rg;                 %鎺у埗杞崲闃舵
        Current_individual=reshape(Positions(i,:),2,dim/2);
        for j=1:dim/2
        teta=RouletteWheelSelection(p);
           if((-1<=R)&&(R<=1))              %R鍊煎湪-1鍒?1涔嬮棿
                Rand_position=abs(rand*BestFit(:,j)-Current_individual(:,j));
                Current_individual(:,j)=BestFit(:,j)-r*Rand_position*cos(teta);
           else                 
                cp=floor(SearchAgents_no*rand()+1);
                CandidatePosition =Positions(cp,:);
                Current_individual(:,j)=r*(CandidatePosition(:,j)-rand*Current_individual(:,j));
            end
        end
        Positions(i,:)=reshape(Current_individual,1,dim);
    end
    t=t+1;
    Convergence_curve(t)=Best_Score;
    f1_zong(t)=f1_best;
    f2_zong(t)=f2_best;
    f3_zong(t)=f3_best;
    f4_zong(t)=f4_best;
    disp(['SCSO代数 = ', num2str(t)]);
end
BestFit_hang=reshape(BestFit,1,dim);
end