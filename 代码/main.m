clear;clc;
close all;
tic
 warning off
 % parpool('local',4)
N = 16; % 16个阅读器天线
dim = N*2;
SearchAgents_no = 100; % 子代个数
Max_iteration = 500;
ub1 = [30 30];
lb1 = [0 0];
ub1 = repmat(ub1,[1,N])';
lb1 = repmat(lb1,[1,N])';
ub2 = 30;
lb2 = 0;

% 参数结构体
all_power_param.z_t = 1;% 标签位置
all_power_param.z = 3.5;% 阅读器位置
all_power_param.z_a = 1; % 终端位置
all_power_param.G_T = 1.58;
all_power_param.G_R = 4;
all_power_param.G_A = 4;%不确定
all_power_param.G_W_R = 4;
all_power_param.miu_T = 1;
all_power_param.miu_A = 0.3;
all_power_param.rho_1 = 0.5;
all_power_param.rho_2 = 0.5;
all_power_param.lamda_U = 1/3;
all_power_param.gamma_1 = -real(sqrt(0.1));
all_power_param.tao = 0.5;
all_power_param.Pt_Pow_Reader2Tag = 0.5;
all_power_param.Pt_Reader2Reader = 0.5;
all_power_param.Pt_ReaderwithAppterminal = 0.5;
all_power_param.S1 = 10^-5;
all_power_param.S2 = 10^-14;
all_power_param.S3 = 10^-5.5;
all_power_param.S4 = 10^-5.5;


% 场景1：EMPCSOI45（无序，45个标签）
A = [3*ones(1,9);7:2:23];
B = [5:2:27 ;23*ones(1,12)];
D = [19.5:2:26;11*ones(1,4)];
H = [19.5:2:26;17*ones(1,4)];
C = [18*ones(1,8);7:2:22];
F = [27*ones(1,8);7:2:22];
pos_tag = [A B C D F H ]'; % 45×2矩阵

% 场景2：EMPCSOR45（有序，45个标签）
%A = [3*ones(1,9);7:2:23];
%D = [11*ones(1,9);7:2:23];
%C = [19*ones(1,9);7:2:23];
%F = [27*ones(1,9);7:2:23];
%B = [5:2:10 ;23*ones(1,3)];
%E = [13:2:18 ;7*ones(1,3)];
%G = [21:2:26 ;23*ones(1,3)];
%pos_tag = [A B C D E F G]'; % 45×2矩阵

% 场景3：EMPCSOI59（无序，59个标签）
%A = [3*ones(1,11);7:1.5:23];
%B = [3.5:1.5:27 ;23*ones(1,16)];
%D = [19.5:1.5:26;11*ones(1,5)];
%H = [19.5:1.5:26;17*ones(1,5)];
%C = [18*ones(1,11);22:-1.5:7];
%F = [27*ones(1,11);22:-1.5:7];
%pos_tag = [A B C D H F ]'; % 59×2矩阵

% 场景4：EMPCSOR59（有序，59个标签）
%A = [3*ones(1,11);7:1.5:23];
%D = [11*ones(1,11);23:-1.5:7];
%C = [19*ones(1,11);7:1.5:23];
%F = [27*ones(1,11);22.5:-1.5:7];
%B = [3.5:1.5:10 ;23*ones(1,5)];
%E = [11.5:1.5:18 ;7*ones(1,5)];
%G = [19.5:1.5:26 ;23*ones(1,5)];
%pos_tag = [A B C D E F G]'; % 59×2矩阵

% 场景5：EMPCSOI89（无序，89个标签）
%A = [3*ones(1,17);7:1:23];
%D = [19:1:26;11*ones(1,8)];
%H = [19:1:26;17*ones(1,8)];
%C = [18*ones(1,17);7:1:23];
%F = [27*ones(1,17);7:1:23];
%B = [4:1:11 ;23*ones(1,8)];
%E = [12:1:17 ;23*ones(1,6)];
%G = [19:1:26 ;23*ones(1,8)];
%pos_tag = [A B C D E F H G]'; % 89×2矩阵

% 场景6：EMPCSOR89（有序，89个标签）
%A = [3*ones(1,17);7:1:23];
%D = [11*ones(1,17);7:1:23];
%C = [19*ones(1,17);7:1:23];
%F = [27*ones(1,17);7:1:23];
%B = [4:1:10 ;23*ones(1,7)];
%E = [12:1:18 ;7*ones(1,7)];
%G = [20:1:26 ;23*ones(1,7)];
%pos_tag = [A B C D E F G]'; % 89×2矩阵

% % 初始化时间记录数组

% 
% time_PSO = zeros(1,50);
% time_CSO = zeros(1,50);
% time_GWO = zeros(1,50);
% time_EMPCSO = zeros(1,50);
% time_SCSO = zeros(1,50);
% time_MPCSO = zeros(1,50);
% time_BSA = zeros(1,50);
time_TLRCSO = zeros(1,50);
% time_EMPCSO1 = zeros(1,50);
% time_EMPCSO2 = zeros(1,50);
% time_EMPCSO3 = zeros(1,50);
% time_EMPCSO12 = zeros(1,50);
% time_EMPCSO13 = zeros(1,50);
% time_EMPCSO23 = zeros(1,50);
% 

% 开启并行池（关键代码）
if isempty(gcp('nocreate'))  % 检查是否已存在并行池，不存在则创建
    parpool('local');  % 'local'表示使用本地计算机的核心
end

parfor i = 1:50
%     t_start = tic;
%     [PSO_score(i),PSO_cg_curve(i,:),PSO_BestFit(i,:),PSO_f1(i,:),PSO_f2(i,:),PSO_f3(i,:),PSO_f4(i,:)] = PSO(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param); % run PSO to compare to results
%     time_PSO(i) = toc(t_start);
%     
%     t_start = tic;
%     [CSO_score(i),CSO_cg_curve(i,:),CSO_BestFit(i,:),CSO_f1(i,:),CSO_f2(i,:),CSO_f3(i,:),CSO_f4(i,:)] = CSO(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_CSO(i) = toc(t_start);
%     
%      t_start = tic;
%     [GWO_score(i),GWO_cg_curve(i,:),GWO_BestFit(i,:),GWO_f1(i,:),GWO_f2(i,:),GWO_f3(i,:),GWO_f4(i,:)] = GWO(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_GWO(i) = toc(t_start);
%     
%      t_start = tic;
%     [EMPCSO_score(i), EMPCSO_cg_curve(i,:), EMPCSO_BestFit(i,:), EMPCSO_f1(i,:), EMPCSO_f2(i,:), EMPCSO_f3(i,:), EMPCSO_f4(i,:)] = EMPCSO(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_EMPCSO(i) = toc(t_start);
%     
%      t_start = tic;
%     [SCSO_score(i),SCSO_cg_curve(i,:),SCSO_BestFit(i,:),SCSO_f1(i,:),SCSO_f2(i,:),SCSO_f3(i,:),SCSO_f4(i,:)] = SCSO(SearchAgents_no,dim,Max_iteration,ub2,lb2,pos_tag,all_power_param);
%     time_SCSO(i) = toc(t_start);
%     
%     t_start = tic;
%     [MPCSO_score(i),MPCSO_cg_curve(i,:),MPCSO_BestFit(i,:),MPCSO_f1(i,:),MPCSO_f2(i,:),MPCSO_f3(i,:),MPCSO_f4(i,:)] = MPCSO(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_MPCSO(i) = toc(t_start);
%     
%     t_start = tic;
%     [BSA_score(i),BSA_cg_curve(i,:),BSA_BestFit(i,:),BSA_f1(i,:),BSA_f2(i,:),BSA_f3(i,:),BSA_f4(i,:)] = BSA(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_BSA(i) = toc(t_start);

   t_start = tic;
    [TLRCSO_score(i), TLRCSO_cg_curve(i,:), TLRCSO_BestFit(i,:), TLRCSO_f1(i,:), TLRCSO_f2(i,:), TLRCSO_f3(i,:), TLRCSO_f4(i,:)] = TLRCSO(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
    time_TLRCSO(i) = toc(t_start);
    
%      t_start = tic;
%      [EMPCSO1_score(i), EMPCSO1_cg_curve(i,:), EMPCSO1_BestFit(i,:), EMPCSO1_f1(i,:), EMPCSO1_f2(i,:), EMPCSO1_f3(i,:), EMPCSO1_f4(i,:)] = EMPCSO1(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_EMPCSO1(i) = toc(t_start);
%     
%       t_start = tic;
%    [EMPCSO2_score(i), EMPCSO2_cg_curve(i,:), EMPCSO2_BestFit(i,:), EMPCSO2_f1(i,:), EMPCSO2_f2(i,:), EMPCSO2_f3(i,:), EMPCSO2_f4(i,:)] = EMPCSO2(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_EMPCSO2(i) = toc(t_start);
%     
%       t_start = tic;
%       [EMPCSO3_score(i), EMPCSO3_cg_curve(i,:), EMPCSO3_BestFit(i,:), EMPCSO3_f1(i,:), EMPCSO3_f2(i,:), EMPCSO3_f3(i,:), EMPCSO3_f4(i,:)] = EMPCSO3(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_EMPCSO3(i) = toc(t_start);
%     
%       t_start = tic;
%     [EMPCSO12_score(i), EMPCSO12_cg_curve(i,:), EMPCSO12_BestFit(i,:), EMPCSO12_f1(i,:), EMPCSO12_f2(i,:), EMPCSO12_f3(i,:), EMPCSO12_f4(i,:)] = EMPCSO12(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_EMPCSO12(i) = toc(t_start);
%     
%      t_start = tic;
%      [EMPCSO13_score(i), EMPCSO13_cg_curve(i,:), EMPCSO13_BestFit(i,:), EMPCSO13_f1(i,:), EMPCSO13_f2(i,:), EMPCSO13_f3(i,:), EMPCSO13_f4(i,:)] = EMPCSO13(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_EMPCSO13(i) = toc(t_start);
%     
%       t_start = tic;
%    [EMPCSO23_score(i), EMPCSO23_cg_curve(i,:), EMPCSO23_BestFit(i,:), EMPCSO23_f1(i,:), EMPCSO23_f2(i,:), EMPCSO23_f3(i,:), EMPCSO23_f4(i,:)] = EMPCSO23(SearchAgents_no,dim,Max_iteration,ub1,lb1,pos_tag,all_power_param);
%     time_EMPCSO23(i) = toc(t_start);
%     
    end
%    PSO_cg_curve_mean=real(mean(PSO_cg_curve,1));
%    CSO_cg_curve_mean=real(mean(CSO_cg_curve,1));
%    GWO_cg_curve_mean=real(mean(GWO_cg_curve,1));
%    EMPCSO_cg_curve_mean=real(mean(EMPCSO_cg_curve,1));
%    SCSO_cg_curve_mean=real(mean(SCSO_cg_curve,1));
%    MPCSO_cg_curve_mean=real(mean(MPCSO_cg_curve,1));
%    BSA_cg_curve_mean=real(mean(BSA_cg_curve,1));

                TLRCSO_cg_curve_mean=real(mean(TLRCSO_cg_curve,1));
%    EMPCSO1_cg_curve_mean=real(mean(EMPCSO1_cg_curve,1));
%        EMPCSO2_cg_curve_mean=real(mean(EMPCSO2_cg_curve,1));
%       EMPCSO3_cg_curve_mean=real(mean(EMPCSO3_cg_curve,1));
%          EMPCSO12_cg_curve_mean=real(mean(EMPCSO12_cg_curve,1));
%    EMPCSO13_cg_curve_mean=real(mean(EMPCSO13_cg_curve,1));
%        EMPCSO23_cg_curve_mean=real(mean(EMPCSO23_cg_curve,1));




figure(1)
   hold on
%   semilogy(PSO_cg_curve_mean, "Color", "#1E90FF", "LineWidth", 1.5);   % 道奇蓝 - PSO
% semilogy(CSO_cg_curve_mean, "Color", "#FF0000", "LineWidth", 1.5);          % 纯红 - CSO
% semilogy(GWO_cg_curve_mean, "Color", "#228B22", "LineWidth", 1.5);          % 森林绿 - GWO
% semilogy(EMPCSO_cg_curve_mean, "Color", "#32CD32", "LineWidth", 1.5);       %  lime绿 - EMPCSO
% semilogy(SCSO_cg_curve_mean, "Color", "#9370DB", "LineWidth", 1.5);         % 中紫色 - SCSO
% semilogy(MPCSO_cg_curve_mean, "Color", "#FF6347", "LineWidth", 1.5);        % 番茄红 - MPCSO
% semilogy(BSA_cg_curve_mean, "Color", "#FFD700", "LineWidth", 1.5);          % 金色 - BSA
semilogy(TLRCSO_cg_curve_mean, "Color", "#000000", "LineWidth", 1.5);        % 黑色 - ESCSO
% semilogy(EMPCSO1_cg_curve_mean, "Color", "#4682B4", "LineWidth", 1.5);   % 钢青蓝 - EMPCSO_new
% semilogy(EMPCSO2_cg_curve_mean, "Color", "#006400", "LineWidth", 1.5);    % 墨绿 - MPCSO_new
%  semilogy(EMPCSO3_cg_curve_mean, "Color", "#DC143C", "LineWidth", 1.5);    % 猩红 - ESCSO_new
%   semilogy(EMPCSO12_cg_curve_mean, "Color", "#1E90FF", "LineWidth", 1.5);   % 道奇蓝 - PSO
% semilogy(EMPCSO13_cg_curve_mean, "Color", "#FF0000", "LineWidth", 1.5);          % 纯红 - CSO
% semilogy(EMPCSO23_cg_curve_mean, "Color", "#228B22", "LineWidth", 1.5);          % 森林绿 - GWO  semilogy(PSO_cg_curve_mean, "Color", "#1E90FF", "LineWidth", 1.5);   % 道奇蓝 - PSO
 title('收敛曲线')
xlabel('迭代次数');
ylabel('适应度值');
 h=legend('PSO','CSO','GWO','EMPCSO','SCSO','MPCSO','BSA','TLRCSO','EMPCSO1','EMPCSO2','EMPCSO3','EMPCSO12','EMPCSO13','EMPCSO23');

%  h=legend('PSO','CSO','GWO','EMPCSO','SCSO','MPCSO','BSA','TLRCSO','EMPCSO1','EMPCSO2','EMPCSO3');
set(h,'FontName','宋体','FontSize',13);
hold off

% 显示各算法的平均运行时间
% fprintf('各算法平均运行时间（秒）：\n');
% fprintf('PSO: %.4f\n', mean(time_PSO));
% fprintf('CSO: %.4f\n', mean(time_CSO));
% fprintf('GWO: %.4f\n', mean(time_GWO));
% fprintf('EMPCSO: %.4f\n', mean(time_EMPCSO));
% fprintf('SCSO: %.4f\n', mean(time_SCSO));
% fprintf('MPCSO: %.4f\n', mean(time_MPCSO));
% fprintf('BSA: %.4f\n', mean(time_BSA));
% fprintf('TLRCSO: %.4f\n', mean(time_TLRCSO));
% fprintf('EMPCSO1: %.4f\n', mean(time_EMPCSO1));
% fprintf('EMPCSO2: %.4f\n', mean(time_EMPCSO2));
% fprintf('EMPCSO3: %.4f\n', mean(time_EMPCSO3));
% fprintf('EMPCSO12: %.4f\n', mean(time_EMPCSO12));
% fprintf('EMPCSO13: %.4f\n', mean(time_EMPCSO13));
% fprintf('EMPCSO23: %.4f\n', mean(time_EMPCSO23));
% 
% 

% 可选：将时间结果保存到表格
% time_results = table( ...
%     mean(time_PSO), mean(time_CSO), mean(time_GWO), mean(time_EMPCSO), ...
%     mean(time_SCSO), mean(time_MPCSO), mean(time_BSA), ...
%     mean(time_TLRCSO), mean(time_EMPCSO1), mean(time_EMPCSO2), mean(time_EMPCSO3), ...
%     'RowNames', {'平均时间(秒)'}, ...
%     'VariableNames', {'PSO','CSO','GWO','EMPCSO','SCSO','MPCSO','BSA','TLRCSO','EMPCSO1','EMPCSO2','EMPCSO3'});
% disp(time_results);
% time_results = table( ...
%      mean(time_PSO), mean(time_CSO), mean(time_GWO), mean(time_EMPCSO), ...
%     mean(time_SCSO), mean(time_MPCSO), mean(time_BSA), ...
%     mean(time_TLRCSO), mean(time_EMPCSO1), mean(time_EMPCSO2), ...
%     mean(time_EMPCSO3), mean(time_EMPCSO12), mean(time_EMPCSO13), mean(time_EMPCSO23), ...
%     'RowNames', {'平均时间(秒)'}, ...
%     'VariableNames', {'PSO','CSO','GWO','EMPCSO','SCSO','MPCSO','BSA','TLRCSO','EMPCSO1','EMPCSO2','EMPCSO3','EMPCSO12','EMPCSO13','EMPCSO23'});
% disp(time_results);
% total_time = toc;
% fprintf('总运行时间：%.4f秒\n', total_time);
% 
