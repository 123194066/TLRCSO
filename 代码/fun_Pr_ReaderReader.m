function Pr_ReaderwithReader = fun_Pr_ReaderReader(pos_reader,all_param)
%FUN_PR_READERREADER 此处显示有关此函数的摘要
%   此处显示详细说明
% 参数(和论文一致)
% z_t =     all_param.z_t;% 标签位置1
% z_a =     all_param.z_a;% 终端位置1
z =       all_param.z; % 阅读器位置3.5
G_A       = all_param.G_A;    
G_W_R     = all_param.G_W_R;  
miu_A     = all_param.miu_A;  
rho_2     = all_param.rho_2;  
lamda_U   = all_param.lamda_U;
gamma_1   = all_param.gamma_1;
Pt_Reader2Reader = all_param.Pt_Reader2Reader;

pos_reader_2 = pos_reader;

LOS_dis_reader_reader = zeros(size(pos_reader_2,1),size(pos_reader,1)); % 阅读器到阅读器直达距离（每一行是一个阅读器到所有阅读器的距离）
for i = 1:size(pos_reader_2,1)
    
    LOS_dis_reader_reader(i,:) = vecnorm(pos_reader_2(i,:)-pos_reader,2,2)';
end

FLOOR_dis_reader_reader = zeros(size(pos_reader_2,1),size(pos_reader,1)); % 阅读器到阅读器地板反射距离，（多径）（每一行是一个阅读器到所有阅读器的距离）
for i = 1:size(pos_reader_2,1)
    theta = atan(2*z./LOS_dis_reader_reader(i,:)); % 多径出发角
    FLOOR_dis_reader_reader(i,:) = LOS_dis_reader_reader(i,:)./cos(theta);
end

% 计算L
L = (lamda_U./4./pi./LOS_dis_reader_reader).^2 ...
     .*(abs(1+gamma_1.*LOS_dis_reader_reader./FLOOR_dis_reader_reader.* ...
     exp(-1j.*2.*pi.*(FLOOR_dis_reader_reader-LOS_dis_reader_reader)./lamda_U) )).^2;

% 计算接收功率
Pr_ReaderwithReader = Pt_Reader2Reader.*G_A.*G_W_R.*miu_A.*rho_2.*L;


end


