function Pr_ReaderwithAppterminal = fun_Pr_ReaderAppterminal(pos_reader,pos_appterminal,all_param)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
% 参数
z_a =     all_param.z_a;% 终端位置1
z =       all_param.z; % 阅读器位置3.5
G_A     = all_param.G_A;   
G_W_R   = all_param.G_W_R;  
miu_A   = all_param.miu_A;  
rho_2   = all_param.rho_2;  
lamda_U = all_param.lamda_U;
gamma_1 = all_param.gamma_1;
Pt_ReaderwithAppterminal = all_param.Pt_ReaderwithAppterminal;

z_diff = z - z_a;
LOS_dis_appteiminal_reader_z = zeros(size(pos_appterminal,1),size(pos_reader,1)); % 终端到阅读器直达距离（每一行是一个终端到所有阅读器的距离）
LOS_dis_appteiminal_reader = zeros(size(pos_appterminal,1),size(pos_reader,1)); % 终端到阅读器水平距离（每一行是一个终端到所有阅读器的距离）
for i = 1:size(pos_appterminal,1)
    xy_dist =vecnorm(pos_appterminal(i,:)-pos_reader,2,2)';
    LOS_dis_appteiminal_reader(i,:) = xy_dist;
    LOS_dis_appteiminal_reader_z(i,:) = sqrt(xy_dist.^2 + z_diff^2)';
end

FLOOR_dis_appterminal_reader = zeros(size(pos_appterminal,1),size(pos_reader,1)); % 终端到阅读器地板反射距离，（多径）（每一行是一个终端到所有阅读器的距离）
for i = 1:size(pos_appterminal,1)
    theta = atan((z + z_a)./LOS_dis_appteiminal_reader(i,:)); % 多径出发角
    FLOOR_dis_appterminal_reader(i,:) = LOS_dis_appteiminal_reader(i,:)./cos(theta);
end

% 计算L
L = (lamda_U./4./pi./LOS_dis_appteiminal_reader_z).^2 ...
     .*(abs(1+gamma_1.*LOS_dis_appteiminal_reader_z./FLOOR_dis_appterminal_reader.* ...
     exp(-1j.*2.*pi.*(FLOOR_dis_appterminal_reader-LOS_dis_appteiminal_reader_z)./lamda_U) )).^2;

% 计算接收功率
Pr_ReaderwithAppterminal = Pt_ReaderwithAppterminal.*G_A.*G_W_R.*miu_A.*rho_2.*L;


end

