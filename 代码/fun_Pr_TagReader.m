function [Pr_TagFromReader, Pr_ReaderFromTag] = fun_Pr_TagReader(pos_reader,pos_tag,all_param)
% 计算标签与阅读器之间的双向接收功率
% pos_reader,pos_tag 每一行是一个二维坐标
% P为一个矩阵，每一个元素是对应行标签到对应列的接收到的功率

% 参数
z_t =     all_param.z_t;% 标签位置1
z =       all_param.z; % 阅读器位置3.5
G_T =     all_param.G_T;
G_R =     all_param.G_R;
miu_T =   all_param.miu_T;
rho_1 =   all_param.rho_1;
lamda_U = all_param.lamda_U;
gamma_1 = all_param.gamma_1;
tao =     all_param.tao;
Pt_Pow_Reader2Tag = all_param.Pt_Pow_Reader2Tag;

z_diff = z - z_t;
LOS_dis_tag_reader_z = zeros(size(pos_tag,1),size(pos_reader,1)); % 标签到阅读器直达距离（每一行是一个标签到所有阅读器的距离）
LOS_dis_tag_reader = zeros(size(pos_tag,1),size(pos_reader,1));% 标签到阅读器水平距离
for i = 1:size(pos_tag,1)
    xy_dist = vecnorm(pos_tag(i,:)-pos_reader,2,2);
    LOS_dis_tag_reader(i,:) = xy_dist;
    LOS_dis_tag_reader_z(i,:) = sqrt(xy_dist.^2 + z_diff^2)';
end

FLOOR_dis_tag_reader = zeros(size(pos_tag,1),size(pos_reader,1)); % 标签到阅读器地板反射距离，（多径）（每一行是一个标签到所有阅读器的距离）
for i = 1:size(pos_tag,1)
    theta = atan((z + z_t)./LOS_dis_tag_reader(i,:)); % 多径出发角
    FLOOR_dis_tag_reader(i,:) = LOS_dis_tag_reader(i,:)./cos(theta);
end

% 计算L
L = (lamda_U./4./pi./LOS_dis_tag_reader_z).^2 ...
     .*(abs(1+gamma_1.*LOS_dis_tag_reader_z./FLOOR_dis_tag_reader.* ...
     exp(-1j.*2.*pi.*(FLOOR_dis_tag_reader-LOS_dis_tag_reader_z)./lamda_U) )).^2;

% 计算接收功率
Pr_TagFromReader = Pt_Pow_Reader2Tag.*G_T.*G_R.*miu_T.*rho_1.*L;

Pr_ReaderFromTag = Pt_Pow_Reader2Tag.*G_T.^2.*G_R.^2.*tao.*miu_T.*rho_1.*L.^2;

end

