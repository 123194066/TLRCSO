pos_reader = reshape(ESCSO_BestFit,2,[]);


% % % % % 看链路类型在linktable里面打断点看中间变量
pos_reader = pos_reader';
% pos_reader = ESCSO_BestFit;
all_param = all_power_param;
pos_tag = pos_tag;
[F,f1,f2,f3,f4,alpha_m]=Link_Table(pos_reader,pos_tag,all_param);