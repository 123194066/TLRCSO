% function F=Link_Table(pos_reader_hang)
function [F,f1,f2,f3,f4,alpha_m]=Link_Table(pos_reader,pos_tag,all_param)




%�����ն�λ����Ϣ
pos_appterminal=[35,20]; %һ������

Num_t = size(pos_tag,1); % ��ǩ����
Num_r = size(pos_reader,1); % �Ķ�������

%������������ֵ�����������һ�£�
S1=all_param.S1; % 
S2=all_param.S2;
S3=all_param.S3;
S4=all_param.S4;

% �����ǩ���Ķ���֮���˫����չ���
[Pr_TagFromReader,Pr_ReaderFromTag] = fun_Pr_TagReader(pos_reader,pos_tag,all_param);

% ����Ӧ���ն����Ķ���֮���˫����չ���
Pr_ReaderwithAppterminal = fun_Pr_ReaderAppterminal(pos_reader,pos_appterminal,all_param);

% �����Ķ������Ķ���֮���˫����չ���
Pr_ReaderwithReader = fun_Pr_ReaderReader(pos_reader,all_param);

% �����Ƿ�������ֵ�ı�־λ
flag_S1 = Pr_TagFromReader >= S1;
flag_S2 = Pr_ReaderFromTag >= S2;
flag_S3 = Pr_ReaderwithAppterminal >= S3;
flag_S3 = repmat(flag_S3,[Num_t 1]);
flag_S4 = Pr_ReaderwithReader >= S4;

% Ϊ�Ķ�������
flag_DI_Reader = flag_S1 & flag_S2 & flag_S3;
flag_TF_Reader = flag_S1 & flag_S2 & (~flag_S3);
flag_AF_Reader = ((~flag_S1) | (~flag_S2)) & flag_S3;
flag_left_Reader = ~(flag_DI_Reader | flag_TF_Reader | flag_AF_Reader); % ʣ���Ķ�����ΪѰ��RC��׼��

% Ѱ����·
Link_Type_1 = cell(Num_t,1); % ��һ��
Link_Type_2 = cell(Num_t,1); % �ڶ���
Link_Type_3 = cell(Num_t,1); % ������

% ��¼������Ӧ���������
all_single_tag_link_num = zeros(Num_t,1); % ÿ��Tag������·��Ŀ�����ڼ���f1
all_GDOP_raeder = cell(Num_t,1); % ��¼ÿ��Tag��������·�ĵ�һ��reader���
all_link1_num = 0;
all_link2_num = 0; % type2��·����
all_link3_num = 0; % type3��·����
all_uniq_reader = []; % ������·���Ķ�����ż���
all_forward_reader = []; % ����ǰ��ͨ·���Ķ�����ż���

% ��ʼ�����׶� 0�����޹���
alpha_m = zeros(1,Num_r);


for i = 1:Num_t
    DI_Reader = find(flag_DI_Reader(i,:) == 1);
    TF_Reader = find(flag_TF_Reader(i,:) == 1);
    AF_Reader = find(flag_AF_Reader(i,:) == 1);
    left_reader = find(flag_left_Reader(i,:) == 1);

    % Ѱ�ҵ�һ��
    Temp_link1 = DI_Reader';
    Link_Type_1{i} = Temp_link1;

    % Ѱ�ҵڶ���(Ϊÿһ��TF����·)
    Temp_link2 = zeros(0,1); % ÿһ��Ϊһ�����е�T2��·
    for j = 1:length(TF_Reader)
        possible_AF_reader = intersect(find(flag_S4(TF_Reader(j),:) == 1) , AF_Reader); % ��¼��ѡAF
        if isempty(possible_AF_reader);continue;end % û�ҵ�����
        % ��Q����
        Q = zeros(1,length(possible_AF_reader));
        for k = 1:length(possible_AF_reader)
            Q(k) = min([Pr_ReaderwithReader(TF_Reader(j),possible_AF_reader(k)) , Pr_ReaderwithAppterminal(possible_AF_reader(k))]);
        end
        best_AF_reader = possible_AF_reader(Q == max(Q)); % ����AF        
        Temp_link2 = [Temp_link2; [TF_Reader(j),best_AF_reader(1)]]; % ��¼       
    end
    Link_Type_2{i} = Temp_link2; % ��¼����ӦTag��

    TF_Reader_for_Type3 = setdiff(TF_Reader , Temp_link2(:,1)'); % ȥ���ѳ�ΪT2��·���Ķ�����ΪT3��ѡ����׼����

    % Ѱ�ҵ�����
    Temp_link3 = zeros(0,3);
    for j = 1:length(TF_Reader_for_Type3)
        % ��¼��ǰTF_reader����link3��·
        local_best_link3 = [];
        local_best_Q = [];
        for k = 1:length(AF_Reader)
            % ��ʼ������ÿһ�д洢����ÿ��AF������·���Ͷ�Ӧ��Qֵ
            % ��͵�ǰTF����ͬʱ�͵�ǰAF�����ı�ѡRC
            possible_RC_reader = intersect(intersect( find(flag_S4(TF_Reader_for_Type3(j),:) == 1) , find(flag_S4(AF_Reader(k)) == 1) ) , left_reader);
            if isempty(possible_RC_reader);continue;end % û�ҵ�ֱ������
            % �ҵ�ǰ���Q����
            Q = zeros(1,length(possible_RC_reader));
            for kk = 1:length(possible_RC_reader)
                Q(kk) = min([Pr_ReaderwithReader(TF_Reader_for_Type3(j),possible_RC_reader(kk)),...
                            Pr_ReaderwithReader(possible_RC_reader(kk),AF_Reader(k)), ...
                            Pr_ReaderwithAppterminal(AF_Reader(k))]);
            end
            % ��¼��ǰAF�ľֲ�����
            local_best_RC = possible_RC_reader(Q == max(Q));
            local_best_link3 = [local_best_link3; [TF_Reader_for_Type3(j),local_best_RC(1),AF_Reader(k)]];
            local_best_Q = [local_best_Q; max(Q)];

        end
        % ÿһ�м�¼����һ��TF������·��
        Temp_link3 = [Temp_link3;local_best_link3(local_best_Q == max(local_best_Q),:)];
    end
    Link_Type_3{i} = Temp_link3;

    all_link1_num = all_link1_num + size(Temp_link1,1);
    all_link2_num = all_link2_num + size(Temp_link2,1);
    all_link3_num = all_link3_num + size(Temp_link3,1);
    all_single_tag_link_num(i) = size(Temp_link1,1) + size(Temp_link2,1) + size(Temp_link3,1);
    all_GDOP_raeder{i} = [Temp_link1(:,1) ; Temp_link2(:,1) ; Temp_link3(:,1)];
    all_forward_reader = union(all_forward_reader,[Temp_link2(:,1)', Temp_link3(:,1)', Temp_link3(:,2)']);
    all_uniq_reader = union(all_uniq_reader,[Temp_link1' reshape(Temp_link2,1,[]),reshape(Temp_link3,1,[])]);
end

all_uniq_reader = all_uniq_reader'; % ���Ķ�����ŷ�������

% ����f1
Cn = all_single_tag_link_num >= 3;
f1 = 1 - 1./Num_t.*sum(Cn);

% ����f2 GDOP
f2 = fun_aver_GDOP(all_GDOP_raeder,pos_reader,pos_tag,Cn,all_param);
if isempty(f2); f2 = 1000; end

% ����f3
k1 = 0.5;
k2 = 0.5;
I_T = sum( sum((1-flag_S1).*Pr_TagFromReader,2)./sum(Pr_TagFromReader,2) )./Num_t;
I_R = sum( sum((1-flag_S2).*Pr_ReaderFromTag,2)./sum(Pr_ReaderFromTag,2) )./Num_t;
f3 = k1.*I_T + k2.*I_R;

% ����f4
f4 = (2.*all_link2_num + 3.*all_link3_num)./(all_link1_num + 2.*all_link2_num + 3.*all_link3_num);
if isnan(f4); f4 = 1; end

% ���㹱�׶�
alpha_m(all_uniq_reader) = 1;



% Ȩ��
w1 =1800;
w2 = 1;
w3 = 360;
w4 = 225;%f4Ȩ��0.6��0.8


F=w1*f1+w2*f2+w3*f3+w4*f4;
end