clear;clc;
close all;
% load('GZ59_16(0825)删除EMTR中最好的三组数据')
load('C:\Users\Administrator\Desktop\ALL_四个场景\场景1密度2\changjing1midu2-f1.mat')
format short

a=pos_tag;
dimension=32;
for qwe=1:5
    A=[];
    AAA=[];
    if qwe==1
        xxxx=PSO_BestFit;
    end
    if qwe==2
        xxxx=CSO_BestFit;
    end
    if qwe==3
        xxxx=GWO_BestFit;
    end
    if qwe==4
        xxxx=GA_BestFit;
    end
    if qwe==5
%         [~,I]=min(SCSO); %ECSO
        xxxx=MNSCSO_BestFit;
    end
    
    A=reshape(xxxx,2,dimension/2)';
    
    b=A;
    Cn=Link_Table(a,A);
    idx_Cn=[];aaa=[];
    idx_Cn=find(Cn==1);
    aaa=a(idx_Cn,:);
    figure(qwe)
    plot(a(:,1),a(:,2),'o')
    hold on
    plot(aaa(:,1),aaa(:,2),'*')
    hold off
end



