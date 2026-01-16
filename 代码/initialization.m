function Positions=initialization(SearchAgents_no,dim)

X=[3,11,19,28];
Y=[6,12,18,24];
count=1;
for i=1:size(X,2)
    for j=1:size(Y,2)
        uniform_point(count,:)=[X(i),Y(j)];
        count=count+1;
    end
end
positions=[];
for i=1:dim/2
    RA=rand(1,SearchAgents_no);
    RB=rand(1,SearchAgents_no);
    position1=uniform_point(i,1)+RA*3.*cos(2*pi*RB);
    position2=uniform_point(i,2)+RA*3.*sin(2*pi*RB);
    positions=[positions;position1;position2];
end
Positions=positions;

end