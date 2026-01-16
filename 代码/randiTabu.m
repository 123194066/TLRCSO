function value = randiTabu( min, max, tabu, dim )
value = ones( dim, 1 ) .* max .* 2;
num = 1;
while ( num <= dim )
    temp = randi( [min, max], 1, 1 );
    if( length( find( value ~= temp ) ) == dim && temp ~= tabu )
        value( num ) = temp;
        num = num + 1;
    end
end
%生成一个长度为 dim 的随机整数向量，其中的每个元素都在 [min, max] 区间内，
%并且每个元素都不能等于给定的禁忌值 tabu，同时各个元素之间互不相同