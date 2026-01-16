function result = randpermF( range, dim )
% The original function "randpermD" in Matlab is only confined to the
% situation that dimension is no bigger than dim. This function is 
% applied to solve that situation.
%Matlab 中原有的 randpermD 函数存在维度限制，只能处理维度不大于某个特定值（可能与这里的 dim 相关）的情况，
%而 randpermF 函数就是为了突破这个限制，实现在更宽泛维度要求下生成特定随机排列的功能。

temp = randpermD( range, range );
temp2 = randi( range, dim, 1 );
index = randpermD( dim, ( dim - range ) );
result = [ temp, temp2( index )' ];
%相当于构建了一部分满足随机排列要求的基础数据