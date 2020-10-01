function [M,Source,binX,binY] = OrderRandXY(X,Y, d)
% reorder X,Y into matrix with resolution d
%    [M,Source,binX,binY] = OrderRandXY(X,Y, d);
%
% X - array with N x-coordinates
% Y - array with N y-coordinates
% d - resolution for placing (X,Y) into matrix
%
% M
% Source - the indices where point ii of(X,Y) is in M
% binX   - the division of X into bins
% binY   - the division of Y into bins

% Jan-2012  MA
%    Updates
% Feb-2020 range of binX, binY extended  MA

%% Intialize
X = double(X);
Y = double(Y);
d = double(d);
minX = min(X);
maxX = max(X);
minY = min(Y);
maxY = max(Y);
n = ceil((maxX-minX)/d) +20;
m = ceil((maxY-minY)/d) +20;
n = max([m,n]);
M = zeros(n,n);
Source = zeros(length(X),2);
x0 = minX-d/2;
y0 = minY-d/2;
% boundaries on X
binX = x0:d:maxX+100*d;
binY = y0:d:maxY+100*d;

%% fill in the matrix
for ii = 1:length(binX)-1
    x1 = binX(ii);
    x2 = binX(ii+1);
    I = X>=x1 & X<x2;
    YX = Y(I);
%     [n,~] = hist(YX,binY);
%     J = find(n>0);
    II = find(I);
    for jj = 1:length(II)
        ij = II(jj);
        j0 = find(YX(jj)>=binY,1,'Last');
        M(j0,ii) = M(j0,ii) +1;
        Source(ij,:) = [j0,ii];
    end
end

return

