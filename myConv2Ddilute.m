function MC = myConv2Ddilute(M,Bin)
% convolve a diluted matrix M by a matrix Bin


% Jan-2020  MA

%% initialize
[rawsData, colsData] = size(M);
[rawsBin,colsBin] = size(Bin);
if Iseven(rawsBin)
    BB = [zeros(1,colsBin);Bin];
    Bin = BB;
    [rawsBin,colsBin] = size(Bin);
end
if Iseven(colsBin)
    BB = [zeros(rawsBin,1),Bin];
    Bin = BB;
    [rawsBin,colsBin] = size(Bin);
end
halfRB = ceil(rawsBin/2);
halfCB = ceil(colsBin/2);
mHalfRB = halfRB-1;
mHalfCB = halfCB-1;
% extend M with zeros on all sides
Tmp = zeros(rawsData+2*halfRB, colsData+2*halfCB);
Tmp(halfRB+1:end-halfRB,halfCB+1:end-halfCB) = M;
Rslt = zeros(size(Tmp));

%% convolve
[I,J] = find(Tmp);  % check Up-Down raws,columns
for jj = 1:length(I)
    ix = I(jj);
    i0 = ix-mHalfRB;
    i1 = i0+rawsBin -1;
%     k0 = i0-halfRB;
%     k1 = i1-halfRB;
    jy = J(jj);
    j0 = jy-mHalfCB;
    j1 = j0+colsBin -1;
%     m0 = j0-halfCB;
%     m1 = j1-halfCB;
    D = Tmp(ix,jy)*Bin;
    Rslt(i0:i1,j0:j1) = Rslt(i0:i1,j0:j1) + D;
end

%% restore the center of Rslt
MC = Rslt(halfRB+1:end-halfRB, halfCB+1:end-halfCB);
return

