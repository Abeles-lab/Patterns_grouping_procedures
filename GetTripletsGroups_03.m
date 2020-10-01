function NoList = GetTripletsGroups_03(Details,splay, minRepeats, strtI, endI)
% find groups yhat repeated witthin bin at least minRepeat
%  V02 - different clustering algorithms
%  V03 - more efficient code by Ohad
%
% Details    _ structure like DetailsOIsimple
% splay      - intervals may vary in +/-splay samples
% minRepeats - do not keep triplets that repeated in less trials
%
% NoList - <1xN> cell array.  In celln} is the list of items in Details
%          that belong to one group

% Jan-2020  MA

%% internal params
maxExpected = 30000000;   % Max triplets with similar delays
numReport = 10000;         % print when somany triplets were found
maxIter = 3000000;        % stop when somany triplets from Details analyzed
iiReport = 20000;         % print when ii is there

%% initialize
if ~exist('strtI' , 'var'), strtI =[]; end
if isempty(strtI), strtI = 1; end
maxIter = strtI+maxIter -1;
splay2 = 2*splay;
sqrt2splay = sqrt(2)*splay;
No = Details.No;
[rL, rId] = findRuns(No);
cached_mat = thisID_cache_mat_OF( rL, rId );
numRuns = sum(rL>=minRepeats);
if ~exist('endI' , 'var'), endI =[]; end
if isempty(endI), endI = numRuns; end
if endI>numRuns
    endI = numRuns;
end

[SrL,sI] = sort(rL,'Descend'); % sort the runs in decreasing order
SrId = rId(sI);

Bin = ones(2*splay+1,2*splay+1);
NoList = cell(1,maxExpected);
kk = 1;
warning('off')

%% for each run search for groups
tic
for ii = strtI:endI
    thisRun = SrL(ii);
    thisId  = SrId(ii);
    N = cached_mat(thisId,2):cached_mat(thisId,3)';
%     N = find(No==thisId);
    T = Details.TimesInFile(N,:);
    T1 = T(:,2)-T(:,1);
    T2 = T(:,3)-T(:,2);
    T3 = T(:,1)-T(:,3);
    if thisRun==minRepeats   % test directly
        allOK = testOK(T1,T2,T3, splay);
        if allOK
            NoList{kk} = N;
            kk = kk+1;
        end
    elseif thisRun== minRepeats+1 % test directly
        % test if all whithin splay
        allOK = testOK(T1,T2,T3, splay);
        if allOK
            NoList{kk} = N;
            kk = kk+1;
        else
            for jj = 1:minRepeats+1
                elOne = true(1,minRepeats+1);
                elOne(jj) = false; % eliminate one
                Tt1 = T1(elOne);
                Tt2 = T2(elOne);
                Tt3 = T3(elOne);
                allOK = testOK(Tt1,Tt2,Tt3, splay);
                if allOK
                    NoList{kk} = N;
                    kk = kk+1;
                end
            end
        end  %  end of analyzinf min or min+1 repeats
    else                           % use clustering
        [X,Y] = EinthovenCoord(T1,T2);
        X = double(X);
        Y = double(Y);
        L = length(X);
        % eliminate points with no neighbors
        DD = DcityBlock(X,Y);
        toDelete = false(1,L);
        for nn = 1:L
            s2= sum(DD(nn,:)<splay2);
            if s2<minRepeats -1;
                toDelete(nn) = true;
            end
        end
        if sum(~toDelete)>=minRepeats    % try to cluster
            Xin = X(~toDelete);
            Yin = Y(~toDelete);
            [M,~,binX,binY] = OrderRandXY(Xin,Yin, 1);
            MC = myConv2Ddilute(M, Bin);
            MP = findpeaks2D(MC,minRepeats);
            MPmark = zeros(size(MP));
            MPmark(MP>=minRepeats)=1;
            % search the center peaks and find the n-tuplets around them
            [~,rStrt] = find(MPmark'==1,1);
            [~,rEnd] = find(MPmark'==1,1,'Last');
            for jj = rStrt:rEnd
                Raw = MPmark(jj,:);
                OKraw = false(size(Raw));
                [rrL, rrId] = findRuns(Raw);
                rPos = cumsum(rrL) +1;       % position of end of run
                % rPos = cumsum(rrL);
                for rr = 1:length(rrL)
                    if rrId(rr) ==1
                        if rrL(rr)==1        % only one point
                            r1 = rPos(rr)-1; % ??-1??
                            OKraw(r1) = true;
                        else                 % several points
                            r1 = rPos(rr)-1;   % end of this run in Raw
                            r0 = r1-rrL(rr)+1;   % where run starts inRaw
                            mid = round((r0+r1)/2);
                            OKraw(mid) = true;
                        end
                    end  % a peak was found
                end  % end of treating plataus in one raw
                MPmark(jj,:) = OKraw;
            end  % end of going through all raws of MP
            
            % work along the columns
            [~,rStrt] = find(MPmark==1,1);
            [~,rEnd] = find(MPmark==1,1,'Last');
            for jj = rStrt:rEnd
                Raw = MPmark(:,jj);
                OKraw = false(size(Raw));
                [rrL, rrId] = findRuns(Raw);
                rPos = cumsum(rrL) +1;       % position of end of run
                
                for rr = 1:length(rrL)
                    if rrId(rr) ==1
                        if rrL(rr)==1        % only one point
                            r1 = rPos(rr)-1; % ??-1??
                            OKraw(r1) = true;
                        else                 % several points
                            r1 = rPos(rr)-1;   % end of this run in Raw
                            r0 = r1-rrL(rr)+1;   % where run starts inRaw
                            mid = round((r0+r1)/2);
                            OKraw(mid) = true;
                        end
                    end  % a peak was found
                end  % end of treating plataus in one raw
                MPmark(:,jj) = OKraw;
            end  % end of going through all columns of MP
            
            % for neighboring peaks leave only the larger ons
            [Ip,Jp] = find(MPmark);
            PD = NaN(length(Ip),length(Jp));
            for nn = 1:length(Ip)-1
                i0 = Ip(nn);
                j0 = Jp(nn);
                for mm = nn+1:length(Jp)
                    i1 = Ip(mm);
                    j1 = Jp(mm);
                    PD(nn,mm) = abs(i0-i1)+abs(j0-j1);
                    PD(mm,nn) = PD(nn,mm);
                end
            end
            eliminatePeak = false(1,length(Ip));
            for nn = 1:length(Ip)-1
                i0 = Ip(nn);
                j0 = Jp(nn);
                p0 = MP(i0,j0);
                pd = PD(nn,:);
                K0 = find(pd<=sqrt2splay);
                if ~isempty(K0)
                    P1 = zeros(1,length(K0));
                    for mm = 1:length(K0)
                        k0 = K0(mm);
                        P1(mm) = MP(Ip(k0),Jp(k0));
                    end
                    if all(P1<p0)  % The nn-th peak is the largest
                        toEliminate = K0;
                    else           % find the largest
                        %                     toEliminate = nn;
                        Imax = find(P1==max(P1),1);
                        K0(Imax) = [];
                        toEliminate = [nn,K0];
                    end
                    eliminatePeak(toEliminate) = true;
                    for mm = 1:length(toEliminate)
                        k0 = toEliminate(mm);
                        PD(k0,:) = NaN;
                        PD(:,k0) = NaN;
                    end   %% end of eliminating smallew=r peaks
                end   %  end of some neighbors found
            end   %  end of going over all peaks and finding near by
            % delete the surplouse
            eliminate = find(eliminatePeak);
            for nn = 1:sum(eliminatePeak)
                mm = eliminate(nn);
                MPmark(Ip(mm),Jp(mm)) = 0;
            end
            
            % delete data which is not within splay of the centroids
            Mask = false(size(M));
            n = size(Mask,1);
            [I,J] = find(MPmark);  % check Up-Down raws,columns
            for jj = 1:length(I)
                i0 = I(jj)-splay;  % position on Y axis
                if i0<1, i0=1; end
                i1 = I(jj) + splay;
                if i1>n, i1=n; end
                j0 = J(jj)-splay;  % position on X axis
                if j0<1, j0=1; end
                j1 = J(jj) + splay;
                if j1>n, j1=n; end
                % Mask(j0:j1,i0:i1) = true;
                Mask(i0:i1,j0:j1) = true;
            end
            Mpruned = zeros(size(M));
            MleftOut = Mpruned;
            for mm = 1:size(M,1)
                Ip = find(Mask(mm,:));
                Mpruned(mm,Ip) = M(mm,Ip);
                Io = find(~Mask(mm,:));
                MleftOut(mm,Io) = M(mm,Io);
            end
            
            % eliminate the left-out points from X and Y
            % toDelete = false(size(X));
            [MM,NN] = find(MleftOut);
            % NOTE: N is the bin on the X-axis and M on the Y-axis
            for bb = 1:length(MM)
                mm = MM(bb);
                nn = NN(bb);
                m0 = find(binY(mm)<=Y & Y< binY(mm+1));
                n0 = find(binX(nn)<=X & X< binX(nn+1));
                % find shared indices in both m0 an n0
                if length(m0)>length(n0)
                    % reverse the order
                    tmp = m0;
                    m0 = n0;
                    n0 = tmp;
                end % now m0 is shorter
                InBoth = false(1,length(m0));
                for xx = 1:length(m0)
                    x0 = m0(xx);
                    if any(n0==x0) % a match
                        InBoth(xx) = true;
                    end
                end
                toDelete(m0(InBoth)) = true;
                if sum(InBoth)==0
                    warning('Matlab:evnt_tbl:Nomatch',...
                        'Nothing to delete in %d of MM',bb)
                end
            end
            if sum(~toDelete)>=minRepeats
                Xin = X(~toDelete);
                Yin = Y(~toDelete);
                XY = [Xin,Yin];
                %%% does not work? Mpruned = M(Mask);
                % explanation
                % I indicates the Y-axis and its place is binY(I+1)
                % J ondicates the X-axis and its place is binX(J+1)
                K = sum(sum(MPmark));
                StrtMat = zeros(K,2);
                for mm = 1:length(I);
                    ix = J(mm);
                    iy = I(mm);
                    bx = binX(ix+1);
                    by = binY(iy+1);
                    StrtMat(mm,:) = [bx,by];
                end
                [IDC,C] = kmeans(XY,K, 'start',StrtMat , 'emptyaction','drop',...
                    'Distance','cityblock');
                StrtMat = C;
                DD = DcityBlock(StrtMat(:,1),StrtMat(:,2));
                
                while true
                    % for each cluster find if all are within +/- splay
                    % then retore their place in Details and add to NoList
                    Kold = K;
                    [newCenters,DeletedData, DeletedCenters] = ...
                        mergeCenters(Xin, Yin, IDC, DD, StrtMat, splay, minRepeats);
                    StrtMat = newCenters(~DeletedCenters,:);
                    DD(DeletedCenters,:) = [];
                    DD(:,DeletedCenters) = [];
                    Ileft = find(~toDelete);
                    Ileft(DeletedData) = [];
                    toDelete = true(size(toDelete));
                    toDelete(Ileft) = false;
                    if length(Ileft)<minRepeats
                        break
                    end
                    Xin = X(~toDelete);
                    Yin = Y(~toDelete);
                    XY = [Xin,Yin];
                    % IDC(DeletedData) = [];
                    K = size(StrtMat,1);
                    if K==Kold
                        break
                    end
                    [IDC,StrtMat] = kmeans(XY,K, 'start',StrtMat , 'emptyaction','drop',...
                        'Distance','cityblock');
                end       % end of uniting centers
                %  once more
                K = size(StrtMat,1);
                if K>1
                    [IDC,~] = kmeans(XY,K, 'start',StrtMat , 'emptyaction','drop',...
                        'Distance','cityblock');
                end
                % purge centers with less than minRepeats and save the others
                newCenters = unique(IDC);
                RemoveCenters = false(1,length(newCenters));
                RemoveData = false(1,length(IDC));
                for nn = 1:length(newCenters)
                    InHere = find(IDC==newCenters(nn));
                    if length(InHere)<minRepeats
                        RemoveData(InHere) = true;
                        RemoveCenters(nn)  = true;
                    else
                        Group = zeros(1,20*length(InHere));
                        nn0 = 1;
                        for mm = 1:length(InHere)
                            x0 = Xin(InHere(mm));
                            y0 = Yin(InHere(mm));
                            Iadd = find(X==x0 & Y == y0);
                            nn1 = nn0 +length(Iadd)-1;
                            Group(nn0:nn1) = Iadd;
                            nn0 = nn1 +1;
                        end
                        Group = unique(Group);
                        Group(Group==0) = [];
                        NoList{kk} = N(Group);
                        kk = kk+1;
                        if mod(kk,numReport)==0
%                             toc
                            disp (['Reached ' num2str(kk) ' triplets'...
                                ' and ii=' num2str(ii),'out of ', num2str(endI),' at ',datestr(now,'HH.MM.SS_dd.mm.yyyy')])
%                             tic
                        end
                    end
                end
            end        % end of more than minRepeats
        end   % at least minrepets left
        if mod(ii,iiReport)==0
%             toc
            disp (['Reached ' num2str(kk) ' triplets'...
                ' and ii=' num2str(ii)])
%             tic
        end  %  end of running over one run
    end
    if ii==maxIter
        break
    end
end  % end of runing over all runs

%% wrap up
toc
NoList(kk:end)= [];
warning('ON')
return
end

%%%%%%%%%%%%%%%%%   Internal procedures    %%%%%%%%%%%%%%%%%%%%%%
%% allOK
function allOK = testOK(T1,T2,T3, splay)
%% initialise
n = length(T1);
allOK = false;

%% test for T1
% DT = zeros(n,n);
for ii = 1:n-1
    t1 = T1(ii);
    for jj = ii+1:n
        t2 = T1(jj);
        if abs(t2-t1)>splay
            return
        end
%         DT(ii,jj) = abs(t1-t2);
    end
end
% if any(any(DT>splay))
%     return
% end
% 
%% test for T2
% DT = zeros(n,n);
for ii = 1:n-1
    t1 = T2(ii);
    for jj = ii+1:n
        t2 = T2(jj);
        if abs(t2-t1)>splay
            return
        end
%         DT(ii,jj) = abs(t1-t2);
    end
end
% if any(any(DT>splay))
%     return
% end

%% test for T3
% DT = zeros(n,n);
for ii = 1:n-1
    t1 = T3(ii);
    for jj = ii+1:n
        t2 = T3(jj);
        if abs(t2-t1)>splay
            return
        end
%         DT(ii,jj) = abs(t1-t2);
    end
end
% if any(any(DT>splay))
%     return
% end

allOK = true;
return
end

%% DcityBlock
function DD = DcityBlock(X,Y)
% all to all city-block distances
L = length(X);
DD = NaN(L,L);
for ii = 1:L-1
    xy0 = [X(ii),Y(ii)];
    for jj = ii+1:L
        xy1 = [X(jj),Y(jj)];
        dx = xy1-xy0;
        DD(ii,jj) = abs(dx(1)) + abs(dx(2));
        DD(jj,ii) = DD(ii,jj);
    end
end

return
end

%% DwithinLimit
function DD = DwithinLimit(X,Y,splay)
% all to all city-block distances
L = length(X);
DD = true(L,L);
for ii = 1:L-1
    xy0 = [X(ii),Y(ii)];
    for jj = ii+1:L
        xy1 = [X(jj),Y(jj)];
        dx = abs(xy1-xy0);
        DD(ii,jj) = all(dx<=splay);
        DD(jj,ii) = DD(ii,jj);
    end
end

return
end


%%  mergeCenters
function [newCenters,DeleteData, DeleteCenters, IDC] = mergeCenters(Xin, Yin, IDC, DD, StrtMat, splay, minRepeats)
% test if by merging a pair of neighboring centers we get lerger groups

%% initialize
K = size(StrtMat,1);
DeleteCenters = false(1,K);
DeleteData = false(1,length(Xin));
newCenters = StrtMat;
splay2 = 2*splay;

%% go over all pairs
for cc = 1:K-1
    InHere = find(IDC==cc);
    if ~isempty(InHere)
        % n0 = length(InHere);
        xx = Xin(InHere);
        yy = Yin(InHere);
        [Stay0, Remove0] = LintOutliers (xx,yy, splay2);
        nn0 = sum(Stay0);
        % repeat for its next neigbor
        dd = DD(cc,:);
        NextNeighbor = find(dd == min(dd),1);
        if isempty(NextNeighbor)    % rearrange InHere
            Accept = InHere(Stay0);
            Reject = InHere(Remove0);
            IDC(Reject) = NaN;
            DeleteData(Reject) = true;
            xx = Xin(Accept);
            yy = Yin(Accept);
            newCenters(cc,:) = [mean(xx), mean(yy)];
            DD(cc,:) = NaN;
            DD(:,cc) = NaN;
        else
        InThere = find(IDC==NextNeighbor);
        % n1 = length(InThere);
        xx = Xin(InThere);
        yy = Yin(InThere);
        [Stay1, ~] = LintOutliers (xx,yy, splay2);
        nn1 = sum(Stay1);
        allNew = [InThere; InHere];
        % n2 = n0+n1;
        xx = Xin(allNew);
        yy = Yin(allNew);
        [Stay2, ~] = LintOutliers (xx,yy, splay2);
        nn2 = sum(Stay2);
        % distiguish among 8 situations
        if nn0>=minRepeats && nn1>=minRepeats
            if nn1<nn2 && nn0<nn2
                noLost1 = length(InHere)-nn0 + length(InThere)-nn1;
                noLost2 = length(allNew) - nn2;
                if noLost1>=noLost2
                    Case=1;
                else
                    Case=2;
                end
            else
                Case=2;
            end
        elseif nn0>=minRepeats && nn1<minRepeats
            if nn2>nn0
                Case=3;
            else
                Case=4;
            end
        elseif nn0<minRepeats && nn1>=minRepeats
            if nn2>nn1
                Case=5;
            else
                Case=6;
            end
        elseif nn0<minRepeats && nn1<minRepeats
            if nn2>=minRepeats
                Case=7;
            else
                Case=8;
            end
        end
        % switch according to condition.  After this stage, cc may be modified
        % or eliminated all together, its neighbor may be modified.  But, its
        % final membership will be determined only after the do loop reaches it
        switch Case
            case 1  % both 0 and 1 are OK but 2 has more members
                % Remove part of IDC==cc. mark the rest as neighbor.  Remove
                % the center cc. Recompute the center for neighbor
                Accept = unique([InThere; allNew(Stay2)]);
                IDC(allNew) = NaN;
                IDC(Accept) = NextNeighbor;
                DeleteData(allNew) = true;
                DeleteData(Accept) = false;
                xx = Xin(Accept);
                yy = Yin(Accept);
                newCenters(NextNeighbor,:) = [mean(xx), mean(yy)];
                DD(cc,:) = NaN;
                DD(:,cc) = NaN;
            case 2  % both 0 and 1 are ok but 2 has less members
                % remove excess of data with IDC==cc. recompute the center of
                % cc.
                Accept = InHere(Stay0);
                Reject = InHere(Remove0);
                IDC(Reject) = NaN;
                DeleteData(Reject) = true;
                xx = Xin(Accept);
                yy = Yin(Accept);
                newCenters(cc,:) = [mean(xx), mean(yy)];
                DD(cc,:) = NaN;
                DD(:,cc) = NaN;
            case 3  % 0 is OK and 1 is not and 2 has more members than 0
                % remove part of data with IDC==cc mark the other as neighbor.
                % Remove center cc.  recompute the center for neighbor
                Accept = unique([InThere; allNew(Stay2)]);
                IDC(allNew) = NaN;
                IDC(Accept) = NextNeighbor;
                DeleteData(allNew) = true;
                DeleteData(Accept) = false;
                xx = Xin(Accept);
                yy = Yin(Accept);
                newCenters(NextNeighbor,:) = [mean(xx), mean(yy)];
                DD(cc,:) = NaN;
                DD(:,cc) = NaN;
            case 4  % 0 is OK and 1 is not and 2 has not more than 0
                % Remove ecees of data with IDC==cc (if any excess)
                Accept = InHere(Stay0);
                Reject = InHere(Remove0);
                IDC(Reject) = NaN;
                DeleteData(Reject) = true;
                xx = Xin(Accept);
                yy = Yin(Accept);
                newCenters(cc,:) = [mean(xx), mean(yy)];
                DD(cc,:) = NaN;
                DD(:,cc) = NaN;
            case 5  % 1 is OK and 0 is not and 2 has more than 1
                % remove part of data for IDC==cc, mark the rest as neighbor.
                % Recompute the center for neighbor
                IDC(allNew) = NaN;
                Accept = unique([InThere; allNew(Stay2)]);
                IDC(Accept) = NextNeighbor;
                DeleteData(allNew) = true;
                DeleteData(Accept) = false;
                % new center for NextNeighbor
                xx = Xin(Accept);
                yy = Yin(Accept);
                newCenters(NextNeighbor,:) = [mean(xx), mean(yy)];
                DeleteCenters(cc) = true;
                newCenters(cc,:) = NaN;
                DD(cc,:) = NaN;
                DD(:,cc) = NaN;
            case 6  % 1 is OK and 0 is not and 2 has not more than 1
                % remove data with IDC==cc, remove center CC.
                IDC(InHere) = NaN;
                DeleteData(InHere) = true;
                DeleteCenters(cc) = true;
                newCenters(cc,:) = NaN;
                DD(cc,:) = NaN;
                DD(:,cc) = NaN;
            case 7  % both 0 and 1 are not OK but 2 is OK
                % remove part of data with IDC==cc, mark the rest as neighbot
                % remove center cc, recompute the center of neighbor
                IDC(allNew) = NaN;
                Accept = unique([InThere; allNew(Stay2)]);
                IDC(Accept) = NextNeighbor;
                DeleteData(allNew) = true;
                DeleteData(Accept) = false;
                % new center for NextNeighbor
                xx = Xin(Accept);
                yy = Yin(Accept);
                newCenters(NextNeighbor,:) = [mean(xx), mean(yy)];
                DeleteCenters(cc) = true;
                newCenters(cc,:) = NaN;
                DD(cc,:) = NaN;
                DD(:,cc) = NaN;
            case 8  % both 0, 1, and 2 are not OK
                % remove data of IDC==cc, removecenter cc
                DeleteData(InHere) = true;
                DeleteCenters(cc) = true;
                newCenters(cc,:) = NaN;
                IDC(InHere) = NaN;
                DD(cc,:) = NaN;
                DD(:,cc) = NaN;
        end
        end  % end of not empty NextNeighbor
    end
end

return
end

%%  Lint outliers
function [Stay, Remove] = LintOutliers (X,Y, maxD)
% comput allto all distance (city-block) and remove outliers

numData = length(X);
Stay = true(1,numData);
Remove = ~Stay;
DD = DwithinLimit(X,Y,maxD);
if all(DD)
    return
end
I = 1:numData;
DD = ~DD;     % out of limits

%% Itteratively find the outliers until OK
for jj = 1:numData
    inHere = I(~isnan(I));
    numHere = length(inHere);
    % find who is the fartest and remove it
    K = zeros(1,numData);
    for ii = 1:numHere
        K(ii) = sum(DD(ii,:));
    end
    if all(K==0)
        break;
    end
    i0 = find(K == max(K),1);
    I(i0) = [];
    toRemove = inHere(i0);
    Remove(toRemove) = true;
    DD(i0,:) = [];
    DD(:,i0) = [];
end

Stay = ~Remove;
return
end
