function DetailsBTgrouped = NoList2details(Details, NoList, Lint)
% build the grouped By Time detailsfrom NoList while linting


% Jan 2020  MA

%% initialize
if ~exist('Lint', 'var'), Lint = []; end
if isempty(Lint), Lint = 0; end
switch Lint
    case 0
        removeDups = false;  % true if to remove triplets in which the 
                             % same POI repeats
        removeMulti = false; % true is to remove triplets that happenned
                             % several times in the same trial
    case 1
        removeDups = true;
        removeMulti = false;
    case 2
        removeDups = false;
        removeMulti = true;
    case 3
        removeDups = true;
        removeMulti = true;
    otherwise
        error('Matlab:evnt_tbl:BadParam' ,...
            'Lint must be 0, or 1,.. or 3.')
end
DetailsBTgrouped = struct('Names',0 , 'Behavior',0,...
                 'Trial',0 , 'No',0 , 'TimesInFile',0 , ...
                 'OrgNo', 0);
numTypes = length(NoList);
% get max number of triplets
s = 0;
for ii = 1:numTypes
    s = s + length(NoList{ii});
end
numTriplets = s;
Names = zeros(numTriplets,3, 'int16');
Behavior = zeros(numTriplets,1, 'int16');
Trial = zeros(numTriplets,1, 'int16');
TimesInFile = zeros(numTriplets,3, 'int32');
No = zeros(numTriplets,1, 'int32');
OrgNo = zeros(numTriplets,1, 'int32');

kk = 1;

%% extract the data from Details
for ii = 1:numTypes
    N = NoList{ii};
    Nlong = length(N);
    Nhere = Details.Names(N,:);
    TIFhere = Details.TimesInFile(N,:);
    Bhere = Details.Behavior(N);
    There = Details.Trial(N);
    NOhere = ii*ones(Nlong,1);
    % check if to lint
    toKeep = true;
    if removeDups     % remove triplets withthe same POI twice
        nn = Nhere(1,:); 
        if nn(1)==nn(2) || nn(1)==nn(3) || nn(2)==nn(3)
            toKeep = false;
        end
    end
    if removeMulti    % remove triplets that apear twice in the same trial
        b1 = Bhere==1;
        b2 = Bhere==2;
        T1 = There(b1);
        T2 = There(b2);
        if length(T1)>length(unique(T1))  % some duplicates
            toKeep = false;
        end
        if length(T2)>length(unique(T2))  % some duplicates
            toKeep = false;
        end
    end
    if toKeep
        kk1 = kk + Nlong -1;
        Names(kk:kk1,:) = Nhere;
        Behavior(kk:kk1) = Bhere;
        Trial(kk:kk1) = There;
        TimesInFile(kk:kk1,:) = TIFhere;
        No(kk:kk1) = NOhere;
        OrgNo(kk:kk1) = N;
        kk = kk1+1;
    end  % end of keeping all cases of this triplet
end

%% wrap up
Names(kk1:end,:) = [];
TimesInFile(kk1:end,:) = [];
Behavior(kk1:end) = [];
Trial(kk1:end) = [];
No(kk1:end) = [];
OrgNo(kk1:end) = [];

DetailsBTgrouped.Names  = Names;
DetailsBTgrouped.TimesInFile  = TimesInFile;
DetailsBTgrouped.Behavior  = Behavior;
DetailsBTgrouped.Trial  = Trial;
DetailsBTgrouped.No  = No;
DetailsBTgrouped.OrgNo  = OrgNo;

return
end