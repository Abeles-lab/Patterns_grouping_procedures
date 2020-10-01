% base script for grouping triplets
sub_nums=[102, 114, 115, 117, 118, 119, 120];
splay_opts = [2, 3, 5, 8, 10, 12, 15];
minRepeats = 4;
%base_path =  place here the working directory without subjects number for example '/home/subject';
%file_name =  place here the name of the input mat file with orderless triplets.

[~, machineName] = system('hostname');
for sub_ind = 1:length(sub_nums)
    sub_num = sub_nums(sub_ind);
    sub_path = [base_path,num2str(sub_num),'/' ,file_name];
    % load relevant orderTriplets file
    load(sub_path);
    if ~exist('DetailsOIgrouped','var')
        DetailsOIgrouped = DetailsBTgrouped;
    end
    for splay_ind = 1:length(splay_opts)
        splay = splay_opts(splay_ind);
        output_file_name = [base_path,num2str(sub_num),'/262_GroupedDetails_maxLen_100_temporalAccuracy_',num2str(splay*2+1)]
        fid = fopen( [output_file_name,'.txt'], 'wt' );
        fprintf( fid, 'Processing is done on machine: %s', machineName);
        fprintf( fid, 'Start time is %s', datestr(now,'HH:MM:SS_dd.mm.yyyy'));
        
        disp(['Starting run ',datestr(now,'HH:MM:SS_dd.mm.yyyy')]);
        NoList = GetTripletsGroups_03(DetailsOIgrouped,splay, minRepeats);
        disp(['Finishing run ',datestr(now,'HH:MM:SS_dd.mm.yyyy')]);
        % Note:  This takes a lot of time.  You may run it in parts like:
        % strtI = 1;
        % endI = 100001;
        % NoList1 = GetTripletsGroups_03(DetailsOIgrouped,splay,minRepeats, strtI, endI);
        % strtI = 100001;
        % endI = 200000;
        % NoList2 = GetTripletsGroups_03(DetailsOIgrouped,splay,minRepeats, strtI, endI);
        % and then merge them in pairs like
%         NoList = mergeNoLists(NoList1,NoList2);
%         save NoListTmp NoList
        Lint=0;
        disp(['Starting NoList2details at ',datestr(now,'HH:MM:SS_dd.mm.yyyy')]);
        DetailsBTgrouped = NoList2details(DetailsOIgrouped, NoList, Lint);
        save([output_file_name,'.mat'],'DetailsBTgrouped','Lint','splay','minRepeats','-v7.3')
        disp([output_file_name,'.mat saved at ',datestr(now,'HH:MM:SS_dd.mm.yyyy')]);
        clear DetailsBTgrouped;
    end
    
end
