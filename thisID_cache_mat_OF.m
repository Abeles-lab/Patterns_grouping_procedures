function [ cached_mat ] = thisID_cache_mat_OF( rL, rId )
%THISID_CACHE_MAT_OF Summary of this function goes here
%   rL, rId are the outputs of findRuns
end_ind = cumsum(rL);
star_ind = [0,end_ind(1:end-1)]+1;
cached_mat = zeros(max(rId),3);
cached_mat(rId,1) = rId;
cached_mat(rId,2) = star_ind;
cached_mat(rId,3) = end_ind;
end

