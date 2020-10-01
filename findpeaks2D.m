function MP = findpeaks2D(M, minHeight)
% find positive peaks whose height is not below minHeight. All points on a
% platau are considered peaks
%    



% May-2017  MA

%% initialize
[numRows, numCols] = size(M);
MP = zeros(numRows, numCols);

% extend M by zeros
ME = zeros(numRows+2, numCols+2);
ME(2:end-1,2:end-1) = M;

% find where high enough
[I,J] = find(ME>=minHeight);


%% search for local maxima
for ii = 1:length(I)
    i0 = I(ii);
    j0 = J(ii);
    vHere = ME(i0,j0);
    vicinity = ME(i0-1:i0+1 , j0-1:j0+1);
    if all(all(vicinity<=vHere))  % may be a peak
        peakOK = true;
        ii0 = i0;
        while peakOK % test to the left
            ii0 = ii0-1;
            v = ME(ii0,j0);
            if v>vHere
                peakOK = false;
                break
            elseif v<vHere
                break
            end
        end
        
        ii0 = i0;
        while peakOK % test to the left
            ii0 = ii0+1;
            v = ME(ii0,j0);
            if v>vHere
                peakOK = false;
                break
            elseif v<vHere
                break
            end
        end
        
        jj0 = j0;
        while peakOK % test to the left
            jj0 = jj0-1;
            v = ME(i0,jj0);
            if v>vHere
                peakOK = false;
                break
            elseif v<vHere
                break
            end
        end
        
        jj0 = j0;
        while peakOK % test to the left
            jj0 = jj0+1;
            v = ME(i0,jj0);
            if v>vHere
                peakOK = false;
                break
            elseif v<vHere
                break
            end
        end
        
        if peakOK
            MP(i0-1, j0-1) = vHere;
            % MP(i0, j0) = vHere;
        end
    end
end

return
