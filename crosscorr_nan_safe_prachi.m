% helper function for CrossCorrelationsOnCombinedData.m

function corrVals = crosscorr_nan_safe_prachi(x, y, maxLag)
    nLags = 2*maxLag + 1;
    corrVals = nan(nLags, 1);
    lagOffsets = -maxLag:maxLag;

    for i = 1:nLags
        lag = lagOffsets(i);

        if lag < 0
            xLag = x(1:end+lag);
            yLag = y(1-lag:end);
        elseif lag > 0
            xLag = x(1+lag:end);
            yLag = y(1:end-lag);
        else
            xLag = x;
            yLag = y;
        end

        validIdx = ~isnan(xLag) & ~isnan(yLag);
        if sum(validIdx) >= 5
            xz = xLag(validIdx);
            yz = yLag(validIdx);
            corrVals(i) = corr(xz, yz);
        end
    end
end
