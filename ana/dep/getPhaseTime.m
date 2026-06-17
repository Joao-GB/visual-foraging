function phaseTime = getPhaseTime(allTimes, phaseLims, trlLims)

    phaseTime = [find(allTimes == phaseLims(1)) - find(allTimes == trlLims(1)) ...
                 find(allTimes == phaseLims(2)) - find(allTimes == trlLims(1))];
end