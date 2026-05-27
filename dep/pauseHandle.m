function [kgBlocks, rBlock, kgTrials, rTrial] = pauseHandle(kgBlocks, rBlock, kgTrials, rTrial, wasRecording, tkP, txP, dpP, drP, prm, pauseMode, debug, foragingMode, ori)
    pauseStart = GetSecs;
    if debug == 0 && foragingMode >= 2
        Eyelink('Message',prm.msg.on.pse);
        Eyelink('StopRecording');
        Eyelink('SetOfflineMode');
    end
    decision = pauseScreen(tkP, dpP, drP, prm, pauseMode);
    
    if strcmp(pauseMode, 'trial')
        if strcmp(decision, 'gotoMenu')
            if debug == 0 && foragingMode >= 2, Eyelink('Message',prm.msg.pse{1}); end
            kgTrials = false;
            kgBlocks = false;
        elseif strcmp(decision, 'restartBlock')
            if debug == 0 && foragingMode >= 2, Eyelink('Message',prm.msg.pse{2}); end
            rBlock = true;
            kgTrials = false;
        elseif strcmp(decision, 'recalibrate')
            if debug == 0 && foragingMode >= 2, Eyelink('Message', prm.msg.pse{3}); end
            if debug == 0 && foragingMode >= 2
                EyelinkDoTrackerSetup(tkP.el);
            else
                disp('Recalibragem solicitada')
            end
            rTrial = true;
        elseif strcmp(decision, 'resume')
            if debug == 0 && foragingMode >= 2, Eyelink('Message',prm.msg.pse{4}); end
            rTrial = true;
        end
        pauseEnd = GetSecs;
        pauseDur = pauseEnd - pauseStart;
        if kgTrials
            if pauseDur > prm.pauseMaxDur
                fakeLoadingScreen(tkP, dpP, drP, prm, 'trial', txP, ori);
                Screen('Flip',dpP.window);
            end
            WaitSecs(.01);
            if wasRecording
                if debug == 0 && foragingMode >= 2, Eyelink('StartRecording'); end
            end
        end
    elseif strcmp(pauseMode, 'block')
        if strcmp(decision, 'gotoMenu')
            if debug == 0 && foragingMode >= 2, Eyelink('Message', prm.msg.pse{1}); end
            kgBlocks = false;
        elseif strcmp(decision, 'recalibrate')
            if debug == 0 && foragingMode >= 2, Eyelink('Message',prm.msg.pse{3}); end
            if debug == 0 && foragingMode >= 2
                EyelinkDoTrackerSetup(tkP.el);
            else
                disp('Recalibragem solicitada')
            end
        elseif strcmp(decision, 'resume')
            if debug == 0 && foragingMode >= 2, Eyelink('Message',prm.msg.pse{4}); end
        end
        if kgBlocks
            if wasRecording
                if debug == 0 && foragingMode >= 2, Eyelink('StartRecording'); end
            end
        end
    end
    if debug == 0 && foragingMode >= 2, Eyelink('Message',prm.msg.off.pse); end
end