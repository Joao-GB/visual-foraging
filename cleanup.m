function cleanup(win)
    if nargin < 1, win = 0; end
    if Eyelink('IsConnected')
        Eyelink('StopRecording');
        WaitSecs(0.1);
        Eyelink('Shutdown');
    end
    sca;
    ShowCursor([], win)
    ListenChar(0);
    Priority(0);
    close all;
end