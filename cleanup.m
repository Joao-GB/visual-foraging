function cleanup(win)
    if nargin < 1, win = 0; end
    if Eyelink('IsConnected')
        Eyelink('Command', 'clear_screen 0');
        Eyelink('StopRecording');
        WaitSecs(0.1);
        Eyelink('Shutdown');
    end
    sca;
    clear Screen;
    clear PsychImaging;
    clear mex;
    
    ShowCursor();
    ListenChar();
    Priority(0);
    close all;
end