AssertOpenGL;

PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1)

screenNumber = max(Screen('Screens'));
bg = [1 1 1]/2;

[win, rect] = PsychImaging('OpenWindow', screenNumber, bg);
[xcenter, ycenter] = RectCenter(rect);
HideCursor;

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
try

    while ~KbCheck

        [mx, my, ~] = GetMouse(win);

        Screen('FillRect', win, bg);
        DrawCalibrationTarget(win, xcenter, ycenter);

        DrawEye(win, [mx my], 15, [50 140 255]/255);

        DrawFormattedText(win, 'Move the mouse. Press any key to quit.', ...
            20, 20, [255 255 255]);

        Screen('Flip', win);

    end

catch ME
    sca;
    ShowCursor;
    rethrow(ME);
end

sca;
ShowCursor;