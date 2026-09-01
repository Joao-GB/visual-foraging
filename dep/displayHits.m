function displayHits(prm, dpP, drP, nHits, nAns)
% displayHits
%
% Displays the number of correct answers at the center of the screen and
% waits until the Escape key is pressed.

%     Screen('FillRect', dpP.window, drP.grey);

    [screenWidth, screenHeight] = Screen('WindowSize', dpP.window);

    Screen('TextSize', dpP.window, prm.textSizeBigger);

    titleText = 'NÚMERO DE ACERTOS';
    hitsText  = [num2str(nHits) '/' num2str(nAns)];

    % Vertical spacing between the two lines
    lineSpacing = round(prm.textSizeBigger * 0.5);

    titleBounds = Screen('TextBounds', dpP.window, titleText);
    hitsBounds  = Screen('TextBounds', dpP.window, hitsText);

    titleHeight = titleBounds(4) - titleBounds(2);
    hitsHeight  = hitsBounds(4) - hitsBounds(2);

    totalHeight = titleHeight + lineSpacing + hitsHeight;

    titleY = (screenHeight - totalHeight) / 2;
    hitsY  = titleY + titleHeight + lineSpacing;

    DrawFormattedText(dpP.window, titleText, 'center', titleY, drP.black);

    DrawFormattedText(dpP.window,  hitsText, 'center',  hitsY, drP.black);

    Screen('TextSize', dpP.window, prm.textSizeNormalish);

    exitText = 'Pressione Esc para sair';

    bottomMargin = round(prm.textSizeNormalish * 1.5);

    exitBounds = Screen('TextBounds', dpP.window, exitText);
    exitHeight = exitBounds(4) - exitBounds(2);

    exitY = screenHeight - bottomMargin - exitHeight;

    DrawFormattedText( ...
        dpP.window, ...
        exitText, ...
        'center', ...
        exitY, ...
        drP.black);

    Screen('Flip', dpP.window);

    escapeKey = KbName('ESCAPE');
    KbReleaseWait;

    while true
        [keyIsDown, ~, keyCode] = KbCheck;

        if keyIsDown && keyCode(escapeKey)
            break;
        end
    end

    KbReleaseWait;
    Screen('Flip', dpP.window);
end