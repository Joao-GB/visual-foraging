function confirmed = confirmationScreen(dpP, drP, prm, keys)

    options  = {'Sim', 'Cancelar'};
    selected = 2;
    blueRect = CenterRectOnPointd([0,0,1.2*dpP.winCenter(1),1.2*dpP.winCenter(2)], dpP.winCenter(1),dpP.winCenter(2));

    %% Posições e tamanhos
    msgY1     = dpP.winCenter(2) * 0.75;
    msgY2     = dpP.winCenter(2) * 0.85;
    btnY      = dpP.winCenter(2) * 1.35;
    btnW      = dpP.winRect(3) * 0.15;
    btnH      = dpP.winRect(4) * 0.08;
    btnGap    = dpP.winRect(3) * 0.1;

    cX = dpP.winCenter(1);

    % Retângulos simétricos em relação ao eixo vertical
    btnRects(:, :, 1) = CenterRectOnPoint([0 0 btnW btnH], cX - btnGap, btnY);
    btnRects(:, :, 2) = CenterRectOnPoint([0 0 btnW btnH], cX + btnGap, btnY);

    KbReleaseWait;
    while true
        
        Screen('FillRect', dpP.window, drP.blue, blueRect);
        Screen('TextSize', dpP.window, prm.textSizeBigger);
        DrawFormattedText(dpP.window, 'Tem certeza que deseja sair?\n', 'center', msgY1, drP.black);
        Screen('TextSize', dpP.window, prm.textSizeBig);
        DrawFormattedText(dpP.window, 'Essa decisão não poderá ser desfeita', 'center', msgY2, drP.blackGrey);

        for i = 1:2
            if i == selected
                frameColor = drP.red;
                frameWidth = 4;
            else
                frameColor = drP.black;
                frameWidth = 2;
            end

            Screen('FrameRect', dpP.window, frameColor, btnRects(:,:,i), frameWidth);

            DrawFormattedText(dpP.window, upper(options{i}), ...
                'center', 'center', drP.black, [], [], [], [], [], btnRects(:,:,i));
        end

        Screen('Flip', dpP.window);

        % Teclado
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(keys{1})     % Seta para esquerda
                selected = 1;

            elseif keyCode(keys{2}) % Seta para direita
                selected = 2;

            elseif keyCode(keys{3}) % Espaço (confirmar)
                confirmed = strcmp(options{selected}, 'Sim');
                Screen('TextSize', dpP.window, prm.textSizeNormal);
                KbReleaseWait;
                return;
            end
            KbReleaseWait;
        end
    end
end