function foragingDrawMain(auxWin, gaborTex, noiseTex, srcRects, dstRects, ori, txP, alphas)

    Screen('BlendFunction', auxWin, GL_ONE, GL_ONE);
    Screen('DrawTextures', auxWin, gaborTex, [], dstRects, ori, [], [], alphas);
    Screen('DrawTextures', auxWin, noiseTex, srcRects, dstRects, ori);

% (j) Desenha a abertura gaussiana
    Screen('BlendFunction', auxWin, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
    Screen('DrawTextures', auxWin, txP.blob.tex, [], dstRects, ori, [], [], [0 0 0 1]', [], [], txP.blob.props);

    Screen('BlendFunction', auxWin, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);