function foragingDrawPedestal(auxWin, noiseTex, srcRects, dstRects, ori, txP)

    dstRects = circshift(dstRects, -1, 2);
    ori = circshift(ori, 1, 2);

    Screen('BlendFunction', auxWin, GL_ONE, GL_ONE);

    Screen('DrawTextures', auxWin, noiseTex, srcRects, dstRects, ori);

    Screen('BlendFunction', auxWin, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
    Screen('DrawTextures', auxWin, txP.blob.tex, [], dstRects, ori, [], [], [0 0 0 1]', [], [], txP.blob.props);

    Screen('BlendFunction', auxWin, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
end