function endScreen(dpP, drP, prm)
    % Funções e parâmetros particulares à forma do olho que escolhi
    % Variáveis x de posição, s era de escala, agora inútil
    uMaxLid = @(x,s) -(0.203267 + 0.0712797*x - 0.888037*x.^2);
    lMaxLid = @(x,s) -(-0.181055 - 0.0429848*x + 0.753893*x.^2 + 1.02328*x.^3);
    lLim = -.545; rLim = .455; dx = .01;
    X = lLim:dx:rLim;
    
    dt = .005;
    blkDur = .07; totalDur = 100; blkFreq = 1/4;
    [blkTimes, blkShape] = blkProps(dt, blkDur, totalDur, blkFreq);
    blkShape = [0 blkShape];
    lBlk = numel(blkShape);
    
    % Para piscar, t \in [0.2, 1]. Para arregalar, [0, 0.2]
    tRelMin = .2; tRelMax = 1; aRel = .85;
    alpha1Rel = @(t) aRel*(tRelMin*(1-t)+t*tRelMax);
    alpha2Rel = @(t) (1-aRel)*(tRelMin*(1-t)+t*tRelMax);
    uRelLid = @(x,t) alpha1Rel(t)*lMaxLid(x,t) + (1-alpha1Rel(t))*uMaxLid(x,t);
    lRelLid = @(x,t) alpha2Rel(t)*uMaxLid(x,t) + (1-alpha2Rel(t))*lMaxLid(x,t);
    
    % Olhos abertos espremidos usam .5; para piscar, o inferior vai no máximo
    % até .5 e o superior a 1.75
    tuSqnMin = .5; tuSqnMax = 1.75; tlSqn = .5;
    aSqn = .4;
    alpha1uSqn = @(t) aSqn*(tuSqnMin*(1-t)+t*tuSqnMax);
    alpha2lSqn = @(t) (1-aSqn)*tlSqn;
    uSqnLid = @(x,t) alpha1uSqn(t)*lMaxLid(x,t) + (1-alpha1uSqn(t))*uMaxLid(x,t);
    lSqnLid = @(x,t) alpha2lSqn(t)*uMaxLid(x,t) + (1-alpha2lSqn(t))*lMaxLid(x,t);
    
    uLid = @(x,a,t) a*uRelLid(x,t) + (1-a)*uSqnLid(x,t);
    lLid = @(x,a,t) a*lRelLid(x,t) + (1-a)*lSqnLid(x,t);
    
    
    eyeCols = [drP.blue; drP.darkGreen; drP.brown; drP.darkBrown; drP.greyBrown; drP.paleBrown]; aux = randsample(1:size(eyeCols, 1), 1);
    eyeColor = eyeCols(aux, :);
    maskTex = Screen('OpenOffscreenWindow', dpP.window, [0 0 0 0], [], 32);
    xc = dpP.winCenter(1); yc = dpP.winCenter(2);
    
    % Parâmetros importantes
    scale = dpP.winRect(3)/10; 
    eyeWindowCenter  = [3*dpP.winRect(3)/8 dpP.winRect(4)/3;
               5*dpP.winRect(3)/8 dpP.winRect(4)/3]';
    eyeBallOffset = dpP.winRect(3)/2 - 1.01*eyeWindowCenter(1,1);
    
    maxRadius = scale; minRadius = scale/3;
    eyeD = scale*.85; eyeR = eyeD/2;
    baseRect = [0 0 eyeD eyeD];
    
    eyeBallRect = zeros(2,4);
    eyeCenter = [dpP.winRect(3)/2-eyeBallOffset, eyeWindowCenter(2,1), 0;
                 dpP.winRect(3)/2+eyeBallOffset, eyeWindowCenter(2,2), 0];
    eyeBallRect(1,:) = CenterRectOnPointd(baseRect, dpP.winRect(3)/2-eyeBallOffset, eyeWindowCenter(2,1));
    eyeBallRect(2,:) = CenterRectOnPointd(baseRect, dpP.winRect(3)/2+eyeBallOffset, eyeWindowCenter(2,2));
    
    irisD = scale*.4; irisR = irisD/2;
    
    pupilD = scale*.2; pupilR = pupilD/2;
    
    Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    tStart = GetSecs;
    
    nBlkT = numel(blkTimes); blkPlayed = zeros(1, nBlkT);
    yawMax   = deg2rad(30); pitchMax = deg2rad(10); zDist = 500;
    
    eyeState = 0; % 0 para blink, 1 para pursuit, 2 para sacada
    waitTime = 2.5; startWait = -1; minDist1 = 100; T = linspace(0, 2*pi, 100);
    yaw = [0 0]; pitch = [0 0]; wAngles = .9; wT0 = 0; wT1 = .50001; wT2 = .5;
    wT = wT1;
    vLim1 = 200; vLim2 = 1000; vLim3 = 4000; eyeVLim = 5;
    w = .9;
    
    alpha = 0; start = true;
    yRange = dpP.winRect(4)*([-.103 .08] + 1/3); 
    yRange1 = dpP.winRect(4)*([-.0145 .0145] + 1/3);
    i = 1;

    grey = repmat(drP.grey, [1, 3]); white = repmat(drP.white, [1, 3]); black = repmat(drP.black, [1, 3]);
    SetMouse(xc, yc/2, dpP.window); HideCursor(dpP.window);
    while ~KbCheck
        [xMouse, yMouse] = GetMouse(dpP.window);
    %     disp(yMouse > yRange(1) && yMouse < yRange(2))
        tu = blkShape(i);
    
        tNow = mod(GetSecs - tStart, totalDur);
        target = [xMouse yMouse zDist];
    
        if start, lastYaw = yaw; lastPitch = pitch; lastPos = [xMouse; yMouse]; lastTu = tu; start = false; lastVel = 0; nextFix = lastPos; end
    
        vel = vecnorm([xMouse; yMouse] - lastPos)/dt;
    
        lastPos = [xMouse; yMouse];
        vel = w*lastVel+(1-w)*vel;
        lastVel = vel;
    
    
        % Sacada: se alvo estiver se mexendo ou olhos se mexendo durante sacada
        if vel > vLim1 && vel < vLim2 || (eyeState == 2 && any(vecnorm([yaw; pitch] - [lastYaw; lastPitch])/dt > eyeVLim))
            % vecnorm(nextFix - [xMouse; yMouse]) > minDist1; vecnorm(nextFix - currFix) < minDist2
            if vecnorm(nextFix - [xMouse; yMouse]) > minDist1
                nextFix = [xMouse; yMouse];
            end
            eyeState = 2;
            wT = wT2;
            startWait = -1;
        elseif vel > vLim2 && vel < vLim3
            eyeState = 1;
            wT = wT1;
            startWait = -1;
        else
            if startWait == -1
                startWait = tNow;
            end
            if tNow - startWait > waitTime
                eyeState = 0;
                wT = wT0;
                yaw = [0 0]; pitch = [0 0];
                nextFix = [xMouse; yMouse];
            end
        end
    %     disp(vel)
    %     disp(eyeState)
    
        mBlkT = find(blkTimes < tNow, 1, 'last');
        if eyeState == 0 || (eyeState == 1 && i ~= 1)
            if ~isempty(mBlkT)
                if blkPlayed(mBlkT) == 0
                    i = i + 1;
                    if i == lBlk+1
                        i = 1;
                        blkPlayed(:) = 0;
                        blkPlayed(mBlkT) = 1;
                    end
                end
            end
        elseif eyeState >= 1
            
            if eyeState == 2
                target = [nextFix' zDist];
            end
    
            v = target - eyeCenter;
                
            yaw = atan2(v(:,1), v(:,3));
            pitch = atan2(v(:,2), v(:,3));
    
            yaw   = max(min(yaw,   yawMax),   -yawMax);
            pitch = max(min(pitch, pitchMax), -pitchMax);
            if exist('C', "var")
%                 disp(C(2))
                tu = -1+2.3*(min(yRange1(2), max(C(2), yRange1(1)))-yRange1(1))/(yRange1(2)-yRange1(1));
            else
%                 disp(-1)
                tu = -1+2.3*(min(yRange(2), max(yMouse, yRange(1)))-yRange(1))/(yRange(2)-yRange(1));
            end
    
%             disp(min(yRange(2), max(yMouse, yRange(1))))
    
        end
        yaw = wAngles*lastYaw + (1-wAngles)*yaw;
        pitch = wAngles*lastPitch + (1-wAngles)*pitch;
        lastYaw = yaw; lastPitch = pitch;
    
        tl = blkShape(i);
        tu = wT*lastTu +(1-wT)*tu;
    
        dist = hypot(xMouse - eyeWindowCenter(1,:), yMouse - eyeWindowCenter(2,:));
        dist = (min(maxRadius, max(dist, minRadius))-minRadius)./(maxRadius-minRadius);
    
        [rWindow, lWindow] = drawEye(X, dist, tu, tl, uLid, lLid, scale, eyeWindowCenter');
        
        % A máscara não deixa passar nada, a não ser que esteja nos buracos
        Screen('FillRect', maskTex, [grey 0]);
        Screen('FillPoly', maskTex, [grey 1], rWindow);
        Screen('FillPoly', maskTex, [grey 1], lWindow);
    
        Screen('FillOval', dpP.window, drP.white, eyeBallRect(1,:));
        Screen('FillOval', dpP.window, drP.white, eyeBallRect(2,:));
    
                
        [rIrisEl, C] = ellipsePoly(irisR, eyeR, yaw(1), pitch(1), alpha, eyeCenter(1,1:2)', T);
%         disp(C(2))
        rPupilEl= ellipsePoly(pupilR, eyeR, yaw(1), pitch(1), alpha, eyeCenter(1,1:2)', T, 1);
    
        lIrisEl = ellipsePoly(irisR, eyeR, yaw(2), pitch(2), alpha, eyeCenter(2,1:2)', T);
        lPupilEl= ellipsePoly(pupilR, eyeR, yaw(2), pitch(2), alpha, eyeCenter(2,1:2)', T, 1);
        
                
        Screen('FillPoly', dpP.window, eyeColor, rIrisEl'); Screen('FillPoly', dpP.window, eyeColor, lIrisEl');
        Screen('FillPoly', dpP.window, drP.black, rPupilEl'); Screen('FillPoly', dpP.window, drP.black, lPupilEl');
    
        Screen('BlendFunction', dpP.window, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
        Screen('DrawTexture', dpP.window, maskTex);
    
        Screen('BlendFunction', dpP.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        Screen('FramePoly', dpP.window, drP.black, rWindow, 4);
        Screen('FramePoly', dpP.window, drP.black, lWindow, 4);

        Screen('TextSize', dpP.window, prm.textSizeEnormous); Screen('TextStyle', dpP.window, 1);
        DrawFormattedText(dpP.window, 'Sessão concluída!', 'center', 2*dpP.winRect(4)/3, drP.black);
    
        msg = [
            'Aperte qualquer tecla para sair, ou\n'...
            'mova o cursor para ver o que acontece.\n'
        ];
        Screen('TextSize', dpP.window, prm.textSizeBig); Screen('TextStyle', dpP.window, 0);
        DrawFormattedText(dpP.window, msg, 'center', 2*dpP.winRect(4)/3 + 55, drP.black);

        Screen('FillOval',  dpP.window, [white .5], [xMouse-3*prm.cursorRadius_px yMouse-3*prm.cursorRadius_px xMouse+3*prm.cursorRadius_px yMouse+3*prm.cursorRadius_px]);
        Screen('FrameOval', dpP.window, [black .5], [xMouse-3*prm.cursorRadius_px yMouse-3*prm.cursorRadius_px xMouse+3*prm.cursorRadius_px yMouse+3*prm.cursorRadius_px], prm.pW1);
        
        Screen('Flip', dpP.window);
        WaitSecs(dt);
    end
end


function [rEye, lEye] = drawEye(X, d, tu, tl, uFunc, lFunc, s, c)
    fX = fliplr(X);
    nX = -X; fnX = fliplr(-X);
    rEye = [nX fnX;uFunc(X, d(1), tu) lFunc(fX, d(1), tl)]'.*s+c(1,:);
    lEye = [X fX;uFunc(X, d(2), tu) lFunc(fX, d(2), tl)]'.*s  +c(2,:);
end


function [bT, bS] = blkProps(tRes, bDur, totDur, bFreq)
    frac    = .2;  % Fração da piscada que é fechada, frac \in [0,1)
    bRefrac = 1;  % Período refratário da piscada
    bDelay  = 1/bFreq;
    nB      = totDur/bDelay;
    alpha = 2; theta = (bDelay-bRefrac)/alpha;

    bT = cumsum(gamrnd(alpha, theta, [1, nB]));
    bT(bT >= totDur) = [];
    bS = min(1, 1/(1-frac)*(1-abs(1-(2/bDur)*(0:tRes:bDur))));
end


function [e, C] = ellipsePoly(r, R, y, p, al, C, T, mode)
            if nargin < 8, mode = 0; end
            if mode == 1
                a = r* cos(1.3*y) ; b = r * cos(1.3*p);
            else
                a = r* cos(1.1*y) ; b = r * cos(1.1*p);
            end
            e = [a * cos(T); b * sin(T)];
            Rot = [cos(al) -sin(al);
                sin(al)  cos(al)];
            if mode == 1
                C = [C(1,:) + R * sin(1.3*y); C(2,:) + R * sin(1.3*p)];
            else
                C = [C(1,:) + R * sin(1.1*y); C(2,:) + R * sin(1.1*p)];
            end
            e = Rot * e + C;
end