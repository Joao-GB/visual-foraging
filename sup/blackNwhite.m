load('s00_01_i.mat')
minDist_px   = dva2pix(prm.screenDist, dpP.monitorW_mm/10, dpP.screenRes.width, prm.minDist_dva+prm.gaborSize_dva);
AssertOpenGL; Screen('Preference','SkipSyncTests',1);
[win,winRect]=Screen('OpenWindow',max(Screen('Screens')),0); winCenter=winRect(3:4)/2;
[currFixCenter,currStimCenter,~]=getStimLocations2_1(winRect(3:4),[winCenter 1],140,minDist_px,txP.gabor.size_px);
s=txP.gabor.size_px/2;
r=[currStimCenter(1,:)-s;currStimCenter(2,:)-s;currStimCenter(1,:)+s;currStimCenter(2,:)+s];
Screen('FillOval',win,WhiteIndex(win),r); Screen('Flip',win);
for k=1:5, KbStrokeWait; end
img=Screen('GetImage',win); imwrite(img,'stimulus_screenshot.png');
Screen('CloseAll');