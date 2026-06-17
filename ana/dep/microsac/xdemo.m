clear;

% Parameters
MINDUR = 3;        % Minimum duration (number of samples)
VTHRES = 6;        % Velocity threshold
SAMPLING = 250;    % Sampling rate     
VTYPE = 2;         % Velocity types (2 = using moving average)
DPP = 0.0675*60;   % Angular minutes per pixel

% Load raw data
d = load('demo.dat');
xl = DPP*d(:,2:3);
xl(:,1) = xl(:,1) - mean(xl(:,1));
xl(:,2) = xl(:,2) - mean(xl(:,2));
xr = DPP*d(:,4:5);
xr(:,1) = xr(:,1) - mean(xr(:,1));
xr(:,2) = xr(:,2) - mean(xr(:,2));

% Compute 2D velocity vectors
vl = vecvel(xl,SAMPLING,VTYPE);
vr = vecvel(xr,SAMPLING,VTYPE);

% Detection of microsaccades
sacl = microsacc(xl,vl,VTHRES,MINDUR);
sacr = microsacc(xr,vr,VTHRES,MINDUR);

% Testing for binocular saccades via temporal overlap
sac = binsacc(sacl,sacr);

% Plot raw data and binocular saccades
plot(xr(:,1),xr(:,2),'b-');
set(gca,'FontSize',16);
axis square;
hold on;
for s=1:size(sac,1)
    idx = sac(s,1):sac(s,2);
    plot(xr(idx,1),xr(idx,2),'r.-','linewidth',2);
end
hold off;
xlim([-30 30]);
ylim([-30 30]);
set(gca,'xtick',[-30 -15 0 15 30]);
set(gca,'ytick',[-30 -15 0 15 30]);
xlabel('x [min arc]');
ylabel('y [min arc]');
