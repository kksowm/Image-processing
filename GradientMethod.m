clear
%close all
% Noise level for simulation
bnoise=0.1;
%Threshold to remove noisy pixels
thresholdA =1;
thresholdP =1;

CornealRI = 1.3376;
AqueousRI = 1.334;
ViterousRI = 1.334;
LensRI = 1.406;

% Length over which Anterior Surface is fit
FitRadius = 1.5; % In mm

% size of image from OCT Device caliper settings
%HorzScanLength = 5.1;
%VerticalScanLength = 4.5;
HorzScanLength = 16.5;
VerticalScanLength = 16.5;

% Image size
xlength = 1151;% 1024;
ylength = 372; %512;

xval = linspace(0,HorzScanLength,xlength);
yval = linspace(0,VerticalScanLength,ylength);

xfactor = HorzScanLength / xlength;
yfactor = VerticalScanLength / ylength;

number_of_radial_scans = 32;
number_of_B_scans = 9;
TotalImages = number_of_radial_scans * number_of_B_scans;

pathf='C:\Krishna Recent\Patents\AI for Premium IOLs DELL\IMAI Codes\OCT Image Processing\RetiEye\';
fnames = fn_findfiles('JPG',pathf);

% Register B Scans in each Radial Scan
k=1;
stFrame = 1;
endFrame = 1;% number_of_B_scans;

iReg1 = 0;
for i = stFrame: endFrame
   i1 = double(imread(fnames{i}));
   imageslices(1:ylength,1:xlength,i)=i1(:,:,1);
   iReg1 = iReg1 + (i1./255);
end
iReg0 = iReg1(:,:,1);

% Noise in image
ss=size(iReg0);
a =0;
r = a + (bnoise-a).*rand(ss(1),ss(2),1);
iReg0 = iReg0 + r;
%[mask2,mu,v,p]=EMSeg(iReg0,3);
%[mask3,mu,v,p]=EMSeg(iReg00,2);
%iReg0=mask3;

% Anterior surface
for k=1:xlength
    d=abs(diff(iReg0(:,k)));
    y=find(d == max(d));
    ylocA(k)=y(1);
    xlocA(k)=k;
    % Find a way to replace those max pixels with a lower value
end
ylocAOrig=ylocA;
xlocAOrig=xlocA;

% Remove outliers 
for t=1:2
Adiff=diff(ylocA);
yOut = find(abs(Adiff) < thresholdA);
m=1;
for k=1:length(yOut)
    xAOutliers(m)=xlocA(yOut(k));
    yAOutliers(m)=ylocA(yOut(k));
    m=m+1;
end

if (t==1)
    ylocA=yAOutliers;
    xlocA=xAOutliers;
    xAOutliers=0;
    yAOutliers=0;
    threshold = 5;   
end

end


set(figure(1), 'Position', [19 150 710 833]);                 
subplot(2,1,1);
imagesc(iReg0);hold on
plot(xlocA,ylocA,'r+');hold off

subplot(2,1,2);
imagesc(iReg0);hold on
plot(xAOutliers,yAOutliers,'r+');hold off
colormap('gray');

% choose set of (x,y) locations randomly and estimaate radius

for w=1:200
a=1;b=length(xAOutliers);
locValues = fix(a + (b-a).*rand(50,1));
xlocFit=xAOutliers(locValues);
ylocFit=yAOutliers(locValues);
% Best fit radius for Anterior surface
[xcA,ycA,Ra,~] = fn_circfit(xlocFit.*xfactor,ylocFit.*yfactor);
AntR(w)=Ra;
end
fprintf('Anterior radius %2.2f mm +/- %2.2f mm\n',mean(AntR),std(AntR));
D1 = (1.376-1)*1000./mean(AntR);
fprintf('Anterior Power %2.2f D\n',D1);



% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Travel down around "offset" pixels and find the gradient peak
% % Posterior surface
% % use maximum offset determined by average human corneal thickness
% % Make sure offset value is much lower than clinically observed
% % thin corneas - not diseased but average eye with thin corneas
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ylocA=ylocAOrig;
% xlocA=xlocAOrig;
% offset=50;
% for k=1:xlength
% d=abs(diff(iReg0(:,k)));
% ylocSecondPeak=d(ylocA(k)+offset:ylength-1);
% y0=find(ylocSecondPeak == max(ylocSecondPeak));
% if (isempty(y0) == 0)
%     ylocP(k)= y0(1)+offset+ylocA(k);
% end
% end
% 
% % Remove outliers 
% 
% for t=1:2
% Pdiff=diff(ylocP);
% yOut = find(abs(Pdiff) < thresholdP);
% m=1;
% for k=1:length(yOut)
%     xPOutliers(m)=xlocA(yOut(k));
%     yPOutliers(m)=ylocP(yOut(k));
%     m=m+1;
% end
% if (t==1)
%     ylocP=yPOutliers;
%     xlocA=xPOutliers;
%     xPOutliers=0;
%     yPOutliers=0;
%     threshold = 5;   
% end
% 
% end
% set(figure(2), 'Position', [731 145 710 833]);   
% 
% subplot(2,1,1);
% imagesc(iReg0);hold on
% plot(xlocA,ylocP,'r*');hold off
% 
% subplot(2,1,2);
% imagesc(iReg0);hold on
% plot(xPOutliers,yPOutliers,'r*');hold off
% 
% colormap('gray');
% % Best fit radius for Posterior surface
% [xcP,ycP,Rp,~] = fn_circfit(xPOutliers.*xfactor,yPOutliers.*yfactor);
% 
% D2 = ((1.336-1.376)*1000)/((Ra*(Rp/Ra)));
% fprintf('Posterior radius %2.2f mm\n',Rp);
% fprintf('Posterior Power %2.2f D\n',D2);
% 
% D12=D1-D2;
% fprintf('Total eye Power %2.2f D\n',D12);
