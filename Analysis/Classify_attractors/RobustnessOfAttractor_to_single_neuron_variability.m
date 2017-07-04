%%% does variability in single neurons correspond to variability between
%%% programs?
% (1) global change in participation
% (2) local change in participation
clear all; close all

addpath ../../Functions
load ParticipationChanges_SamePCs

blnCommon = 0;

%% do comparison for commn axes projections
if blnCommon
    load RecurrenceManifoldStats_CommonAxes_AllPreps
else
    load RecurrenceManifoldStats_AllPreps
end

load RecurrenceStats RcrStats  % to get density data


Robust.allH = []; Robust.allDpart = []; Robust.allHaussD = []; Robust.allHraw = []; Robust.allPchange = []; Robust.allMDensity = []; Robust.allHD_Diff = [];
Robust.allHD_MagCtrl = [];

pbins = 0:0.05:1; % bins for histogram

ixpairs = [1,2; 2,3; 1,3];  % pairs stored sequentially in e.g. Coeffs.dMaxP field

for iPrep = 1:numel(Prep)
    if ~isempty(Prep(iPrep).Dmin)
        Robust.Density(iPrep,:) = [RcrStats(iPrep,:).densAllPeriod];

        % then we have computed the distances between manifolds etc
        
        pbinsRaw = 0:0.05:max(max(Coeffs(iPrep).sumWProg));
        Pchange = zeros(3) + inf;
        for i = 1:3
            for j=i+1:3
                % so compute the distances between the distributions of participation
                p1 = hist(Coeffs(iPrep).normWMax(i,:),pbins); p1 = p1 ./sum(p1);
                p2 = hist(Coeffs(iPrep).normWMax(j,:),pbins); p2 = p2 ./sum(p2);               
                H(i,j) = sqrt(sum((sqrt(p1) - sqrt(p2)).^2)) ./ sqrt(2); 
  
                % distances between raw participation
                p1 = hist(Coeffs(iPrep).sumWProg(i,:),pbinsRaw); p1 = p1 ./sum(p1);
                p2 = hist(Coeffs(iPrep).sumWProg(j,:),pbinsRaw); p2 = p2 ./sum(p2);               
                Hraw(i,j) = sqrt(sum((sqrt(p1) - sqrt(p2)).^2)) ./ sqrt(2); 

                % and the total change
                Dpart(i,j) = sum(abs(Coeffs(iPrep).normWMax(i,:) - Coeffs(iPrep).normWMax(j,:)));
                
                MDensity(i,j) = mean([Robust.Density(iPrep,i),Robust.Density(iPrep,j)]); % mean density of recurrence points in pair
            end
            % proportion of neurons that change greater than threshold on
            % this pair of programs
            Pchange(ixpairs(i,1),ixpairs(i,2)) = sum(Coeffs(iPrep).dMaxP(:,i) >= Coeffs(iPrep).MinVarPart) ./ numel(Coeffs(iPrep).Rank); 
            % keyboard
        end
        
        % collate into vectors for scatter plot
        Robust.allH = [Robust.allH; nonzeros(H)];
        Robust.allHraw = [Robust.allHraw; nonzeros(Hraw)];
        Robust.allDpart = [Robust.allDpart; nonzeros(Dpart)];
        Robust.allPchange = [Robust.allPchange; Pchange(~isinf(Pchange))];
        Robust.allHaussD = [Robust.allHaussD; nonzeros(Prep(iPrep).HaussD)];
        % keyboard
        ctrl = mean(Prep(iPrep).HaussDCtrl,3);
        Robust.allHD_Diff = [Robust.allHD_Diff; nonzeros(Prep(iPrep).HaussD) - nonzeros(ctrl)];  % difference between (mean) control and data distance
        Robust.allHD_MagCtrl = [Robust.allHD_MagCtrl; nonzeros(Prep(iPrep).HaussD) ./nonzeros(ctrl)];  % expressed as magnitude

        Robust.allMDensity = [Robust.allMDensity; nonzeros(MDensity)];
    end
end

%% fit models
tbl = table(Robust.allH,Robust.allHaussD,'VariableNames',{'NormDist','HaussD'});
Robust.lm_NormDist =fitlm(tbl,'linear');  % normalised participation distribution

tbl = table(Robust.allHraw,Robust.allHaussD,'VariableNames',{'RawDist','HaussD'});
Robust.lm_RawDist =fitlm(tbl,'linear');

tbl = table(Robust.allDpart,Robust.allHaussD,'VariableNames',{'ChangePart','HaussD'});
Robust.lm_ChangePart =fitlm(tbl,'linear');
Robust.lm_ChangePartRob =fitlm(tbl,'linear','RobustOpts','bisquare');

% lm_ChangePartQ =fitlm(tbl,'purequadratic'); % ,'RobustOpts','bisquare');

tbl = table(Robust.allPchange,Robust.allHaussD,'VariableNames',{'Pchange','HaussD'});
Robust.lm_Pchange = fitlm(tbl,'linear');
Robust.lm_Pchange_Rob =fitlm(tbl,'linear','RobustOpts','bisquare');

% fit models to normalised distance between programs
% omit difference > 0 as then these are different programs!
tbl = table(Robust.allH(Robust.allHD_Diff < 0),Robust.allHD_Diff(Robust.allHD_Diff < 0),'VariableNames',{'NormDist','HaussD'});
Robust.lm_NormDist_N =fitlm(tbl,'linear');

tbl = table(Robust.allDpart(Robust.allHD_Diff < 0),Robust.allHD_Diff(Robust.allHD_Diff < 0),'VariableNames',{'ChangePart','HaussD'});
Robust.lm_ChangePart_N =fitlm(tbl,'linear');
Robust.lm_ChangePartRob_N =fitlm(tbl,'linear','RobustOpts','bisquare');

tbl = table(Robust.allPchange(Robust.allHD_Diff < 0),Robust.allHD_Diff(Robust.allHD_Diff < 0),'VariableNames',{'Pchange','HaussD'});
Robust.lm_Pchange_N = fitlm(tbl,'linear');
Robust.lm_Pchange_Rob_N =fitlm(tbl,'linear','RobustOpts','bisquare');

% normalised as magnitude of control
tbl = table(Robust.allH(Robust.allHD_MagCtrl < 1),Robust.allHD_MagCtrl(Robust.allHD_MagCtrl < 1),'VariableNames',{'NormDist','HaussD'});
Robust.lm_NormDist_M =fitlm(tbl,'linear');

tbl = table(Robust.allDpart(Robust.allHD_MagCtrl < 1),Robust.allHD_MagCtrl(Robust.allHD_MagCtrl < 1),'VariableNames',{'ChangePart','HaussD'});
Robust.lm_ChangePart_M =fitlm(tbl,'linear');
Robust.lm_ChangePartRob_M =fitlm(tbl,'linear','RobustOpts','bisquare');

tbl = table(Robust.allPchange(Robust.allHD_MagCtrl < 1),Robust.allHD_MagCtrl(Robust.allHD_MagCtrl < 1),'VariableNames',{'Pchange','HaussD'});
Robust.lm_Pchange_M = fitlm(tbl,'linear');
Robust.lm_Pchange_Rob_M =fitlm(tbl,'linear','RobustOpts','bisquare');

if blnCommon
    save RobustnessNeuron_vs_Program_Common Robust
else
    save RobustnessNeuron_vs_Program_Separate Robust
end
%% show  - do any neuron particpation changes predict manifold distance?

% raw distance
x = 0.1:0.01:0.45;
y = Robust.lm_NormDist.Coefficients.Estimate(1) + Robust.lm_NormDist.Coefficients.Estimate(2).*x;

figure
% subplot(211), 
plot(Robust.allH,Robust.allHaussD,'ko'); hold on
% scatter(allH,allHaussD,20*allMDensity,'k'); hold on
plot(x,y,'r')
axis([0.1 0.45 6 26])
ylabel('Distance between programs')
xlabel('Participation (normalised) distance')
text(0.125,24,['R2 = ' num2str(Robust.lm_NormDist.Rsquared.Ordinary)],'Color',[1 0 0])
text(0.125,23,['P = ' num2str(Robust.lm_NormDist.Coefficients.pValue(2))],'Color',[1 0 0]);

% subplot(212), 
% plot(allHraw,allHaussD,'ko')
% ylabel('Distance between programs')
% xlabel('Participation (raw) distance')
% exportPPTfig(gcf,'Distance_vs_distance',[10 15 5 5])

% plot with robust regression line
x = 4:0.1:25;
y = Robust.lm_ChangePartRob.Coefficients.Estimate(1) + Robust.lm_ChangePartRob.Coefficients.Estimate(2).*x;
figure
plot(Robust.allDpart,Robust.allHaussD,'ko'); hold on
plot(x,y,'r')
axis([4 24 6 26])
ylabel('Distance between programs')
xlabel('Total participation change')
text(16,8,['R2 = ' num2str(Robust.lm_ChangePartRob.Rsquared.Ordinary)],'Color',[1 0 0])
text(16,7,['P = ' num2str(Robust.lm_ChangePartRob.Coefficients.pValue(2))],'Color',[1 0 0]);

%exportPPTfig(gcf,'Total_Variation_vs_distance',[10 15 5 5])


x = 0:0.01:0.3;
y = Robust.lm_Pchange.Coefficients.Estimate(1) + Robust.lm_Pchange.Coefficients.Estimate(2).*x;
figure
plot(Robust.allPchange,Robust.allHaussD,'ko'); hold on
plot(x,y,'r')
ylabel('Distance between programs')
xlabel('Proportion of variable neurons')
text(0.12,19,['R2 = ' num2str(Robust.lm_Pchange.Rsquared.Ordinary)],'Color',[1 0 0])
text(0.12,18,['P = ' num2str(Robust.lm_Pchange.Coefficients.pValue(2))],'Color',[1 0 0]);

%exportPPTfig(gcf,'NeuronsVariation_vs_distance',[10 15 5 5])

%% normalised to control
figure
% subplot(211), 
line([0.1 0.45],[0 0],'Color','k'); hold on
plot(Robust.allH,Robust.allHD_Diff,'ko'); hold on
% scatter(allH,allHaussD,20*allMDensity,'k'); hold on
%plot(x,y,'r')
%axis([0.1 0.45 6 26])
ylabel('Norm. Distance between programs')
xlabel('Participation (normalised) distance')
%text(0.125,24,['R2 = ' num2str(lm_NormDist.Rsquared.Ordinary)],'Color',[1 0 0])
%text(0.125,23,['P = ' num2str(lm_NormDist.Coefficients.pValue(2))],'Color',[1 0 0]);

% subplot(212), 
% plot(allHraw,allHaussD,'ko')
% ylabel('Distance between programs')
% xlabel('Participation (raw) distance')
% exportPPTfig(gcf,'Distance_vs_distance',[10 15 5 5])

% plot with robust regression line
%x = 4:0.1:25;
%y = lm_ChangePartRob.Coefficients.Estimate(1) + lm_ChangePartRob.Coefficients.Estimate(2).*x;
figure
line([4 24],[0 0],'Color','k'); hold on
plot(Robust.allDpart,Robust.allHD_Diff,'ko'); hold on
%plot(x,y,'r')
%axis([4 24 6 26])
ylabel('Norm. Distance between programs')
xlabel('Total participation change')
%text(16,8,['R2 = ' num2str(lm_ChangePartRob.Rsquared.Ordinary)],'Color',[1 0 0])
%text(16,7,['P = ' num2str(lm_ChangePartRob.Coefficients.pValue(2))],'Color',[1 0 0]);

%exportPPTfig(gcf,'Total_Variation_vs_distance',[10 15 5 5])


%x = 0:0.01:0.3;
%y = lm_Pchange.Coefficients.Estimate(1) + lm_Pchange.Coefficients.Estimate(2).*x;
figure
line([0 0.3],[0 0],'Color','k'); hold on
plot(Robust.allPchange,Robust.allHD_Diff,'ko'); hold on
%plot(x,y,'r')
ylabel('Norm. Distance between programs')
xlabel('Proportion of variable neurons')
%text(0.12,19,['R2 = ' num2str(lm_Pchange.Rsquared.Ordinary)],'Color',[1 0 0])
%text(0.12,18,['P = ' num2str(lm_Pchange.Coefficients.pValue(2))],'Color',[1 0 0]);

%% normalised to control by magnitude

figure
line([0.1 0.45],[1 1],'Color','k'); hold on
plot(Robust.allH,Robust.allHD_MagCtrl,'ko'); hold on
ylabel('Norm. Distance between programs')
xlabel('Participation (normalised) distance')

figure
line([4 24],[1 1],'Color','k'); hold on
plot(Robust.allDpart,Robust.allHD_MagCtrl,'ko'); hold on
ylabel('Norm. Distance between programs')
xlabel('Total participation change')

figure
line([0 0.3],[1 1],'Color','k'); hold on
plot(Robust.allPchange,Robust.allHD_MagCtrl,'ko'); hold on
ylabel('Norm. Distance between programs')
xlabel('Proportion of variable neurons')

%% does distance relate to number of embedding dimensions?

%% do rates also predict distances?

% figure
% line([0.1 0.45],[0 0],'Color','k'); hold on
% plot(allH,allHD_Diff,'ko'); hold on
% %axis([0.1 0.45 6 26])
% ylabel('Norm. Distance between programs')
% xlabel('Change in rate')
% %text(0.125,24,['R2 = ' num2str(lm_NormDist.Rsquared.Ordinary)],'Color',[1 0 0])
% %text(0.125,23,['P = ' num2str(lm_NormDist.Coefficients.pValue(2))],'Color',[1 0 0]);


