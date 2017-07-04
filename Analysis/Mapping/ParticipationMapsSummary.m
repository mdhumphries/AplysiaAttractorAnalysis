% Angela Bruno & Mark Humphries 31/10/16

clear all; close all

load Participation
load MapCoords
load diode_template_adjusted
load FileTable

LEFT = [1,2,4,10];
RIGHT = [3,5,6,9];

Mscale = 25;  % for scaling MarkerSize to in units of "Participation"
Rscale = 50;

xtiles = 5;  % how many divisions of Participation range?

Group = LEFT;

% check mapping
FileTable{LEFT}
FileTable{RIGHT}

nPreps = numel(Participation);

%% Range and Mean participation. %arranged into quintiles for the sake of making maps
 % Range
for  i = 1:nPreps
    ParticipationRange{i} = max(Participation(i).Normalised) - min(Participation(i).Normalised);
end

%Range
for i = 1:nPreps
    MaxR(i) = max(max(ParticipationRange{1,i}));
end
for i = 1:nPreps
    MinR(i) = min(min(ParticipationRange{1,i}));
end

%%Participation Mean quintiles 
for  i = 1:nPreps
% ParticipationMean{i} = mean(Participation(i).Normalised);
ParticipationMean{i} = max(Participation(i).Normalised);
end
for i = 1:nPreps
    MaxM(i) = max(max(ParticipationMean{1,i}));
end
for i = 1:nPreps
    MinM(i) = min(min(ParticipationMean{1,i}));
end

%% divide into Xtiles, across all pooled preps

tilesR = linspace(min(MinR),max(MaxR),xtiles+1); %quintiles of participation range

for F = 1:nPreps
    PR = ParticipationRange{1,F};  
    for iX = 1:xtiles
        GroupRange{iX,F} = find(PR > tilesR(:,iX) & PR < tilesR(:,iX+1));
    end
end


tilesM = linspace(min(MinM),max(MaxM),xtiles+1); %xtiles of participation mean

for F = 1:nPreps
    PM = ParticipationMean{1,F}; 
    for iX = 1:xtiles
        GroupMean{iX,F} = find(PM >= tilesM(:,iX) & PM < tilesM(:,iX+1));
    end
end 


%% Make Summary Map 
rsize = [3,10,20,40]; %largest marker size to represent the greatest range in participation
msize = [3,8,10,20];   %mean size - make smaller 
cmap = brewermap(xtiles,'OrRd');
diodemap = [1 1 1; 1 1 0.9];

%Make Mean map
figure
imagesc(diode_template); colormap(diodemap); set(gca,'YDir','normal'); hold on
for T = Group; %split files by left or right RIGHT =[3,5,6,9] LEFT =[1,2,4,10]
for iT = 1:xtiles
    pcolor = cmap(iT,:); %one color per Xtile
    if numel(GroupMean{iT,T} > 0);
%      plot(MapX{1,T}(GroupMean{iT,T})*10,MapY{1,T}(GroupMean{iT,T})*10,...
%          'o','MarkerEdgeColor','k',...
%                'MarkerFaceColor',pcolor,...
%                'MarkerSize','MarkerSize',msize(iT)); 
     scatter(MapX{1,T}(GroupMean{iT,T})*10,MapY{1,T}(GroupMean{iT,T})*10,...
         (ParticipationMean{T}(GroupMean{iT,T})*Mscale).^2,...
         'MarkerEdgeColor','k',...
         'MarkerFaceColor',pcolor);   %'MarkerSize',msize(iT)); 

                hold on %4 sizes for quintiles
    end 
end
% title('Normalised Participation Mean Summary')
if Group == LEFT
    axis tight
    axis square
    text(450,490,'Rostral','Fontsize',7);
    text(490,475,'Medial','Rotation',90,'Fontsize',7);
else
    
end
axis off

%% make range map
figure
imagesc(diode_template); colormap(diodemap); set(gca,'YDir','normal'); hold on

for T = Group; %split files by left or right Left =[3,5,6,9] Right=[1,2,4,10]
for iT = 1:4
    pcolor = cmap(iT,:); %one color per Xtile
    if numel(GroupRange{iT,T} > 0);
     scatter(MapX{1,T}(GroupRange{iT,T})*10,MapY{1,T}(GroupRange{iT,T})*10,...
         (ParticipationRange{T}(GroupRange{iT,T})*Rscale).^2,...
         'MarkerEdgeColor','k',...
         'MarkerFaceColor',pcolor); 
        hold on
    end 
end
% title('Normalised Participation Range Summary')
end
axis off