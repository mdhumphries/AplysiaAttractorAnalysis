% Angela Bruno & Mark Humphries 31/10/16

clear all; close all

if ispc
    filepath = 'C:\Users\mqbssmhg.DS\Dropbox\Analyses\My modularity based methods\Collaborations using my methods\AB - Aplysia motor program dynamics\Mapping\';
else
    filepath = '/Users/mqbssmhg/Dropbox/Analyses/My modularity based methods/Collaborations using my methods/AB - Aplysia motor program dynamics/Mapping/';
end

run figure_properties

load([filepath 'Participation']);  % NOTE: these maps use extracted Participation for ease of use - if updating Participation analysis, must update here too
load([filepath 'MapCoords'])
load([filepath 'diode_template_adjusted']);
load([filepath 'FileTable']);

LEFT = [1,2,4,10];
RIGHT = [3,5,6,9];

Mscale = 10;  % for scaling MarkerSize to in units of "Participation"
Rscale = 20;

xtiles = 5;  % how many divisions of Participation range?

% axis scaling to make diode square
ratio = 125.5/131.9;

Group = LEFT;

% check mapping
FileTable{LEFT}
FileTable{RIGHT}

nPreps = numel(Participation);

%% Range and Max participation. %arranged into quintiles for the sake of making maps
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

%%Participation Max quintiles 
for  i = 1:nPreps
% ParticipationMax{i} = mean(Participation(i).Normalised);
ParticipationMax{i} = max(Participation(i).Normalised);
end
for i = 1:nPreps
    MaxM(i) = max(max(ParticipationMax{1,i}));
end
for i = 1:nPreps
    MinM(i) = min(min(ParticipationMax{1,i}));
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
    PM = ParticipationMax{1,F}; 
    for iX = 1:xtiles
        GroupMax{iX,F} = find(PM >= tilesM(:,iX) & PM < tilesM(:,iX+1));
    end
end 


%% Make Summary Map 
rsize = [3,10,20,40]; %largest marker size to represent the greatest range in participation
msize = [3,8,10,20];   %mean size - make smaller 
cmap = brewermap(xtiles,'OrRd');
diodemap = [1 1 1; 1 1 0.9];

%Make Max map
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 5]); hold on
imagesc(diode_template); colormap(diodemap); set(gca,'YDir','normal'); hold on
for T = Group; %split files by left or right RIGHT =[3,5,6,9] LEFT =[1,2,4,10]
    for iT = 1:xtiles
        pcolor = cmap(iT,:); %one color per Xtile
        if numel(GroupMax{iT,T} > 0);
    %      plot(MapX{1,T}(GroupMax{iT,T})*10,MapY{1,T}(GroupMax{iT,T})*10,...
    %          'o','MarkerEdgeColor','k',...
    %                'MarkerFaceColor',pcolor,...
    %                'MarkerSize','MarkerSize',msize(iT)); 
         scatter(MapX{1,T}(GroupMax{iT,T})*10,MapY{1,T}(GroupMax{iT,T})*10,...
             (ParticipationMax{T}(GroupMax{iT,T})*Mscale).^2,...
             'MarkerEdgeColor','k',...
             'MarkerFaceColor',pcolor);   %'MarkerSize',msize(iT)); 
                    hold on %4 sizes for quintiles
        end 
    end
end
axis tight

% axis square  % for when the diode template is truly square...
pbaspect([1 ratio 1])
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1:3) pos(4) * ratio]);

axis off

if Group == LEFT
%     text(400,490,'Rostral','Fontsize',7);
%     text(490,450,'Medial','Rotation',270,'Fontsize',7);
    print(['Participation_Map_Left'],'-depsc')

else
%     text(1,490,'Rostral','Fontsize',7);
%     text(10,375,'Medial','Rotation',90,'Fontsize',7);
    print(['Participation_Map_Right'],'-depsc')
end


%% make range map
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 5 5]); hold on
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
end
axis tight
% axis square
pbaspect([1 ratio 1])

axis off
if Group == LEFT
%     text(400,490,'Rostral','Fontsize',7);
%     text(490,450,'Medial','Rotation',270,'Fontsize',7);
    print(['Range_Map_Left'],'-depsc')

else
%     text(1,490,'Rostral','Fontsize',7);
%     text(10,375,'Medial','Rotation',90,'Fontsize',7);
    print(['Range_Map_Right'],'-depsc')
    
end

