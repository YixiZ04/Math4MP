nsim = 2;

[poptot,pops,nec,act,o2con] = dcode_simfiles(nsim);

%% PET

nsnapshots = size(poptot,2);
snapshot = round(nsnapshots*0.95);
figpos = get(groot, 'ScreenSize');
ventana = figure('Name','In-silico PET slice','Units','pixels','OuterPosition',[100 100 2*figpos(4)/3 2*figpos(4)/3]);
colormap(flipud(gray))
imagesc(act{snapshot}(:,:,40))
axis tight equal
axis off
set(gcf,'color','w');

%% CIRC

nsnapshots = size(poptot,2);
snapshot = round(nsnapshots*0.95);
figpos = get(groot, 'ScreenSize');
ventana = figure('Name','Circ map','Units','pixels','OuterPosition',[100 100 2*figpos(4)/3 2*figpos(4)/3]);
heatmap(o2con{snapshot}(:,:,40))
set(gcf,'color','w');

%% MRI

nsnapshots = size(poptot,2);
snapshot = round(nsnapshots*0.95);
figpos = get(groot, 'ScreenSize');
w = figure('Name','In-silico MRI slice','Units','pixels','OuterPosition',[100 100 2*figpos(4)/3 2*figpos(4)/3]);
colormap(gray)
imagesc(poptot{snapshot}(:,:,40))
axis tight equal
axis off
set(gcf,'color','w');

%% MRI BY POPULATIONS

nsnapshots = size(poptot,2); 
snapshot = round(nsnapshots*0.95);
% Loop through each population
for popIdx = 1:5
    figpos = get(groot, 'ScreenSize');
    w = figure('Name',['In-silico MRI slice of population ', num2str(popIdx)],'Units','pixels','OuterPosition',[100 100 2*figpos(4)/3 2*figpos(4)/3]);
    colormap(gray)
    imagesc(pops{snapshot}(:,:,40, popIdx))
    axis tight equal
    axis off
    set(gcf,'color','w');
end

%% 3D

col1 = [27 161 226]./255; % Cyan 
nsnapshots = size(poptot,2); 
snapshot = round(nsnapshots*0.55); 
popT = poptot{snapshot}; 
figpos = get(groot, 'ScreenSize');
w = figure('Name','3D Tumor rendering','Units','pixels','OuterPosition',[0 0 figpos(4) figpos(4)]);
hold on
data1=smooth3(popT,'gaussian',5);
patch(isocaps(data1,13000),'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(data1,13000),'FaceColor',col1,'EdgeColor','none');
isonormals(data1,p1);
view(3)
hold off
fig = gcf;
fig.Renderer = 'zbuffer';
camlight left
lighting gouraud
material shiny
legend('3d tumor')
set(gca, 'Color', 'none');
axis equal tight off

%% 3D BY POPULTAIONS

col1 = [27 161 226]./255;   % Cyan
col2 = [226 75 75]./255;    % Coral Red
col3 = [100 200 100]./255;  % Mint Green
col4 = [255 200 100]./255;  % Amber
col5 = [180 120 220]./255;  % Lavender

thresholds = [20, 10000, 10000, 10000, 42000];
alphas = [0.2, 0.2, 0.2, 0.2, 0.2]; % opacity

nsnapshots = size(poptot,2);
snapshot = round(nsnapshots*0.95);
figpos = get(groot, 'ScreenSize');
w = figure('Name','3D Tumor rendering','Units','pixels','OuterPosition',[0 0 figpos(4) figpos(4)]);
hold on

% Assuming poptot{snapshot} is a cell array of 5 volumes
for popIdx = 1:5
    popVolume = pops{snapshot}(:,:,:,popIdx); % Adjust indexing as needed
    data = smooth3(popVolume, 'gaussian', 5);
    threshold = thresholds(popIdx); % Adjust based on your cell densities
    
    p = patch(isosurface(data, threshold), ...
              'FaceColor', eval(['col' num2str(popIdx)]), ...
              'EdgeColor', 'none', ...
              'FaceAlpha', alphas(popIdx)); % Transparency helps see overlaps
    isonormals(data, p);
end

legend('pop1', 'pop2', 'pop3', 'pop4', 'pop5')
view(3)
fig = gcf;
fig.Renderer = 'zbuffer';
camlight left
lighting gouraud
material shiny
set(gca, 'Color', 'none');
axis equal tight off