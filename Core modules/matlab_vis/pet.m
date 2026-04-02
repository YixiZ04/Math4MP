nsim = 4;

[poptot,pops,nec,act] = dcode_simfiles(nsim);

%% PET

nsnapshots = size(poptot,2);
snapshot = round(nsnapshots*0.55);
figpos = get(groot, 'ScreenSize');
ventana = figure('Name','In-silico PET slice','Units','pixels','OuterPosition',[100 100 2*figpos(4)/3 2*figpos(4)/3]);
colormap(flipud(gray))
imagesc(act{snapshot}(:,:,40))
axis tight equal
axis off
set(gcf,'color','w');

%% MRI

nsnapshots = size(poptot,2);
snapshot = round(nsnapshots*0.85);
figpos = get(groot, 'ScreenSize');
w = figure('Name','In-silico MRI slice','Units','pixels','OuterPosition',[100 100 2*figpos(4)/3 2*figpos(4)/3]);
colormap(gray)
imagesc(poptot{snapshot}(:,:,40))
axis tight equal
axis off
set(gcf,'color','w');

%% 3D

col1 = [27 161 226]./255; % Cyan 
nsnapshots = size(poptot,2); 
snapshot = round(nsnapshots*0.95); 
popT = poptot{snapshot}; 
figpos = get(groot, 'ScreenSize');
w = figure('Name','3D Tumor rendering','Units','pixels','OuterPosition',[0 0 figpos(4) figpos(4)]);
hold on
data1=smooth3(popT,'gaussian',5);
patch(isocaps(data1,40000),'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(data1,40000),'FaceColor',col1,'EdgeColor','none');
isonormals(data1,p1);
view(3)
hold off
fig = gcf;
fig.Renderer = 'zbuffer';
camlight left
lighting gouraud
material shiny
set(gca, 'Color', 'none');
axis equal tight off