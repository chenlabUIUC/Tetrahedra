% Written by Jiahui Li
% This code is for computing the electric field distribution from FDTD
% simulation. 

clear all;
clf;
modelname ='largeoffset-04-3'; %filename
sizenotmatch = 0; %0: default, 1: adjustment when size doesn't match
wavelength = [570 620]; %electric field distribution at what wavelength
position_only = [-1:1:1]; 
for wl = wavelength
    for position_n = position_only
        switch position_n
            case -1
                position = '-bottom';
            case 0
                position = '-middle';
            case 1
                position = '-top';
        end

        filename_L = [modelname '/' modelname '-LCP-' num2str(wl) '-E' position '.mat'];
        filename_R = [modelname '/' modelname '-RCP-' num2str(wl) '-E' position '.mat'];
        L = load(filename_L);
        R = load(filename_R);
        
        if length(position_only) == 1
            position = [position '-' position '-only'];
        end
        
        if sizenotmatch == 1
            if position_n == 1
                LMR.x = L.lum.x(1,1:199);
                LMR.y = R.lum.y;
                LMR.z = L.lum.z(1:199,1:201)-R.lum.z(1:199,1:201);
            else
                LMR.x = L.lum.x(1,3:201);
                LMR.y = R.lum.y;
                LMR.z = L.lum.z(3:201,1:201)-R.lum.z(3:201,1:201);
            end
        else
            LMR.x = L.lum.x;
            LMR.y = R.lum.y;
            LMR.z = L.lum.z-R.lum.z;
        end
        
        if length(position_only) == 1
            LMR.zall = LMR.z;
        else
            if position_n == -1
                LMR.zall = LMR.z;
                Lall = L;
                Rall = R;
            else
                LMR.zall = LMR.zall + LMR.z;
                Lall.lum.z = Lall.lum.z + L.lum.z;
                Rall.lum.z = Rall.lum.z + R.lum.z;
            end
        end
        
        sum_LMRz = sum(LMR.z, 'all');
        disp([num2str(sum_LMRz), ' ', position])
        
        figure(1)

        p1 = tiledlayout(1,3);
        ax1 = nexttile;
        surfaceplot_gray(L);
        title(['Left ' num2str(wl)]);
        caxis([0 10])
        
        ax2 = nexttile;
        surfaceplot_gray(R);
        title(['Right ' num2str(wl)]);
        caxis([0 10])
        
        ax3 = nexttile;
        surfaceplot_LMR(LMR);
        load('BWR.mat')
        colormap(ax3, BWR)
        caxis([-5 5])
        title(['Left minus Right ' num2str(wl)]);
        ImgName = [modelname '/' 'Rescaled-' modelname '-' num2str(wl) position '-L-minus-R' '.png'];
        exportgraphics(p1, ImgName, 'Resolution', 300);
    end
    
    sum_LMRzall = sum(LMR.zall, 'all');
    disp([num2str(sum_LMRzall), ' all positions'])

    maxLMRzall = max_LMR(LMR.zall);
    
    p2 = figure(2);
    tiledlayout(1,1);
    ax = nexttile;
    surfaceplot_LMR_all(LMR);
    load('BWR.mat') %Colormap
    colormap(ax, BWR)
    caxis([-5 5])
    title(['Left minus Right as total three ' num2str(wl)]);
    ImgName2 = [modelname '/' 'Rescaled-' modelname '-' num2str(wl) position '-L-minus-R-all' '.png'];
    exportgraphics(p2, ImgName2, 'Resolution', 300);
    
    p3 = figure(3);
    tiledlayout(1,2);
    ax = nexttile;
    surfaceplot_gray(Lall);
    title(['Left ' num2str(wl)]);
    caxis([0 10])
        
    ax2 = nexttile;
    surfaceplot_gray(Rall);
    title(['Right ' num2str(wl)]);
    caxis([0 10])
    ImgName3 = [modelname '/' 'Rescaled-' modelname '-' num2str(wl) position '-LandR' '.png'];
    exportgraphics(p3, ImgName3, 'Resolution', 300);
    %masking
    maskname = ['mask/' modelname '.tif'];
    mask = imread(maskname);
    mask = imresize(mask, [length(LMR.y) length(LMR.y)], 'lanczos3');
    maskbi = imbinarize(mask);
    LMR.zmask = LMR.zall;
    for i = 1: length(LMR.x)
        for j = 1: length(LMR.y)
            if maskbi(i, j) == 0
                LMR.zmask(length(LMR.x)-j+1,i) = 0;
            end
        end
    end

    %plot masked figure
    p4 = figure(4);
    tiledlayout(1,1);
    ax = nexttile;
    surfaceplot_LMR_all_mask(LMR);
    load('BWR.mat')
    colormap(ax, BWR)
    maxLMRmasked = max_LMR(LMR.zmask);
    caxis([-40 40])
    title('Left minus Right as total three, masked');

    ImgName3 = [modelname '/' modelname num2str(wavelength) position '-L-minus-R-allmasked' '.png'];
    exportgraphics(p4, ImgName3, 'Resolution', 300)

    % plot circular profile
    figure(4)
    p5 = tiledlayout(1,1);
    ax = nexttile;
    imshow(LMR.zmask)
    load('BWR.mat')
    colormap(ax, BWR)
    maxLMRmasked = max_LMR(LMR.zmask);
    caxis([-40 40])
    title('Left minus Right as total three, masked');
    hold on
    circularprofile(LMR, modelname);
end

function surfp = surfaceplot_gray(data)
    [X,Y]=meshgrid(min(data.lum.y):(max(data.lum.y)-min(data.lum.y))/(length(data.lum.y)-1):max(data.lum.y),min(data.lum.x):(max(data.lum.x)-min(data.lum.x))/(length(data.lum.x)-1):max(data.lum.x));
    surf(X,Y,data.lum.z,'EdgeColor','none');
    axis ij
    colormap('gray')
    colorbar
    axis equal
    view(90,90)
end

function surfp = surfaceplot_LMR(data)
    [X,Y]=meshgrid(min(data.y):(max(data.y)-min(data.y))/(length(data.y)-1):max(data.y),min(data.x):(max(data.x)-min(data.x))/(length(data.x)-1):max(data.x));
    surf(X,Y,data.z,'EdgeColor','none');
    axis ij
    colorbar
    axis equal
    view(90,90)
end

function surfp = surfaceplot_LMR_all(data)
    [X,Y]=meshgrid(min(data.y):(max(data.y)-min(data.y))/(length(data.y)-1):max(data.y),min(data.x):(max(data.x)-min(data.x))/(length(data.x)-1):max(data.x));
    surf(X,Y,data.zall,'EdgeColor','none');
    axis ij
    colorbar
    axis equal
    view(90,90)
end

function surfp = surfaceplot_LMR_all_mask(data)
    [X,Y]=meshgrid(min(data.y):(max(data.y)-min(data.y))/(length(data.y)-1):max(data.y),min(data.x):(max(data.x)-min(data.x))/(length(data.x)-1):max(data.x));
    surf(X,Y,data.zmask,'EdgeColor','none');
    axis ij
    colorbar
    axis equal
    view(90,90)
end

function maxLMR = max_LMR(data)
    if max(data,[],'all')+min(data,[],'all') >=0
        maxLMR = max(data,[],'all');
    else
        maxLMR = -min(data,[],'all');
    end
end 