% Written by Jiahui Li
% This is for calculating the subtracted electric field in PINEM experiment

clear all;
clf;
% File name
for i = 17 %9 10 16 17 19 24 29
    try
    ra=111; %rotation correction to compare with FDTD simulation
    leftx =109; %reorient the center of domain
    rightx = 230;
    topy=78;
    bottomy=239;
    %Read TEM, LCP and RCP image
    ImgName1 = [num2str(i)];
    ImgNameleft = [ImgName1,'_L.tif'];
    ImgNameright = [ImgName1,'_R.tif'];

    TEM = double(imread([ImgName1 '.tif']));
    left = double(imread(ImgNameleft));
    right = double(imread(ImgNameright));
    sub = left-right;

    figure(1)
    p1 = tiledlayout(2,3);
    ax1 = nexttile;
    imagesc(left);
    colormap('gray');
    caxis([0 255]);
    colorbar('eastoutside')
    pbaspect([1 1 1]);
    set(gca,'TickDir','out');

    ax2 = nexttile;
    imagesc(right);
    colormap('gray');
    caxis([0 255]);
    colorbar('eastoutside')
    pbaspect([1 1 1]);
    set(gca,'TickDir','out');

    ax3 = nexttile;
    imagesc(sub);
    load('BWR.mat')
    colormap(ax3,BWR)
    colorbar('eastoutside')
    caxis([-200 200]);
    pbaspect([1 1 1]);
    set(gca,'TickDir','out');

    ax4 = nexttile;
    histogram(left(:),'BinWidth',50)
    xlim([0 255])
    set(gca,'TickDir','out');
    ax5 = nexttile;
    histogram(right(:),'BinWidth',50)
    xlim([0 255])
    set(gca,'TickDir','out');
    ax6 = nexttile;
    histogram(sub(:),'BinWidth',0.2)
    xlim([-200 200])
    set(gca,'TickDir','out');

    ImgName2 = [num2str(i),'-L-minus-R-redraw'];
    print(ImgName2,'-dpng');

    figure(2)
    p2 = tiledlayout(1,3);
    ax1 = nexttile;
    imagesc(TEM);
    colormap('gray');
    pbaspect([1 1 1]);
    set(gca,'TickDir','out');

    ax2 = nexttile;
    ar_TEM = imrotate(TEM,ra);
    imagesc(ar_TEM);
    colormap('gray');
    pbaspect([1 1 1]);
    set(gca,'TickDir','out');

    ra2 = 90;
    ax3 = nexttile;
    ar_TEM2 = imrotate(ar_TEM,ra2);
    imagesc(ar_TEM2);
    colormap('gray');
    pbaspect([1 1 1]);
    set(gca,'TickDir','out');

    ImgName3 = [num2str(i),'TEMrotate'];
    print(ImgName3,'-dpng');

    figure(3)
    ar_left = imrotate(left,ra);
    ar_right = imrotate(right,ra);
    ar_sub = imrotate(sub,ra);
    p3 = tiledlayout(1,3);
    ax1 = nexttile;
    imagesc(ar_left);
    colormap('gray');
    caxis([0 255]);
    colorbar('eastoutside')
    pbaspect([1 1 1]);
    set(gca,'TickDir','out');

    ax2 = nexttile;
    imagesc(ar_right);
    colormap('gray');
    caxis([0 255]);
    colorbar('eastoutside')
    pbaspect([1 1 1]);
    set(gca,'TickDir','out');

    ax3 = nexttile;
    imagesc(ar_sub);
    load('BWR.mat')
    colormap(ax3,BWR)
    colorbar('eastoutside')
    caxis([-200 200]);
    pbaspect([1 1 1]);
    set(gca,'TickDir','out');

    ImgName4 = [num2str(i),'-sub-rotate'];
    print(ImgName4,'-dpng');
    catch
        continue
    end
end

%% Plot circular profile
center = [(leftx+rightx)/2;(topy+bottomy)/2];

ax1 = figure(4);
ar_sub = imrotate(ar_sub, ra2);
imagesc(ar_sub);
load('BWR.mat')
colormap(ax1,BWR)
colorbar('eastoutside')
caxis([-200 200]);
pbaspect([1 1 1]);
set(gca,'TickDir','out');

circularprofile_PINEM(ar_sub, ImgName1, center);

%% 
