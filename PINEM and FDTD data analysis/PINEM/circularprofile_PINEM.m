% Written by Jiahui Li
% This is for calculating the circularprofile for PINEM data with center
% correction.

function circularprofile_PINEM(data,name,center)
    figure(16)
    tiledlayout(1,1);
    ax = nexttile;
    imshow(data)
    load('BWR.mat')
    colormap(ax, BWR)
    caxis([-200 200]);
    title('Left minus Right');
    hold on
    imageSizeX = length(data);

    %plot
    radius = min([imageSizeX-center(1), imageSizeX-center(2), center(1), center(2)],[],'all');
    degree_line = [-pi:0.02:pi];
    int_circle = double.empty;
    for deg = degree_line
        point1.x = center(1);
        point1.y = center(2);
        point2.x = double(cos(deg)*radius+point1.x);
        point2.y = double(sin(deg)*radius+point1.y);
        lin_prof = improfile(data, [point1.x point2.x], [point1.y point2.y] , 'bicubic');
        if deg == -pi
            length_lin = length(lin_prof);
        end
        line([point1.x point2.x], [point1.y point2.y],'Color','red');
        int_circle = [int_circle, sum(lin_prof,'all','omitnan')*length_lin/length(lin_prof)];
    end
    hold off
    print([name '_overlay'], '-dpng');
    
    p8 = figure(17);
    plot(degree_line*180/pi, int_circle, 'Color', 'black');
    xlim([-180 180])
    hl = refline([0 0]);
    hl.Color = 'k';
    set(p8, 'Position', [10 10 500 400]);
    writematrix(int_circle, [name '_circular profile']);
    print([name '_circular profile'], '-dpng');
end