function circularprofile(data, name)
    imageSizeX = length(data.x);
    
%     %benchmark
%     imageSizeX = 201;
%     imageSizeY = 201;
%     [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
%     % Next create the circle in the image.
%     centerX = (imageSizeX +1) / 2; % Wherever you want.
%     centerY = (imageSizeY +1) / 2;
%     radius = 50;
%     circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
    % circlePixels is a 2D "logical" array.
    % Now, display it.
    %plot
    center = (imageSizeX+1)/2;
    radius = center;
    %     degree_line = 0;
    degree_line = [-pi:0.02:pi];
    int_circle = double.empty;
    for deg = degree_line
        point1.x = center;
        point1.y = center;
        point2.x = double(cos(deg)*radius+point1.x);
        point2.y = double(sin(deg)*radius+point1.y);
        lin_prof = improfile(data.zmask, [point1.x point2.x], [point1.y point2.y] , 'bicubic');
        if deg == -pi
            length_lin = length(lin_prof);
        end
        line([point1.x point2.x], [point1.y point2.y],'Color','red');
        int_circle = [int_circle, sum(lin_prof,'all','omitnan')*length_lin/length(lin_prof)];
    end
    hold off
    print([name '_overlay'], '-dpng');
    
    p5 = figure(5);
    plot(degree_line*180/pi, int_circle, 'Color', 'black');
    xlim([-180 180])
    hl = refline([0 0]);
    hl.Color = 'k';
    set(p5, 'Position', [10 10 500 400]);
    print([name '_circular profile'], '-dpng');
    writematrix(int_circle, [name '_circular profile']);
end