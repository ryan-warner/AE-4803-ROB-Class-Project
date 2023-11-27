function saveImage(figure)
%figure = findobj( 'Type', 'Figure', 'Name', figure);
rez=600; %resolution (dpi) of final graphic. 600 is a good number.
figpos=getpixelposition(figure);
resolution=get(0,'ScreenPixelsPerInch');
set(figure,'paperunits','inches','papersize',figpos(3:4)/resolution,...
 'paperposition',[0 0 figpos(3:4)/resolution]);
path='/Users/ryanwarner/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/Senior Year/Fall 2023/Robotics and Autonomy/HW 3/HW3/plots';
if sum(strfind(figure.Name,"/")) > 0
    indices = strfind(figure.Name,"/");
    tempName = figure.Name;
    tempName(indices) = "o";
    name=sprintf("%s.png",tempName); %What should it be called?
else
    name=sprintf("%s.png",figure.Name); %What should it be called?
end
print(figure,fullfile(path,name),'-dpng',['-r',num2str(rez)],'-opengl')
end
