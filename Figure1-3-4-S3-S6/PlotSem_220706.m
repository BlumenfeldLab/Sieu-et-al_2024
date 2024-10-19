% Edited by Jerry. 07/06/2022
function PlotSem_220706(DataX,DataY,se,ax,color)
    %   PLOT_SEM Plot SE/SEM using rectangular drawing function
    %   Currently plots SE; in order to change to SEM just divide SE by the
    %   number of points.
    %   DataX: x values of data on which SEM plot is centered
    %   DataY: y values of data on which SEM plot is centered
    %   se: standard error vector corresponding to each data point.
    %   color: e.g. 'b' - blue, 'g' - green
    if ~isrow(DataX)
        DataX = DataX';
    end
    if ~isrow(DataY)
        DataY = DataY';
    end
    if ~isrow(se)
        se = se';
    end
    X = [DataX,flip(DataX)];
    Y = [DataY+se,flip(DataY-se)];
    fill(ax,X,Y,color,'FaceAlpha',0.2);
end