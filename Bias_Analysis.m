function [m_Val, std_Val, R_val] = Bias_Analysis(X, Y1, Y2, Index_Plot)

y1 = []; y2 = []; x = []; R_val = 0;
for(iT = 2 : length(Y1))
    y1(iT - 1) = Y1(iT); y2(iT - 1) = Y2(iT); x(iT - 1) = X(iT);
end

%frst_Bisec = linspace(0, 1.5*max(y1));

m_Val = mean(y1-y2); std_Val = sqrt(var(y1-y2)) ;
for(iT = 1 : length(y1))
    R_val = R_val + (y1(iT) - mean(y1))*(y2(iT) - mean(y2));
end
R_val = R_val / length(y2) ;
R_val = R_val / (sqrt(var(y1)) * sqrt(var(y2))) ;

Hist_Name = figure;
curent_Name = [pwd '\OutPuts\Validation Of Reconstruction\Histogramme_Rec_' Index_Plot '.png'];
subplot(1,2,1);
hist(y1, 50); grid on; xlabel Values; ylabel Frequenciy; title 'Histogram of the input';
subplot(1,2,2);
hist(y2, 50); grid on; xlabel Values; ylabel Frequenciy; title 'Histogram of the output';
print(Hist_Name, '-dpng', curent_Name);

plot_Name = figure;
curent_Name2 = [pwd '\OutPuts\Validation Of Reconstruction\Variogram_Rec_' Index_Plot '.png'];
subplot(1,2,1);
plot(x,y1,'r','linewidth',1.5); hold on; plot(x, y2, 'b', 'linewidth', 1.5); grid on; xlabel 'Lag-Distance'; ylabel 'g(h)'; title 'Experimental Variogram';
subplot(1,2,2);
plot(y1,y2,'*r','linewidth',1.5); hold on; plot(y1, y1, 'b', 'linewidth', 2.2); grid on; xlabel Input; ylabel Output; title 'Cross validating';
print(plot_Name, '-dpng', curent_Name2);


return