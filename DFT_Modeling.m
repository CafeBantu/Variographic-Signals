function [OutPut_Mdl, Bias_Val, Y_Fit] = DFT_Modeling(IndexPlot, X, Y, Param_List, nature_Sig)

N = length(X);
K = Param_List(1);

dx = max(X)/(N) ;

H2_Name = figure;

%plot(X, Y, '*-r'); %hold on;

%Here we have to compute the value of A and B, the coefficients of the
%Discrete Fourier Serie fitting the variographic signals

for(k = 1 : K)
    A(k) = 0; B(k) = 0;
    for(l = 1 : N)
        A(k) = A(k) + Y(l)*cos(k*X(l)*pi/N)*dx/max(X) ;
        B(k) = B(k) + Y(l)*sin(k*X(l)*pi/N)*dx/max(X) ;
    end
    
    %X_Fit = linspace(min(X), max(X), 100);
    for(l = 1 : length(X))
        Y_Fit(l) = 0.5*sum(Y.*ones(size(X)))/N;
        for(m = 1 : length(A))
            Y_Fit(l) = Y_Fit(l) + A(m)*cos(m*X(l)*pi/N) + B(m)*sin(m*X(l)*pi/N);
        end
    end
    subplot(1, 2, 1);
    plot(X, Y_Fit); grid on; hold on;
    pause(.2);
    
    Bias_Dat = Y - Y_Fit ;
    subplot(1, 2, 2);
    plot(X, Bias_Dat); grid on;
    pause(.2);
end

Bias_Val = Y - Y_Fit ;
%hold on;
subplot(1, 2, 1);
plot(X, Y, '*r'); pause(.2);

White_Dat = zeros(size(X));
subplot(1, 2, 2);
hold on;
plot(X, White_Dat, 'r-*');

% curent_Name = ['C:\Users\Kasys\Documents\MATLAB\13_07_2021_UNK_SYNTH_1\OutPuts\Reconstructed Singals\Stride_Index_' num2str(IndexPlot) nature_Sig '-Recons_Signal.png'];
% print(H2_Name,'-dpng',curent_Name);

OutPut_Mdl = [A; B];
return