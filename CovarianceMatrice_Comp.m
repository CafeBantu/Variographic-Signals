function [outPut, estimation_Val] = CovarianceMatrice_Comp(mdl_Vario, X, Y)
outPut ={}; estimation_Val = {};
x = linspace(min(X), max(X));
estimation_Val_o = zeros(1, length(x));
for(iT = 1 : size(mdl_Vario, 1))
    parameters_Value = mdl_Vario(iT,:);
    y = Y(iT,:);
    for(k = 1 : length(X))
        for(l = 1 : length(X))
            h = abs(X(k) - X(l));
            switch(parameters_Value(1))
                case 1
                    %spherical corelation
                    if(h <= parameters_Value(2))
                        tmp_Matrice(k,l) = parameters_Value(3) + (parameters_Value(4)-parameters_Value(3))*(1.5*h/parameters_Value(2) - 0.5*(h/parameters_Value(2))^3);
                    else
                        tmp_Matrice(k,l) = parameters_Value(4);
                    end
                case 2
                    %exponential corelation
                    tmp_Matrice(k,l) = parameters_Value(3) + (parameters_Value(4) - parameters_Value(3))*(1 - exp(-3*(h/parameters_Value(2))));
                case 3
                    %gaussian corelation
                    tmp_Matrice(k,l) = parameters_Value(3) + (parameters_Value(4) - parameters_Value(3))*(1 - exp(-3*(h/parameters_Value(2))^2));
            end
        end
    end
    
    for(k = 1 : length(x))
        for(l = 1 : length(X))
            h_o = abs(x(k) - X(l));
            switch(parameters_Value(1))
                case 1
                    %spherical corelation
                    if(h_o <= parameters_Value(2))
                        tmp_Matrice_o(l) = parameters_Value(3) + (parameters_Value(4)-parameters_Value(3))*(1.5*h_o/parameters_Value(2) - 0.5*(h_o/parameters_Value(2))^3);
                    else
                        tmp_Matrice_o(l) = parameters_Value(4);
                    end
                case 2
                    %exponential corelation
                    tmp_Matrice_o(l) = parameters_Value(3) + (parameters_Value(4) - parameters_Value(3))*(1 - exp(-3*(h_o/parameters_Value(2))));
                case 3
                    %gaussian corelation
                    tmp_Matrice_o(l) = parameters_Value(3) + (parameters_Value(4) - parameters_Value(3))*(1 - exp(-3*(h_o/parameters_Value(2))^2));
            end
        end
        tmp_Val = tmp_Matrice_o/tmp_Matrice ;
        est_Tmp = 0;
        for(l = 1 : length(tmp_Val))
            est_Tmp = est_Tmp + tmp_Val(l)*y(l) ;
        end
        estimation_Val_o(k) = est_Tmp ;
    end
    estimation_Val = [estimation_Val {estimation_Val_o}];
    outPut = [outPut {tmp_Matrice}];
end

return