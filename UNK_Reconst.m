function [outVal, estimation_Val, trend_Val]   = UNK_Reconst(in_Trend, in_Residual, X, pos2fix, cov_Mat, mdl_Vario)

estimation_Val = []; trend_Val = []; outVal = [];
%In each Line in mdle_Vario
for(iT = 1 : length(cov_Mat))
    %The residual terms
    tempVal = [] ; parameters_Value = mdl_Vario(iT,:); tmp_Matrice = cov_Mat{iT};
    for(k = 1 : length(pos2fix))
        for(l = 1 : length(X))
            h_o = abs(pos2fix(k) - X(l));
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
        est_Tmp = 0; y = in_Residual(iT,:);
        for(l = 1 : length(tmp_Val))
            est_Tmp = est_Tmp + tmp_Val(l)*y(l) ;
        end
        estimation_Val_o(k) = est_Tmp ;
    end
    estimation_Val = [estimation_Val ; estimation_Val_o];
    
    %The trend
    trnd_Tmp = [];
    for(k = 1 : length(pos2fix))
        Nbh_1 = []; Nbh_2 = []; d1 = 1000000000000; d2 = -1000000000000;
        for(l = 1 : length(X))
            Cmp_Dst = pos2fix(k) - X(l);
            if((Cmp_Dst < 0) && (Cmp_Dst > d2))
                Nbh_2(1) = X(l) ; Nbh_2(2) = in_Trend(iT,l);
            end
            if((Cmp_Dst > 0) && (Cmp_Dst < d1))
                Nbh_1(1) = X(l) ; Nbh_1(2) = in_Trend(iT,l);
            end
        end
        trnd_Tmp(k) = ( (Nbh_2(2) - Nbh_1(2))/(Nbh_2(1) - Nbh_1(1)) )*pos2fix(k) + Nbh_2(2) - ( (Nbh_2(2) - Nbh_1(2))/(Nbh_2(1) - Nbh_1(1)) )*Nbh_2(1) ;
    end
    trend_Val = [trend_Val ; trnd_Tmp];
    
    outVal = [outVal ; (trnd_Tmp + estimation_Val_o)];
end

return