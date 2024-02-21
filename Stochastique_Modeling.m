function out_Mdl = Stochastique_Modeling(Stride_Val, inputVal, prf_Pos)
out_Mdl = []; P1 = [1 2 3];

%In tihs function we are going to outcome the linear model linking the
%experimental variogram for eac profile

%Frist of all we step on a stride
 
tol_Dist = 1 ;

for(iT = 1 : size(inputVal, 1))
    
    vario_Val = [];
    for(k = 1 : length(Stride_Val))
        N_H = [];
        for(l = 1 : length(prf_Pos))
            for(p = 1 : length(prf_Pos))
                tmpDist = abs(prf_Pos(l) - prf_Pos(p));
                b_Sup = Stride_Val(k) + tol_Dist ;
                b_Inf = Stride_Val(k) - tol_Dist ;
            
                if(tmpDist >= b_Inf && tmpDist <= b_Sup)
                    N_H = [N_H; inputVal(iT,l) inputVal(iT,p)];
                end
            end
        end
    
        S = 0;
        for(l = 1 : size(N_H, 1))
            S = S + (N_H(l,1) - N_H(l,2))^2 ;
        end
        vario_Val = [vario_Val S*(1/(2*size(N_H,1)))];
    end
    
%     figure;
%     plot(Stride_Val, vario_Val, '*-'); hold on;
%     
%   Automatically fit the experimental model of variogram
    m_Vario = mean(vario_Val); v_Vario = var(vario_Val); s_Vario = sqrt(v_Vario);
    P2 = linspace(round(min(Stride_Val)), round(max(Stride_Val)), 2000); P3 = linspace(0, round(abs(m_Vario-s_Vario)), 2000); P4 = linspace(round(abs(m_Vario)), round(max(vario_Val)), 2000);
    Ftd_Par = GeneticAlgo_Modeling(P1, P2, P3, P4, Stride_Val, vario_Val);
%     plot_By_Param(Ftd_Par(1), Ftd_Par(2), Ftd_Par(3), Ftd_Par(4), Stride_Val);
%     
    out_Mdl = [out_Mdl ; Ftd_Par];

end

return