function outVal = ExperimentalVariogram_Fnc(profileDat, specString)

    outVal = {};
    AB_Val = []; MN_Val = []; Rho_Val = []; Sp_Val = [];
    
    Hst1_Name = figure;
    hist(profileDat.Rho); 
    curent_Name = ['C:\Users\Kasys\Documents\MATLAB\MATLAB\02_04_2021_OK_SYNTH_1\OutPuts\Original Histogram\Histogram - Profiles ' num2str(profileDat.ID) '.png'];
    print(Hst1_Name, '-dpng', curent_Name);
    
    %Retrieve values for AB and MN
%     for(k = 1 : length(profileDat.xA))
%         AB_Val = [AB_Val (profileDat.xA{k} - profileDat.xB{k})];
%         MN_Val = [MN_Val (profileDat.xM{k} + profileDat.xN{k})];
%     end
    MN_Val = profileDat.x ; AB_Val = profileDat.y ;
    figure;
    plot(MN_Val, AB_Val, 'o');
    switch(specString)
        case 'Rho'
%             for(k = 1 : length(profileDat.Rho))
%                 Rho_Val = [Rho_Val (profileDat.Rho{k})];
%             end
            
            Rho_Val = profileDat.Rho ;
            %Compute the experimental varioagram of the profile
            listH = [];
            tol_Dist = 5;

            for k = 1 : length(AB_Val)
                for l = 1: length(AB_Val)
                    tmpDist = sqrt( (MN_Val(l) - MN_Val(k))^2 + (AB_Val(l) - AB_Val(k))^2 );
        
                    if(Check_In_Vec(listH, tmpDist, tol_Dist) == 0)
                        listH = [listH tmpDist];
                    end
                end
            end
            listH = sort(listH, 'ascend');
            disp(listH);

            g_H = [] ;
            for p = 1 : length(listH)
                N_H = [];
                for k = 1 : length(MN_Val)
                    for l = 1: length(MN_Val)
                        tmpDist = sqrt( (MN_Val(l) - MN_Val(k))^2 + (AB_Val(l) - AB_Val(k))^2 );
                        b_Sup = listH(p) + tol_Dist ;
                        b_Inf = listH(p) - tol_Dist ;
        
        
                        if( tmpDist >= b_Inf && tmpDist <= b_Sup)
                            N_H = [N_H; [Rho_Val(k) Rho_Val(l)]];
                        end
                    end
                end

                S = 0;
                for k = 1 : size(N_H, 1)
                    S = S + (N_H(k,1) - N_H(k,2))^2 ;
                end

                g_H = [g_H S*(1/(2*size(N_H,1)))];
            end
            outVal = {listH, g_H};
            H1_Name = figure;
            plot(listH, g_H, '*-');
            curent_Name = ['C:\Users\Kasys\Documents\MATLAB\MATLAB\02_04_2021_OK_SYNTH_1\OutPuts\Experimental Model\Exp_Model - Profiles ' profileDat.ID '.png'];
            print(H1_Name, '-dpng', curent_Name);

        case 'Sp'
            for(k = 1 : length(profileDat.Rho))
                Sp_Val = [Sp_Val (profileDat.Sp{k})];
            end
            %Compute the experimental varioagram of the profile
            listH = [];
            tol_Dist = 2;

            for k = 1 : length(AB_Val)
                for l = 1: length(AB_Val)
                    tmpDist = sqrt( (MN_Val(l) - MN_Val(k))^2 + (AB_Val(l) - AB_Val(k))^2 );
        
                    if(Check_In_Vec(listH, tmpDist, tol_Dist) == 0)
                        listH = [listH tmpDist];
                    end
                end
            end
            listH = sort(listH, 'ascend');
            disp(listH);

            g_H = [] ;
            for p = 1 : length(listH)
                N_H = [];
                for k = 1 : length(MN_Val)
                    for l = 1: length(MN_Val)
                        tmpDist = sqrt( (MN_Val(l) - MN_Val(k))^2 + (AB_Val(l) - AB_Val(k))^2 );
                        b_Sup = listH(p) + tol_Dist ;
                        b_Inf = listH(p) - tol_Dist ;
        
        
                        if( tmpDist >= b_Inf && tmpDist <= b_Sup)
                            N_H = [N_H; [Sp_Val(k) Sp_Val(l)]];
                        end
                    end
                end

                S = 0;
                for k = 1 : size(N_H, 1)
                    S = S + (N_H(k,1) - N_H(k,2))^2 ;
                end

                g_H = [g_H S*(1/(2*size(N_H,1)))];
            end
            outVal = {listH, g_H};
            figure;
            plot(listH, g_H, '*-');
    end
    
end
