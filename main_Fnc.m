 clear, clc;
 Strating_Time = cputime;

%Loading Data from iles
% loodingData;

%Here let us start by creating our synthetic dataset containing the
%tomographies
x_List = 5 : 15 : 750;
y_List = 5 : 15 : 200;

profile_Tomo = {};
for(iT = 1 : 100)
    curent_Tomo.ID = ['Profile ' num2str(iT)];
    curent_Tomo.x = [];
    curent_Tomo.y = [];
    curent_Tomo.Rho = [];
    curent_Tomo.Rho_ExpVariogram = {};
    
    for(k = 1 : length(x_List))
        for(l = 1 : length(y_List))
                curent_Tomo.x = [curent_Tomo.x x_List(k)] ;
                curent_Tomo.y = [curent_Tomo.y y_List(l)] ;
                z_List(k,l) = randi([500 5000], 1, 1);
                curent_Tomo.Rho = [curent_Tomo.Rho z_List(k,l)] ; %Ressitivities are irregularly distributed in the space domain
        end
    end
    
    profile_Tomo = [profile_Tomo {curent_Tomo}];
    
    pFl_Name = figure;
    imagesc(x_List, y_List, z_List); set(gca, 'YDir', 'normal'); 
    curent_Name = ['C:\Users\Administrator\Documents\MATLAB\UNK_SYNTH_1\OutPuts\Profiles\Profiles ' num2str(iT) '.png'];
    print(pFl_Name, '-dpng', curent_Name);    
end
    

%Compute each profile's experimental variogram
for(k = 1 : length(profile_Tomo))
    outPut1 = ExperimentalVariogram_Fnc(profile_Tomo{k}, 'Rho');
%     outPut2 = ExperimentalVariogram_Fnc(profile_Tomo{k}, 'Sp');
    
    profile_Tomo{k}.Rho_ExpVariogram = outPut1;
%     profile_Tomo{k}.Sp_ExpVariogram = outPut2;
end

Stride_ToUse = profile_Tomo{1}.Rho_ExpVariogram{1};
Prf_Pos = 1:100;

%In tihs function we are going to outcome the linear model linking the
%experimental variogram for eac profile

%Frist of all we step on a stride
listH = [] ;
tol_Dist = 2 ;

for(k = 1 : length(Stride_ToUse))
    for(l = 1 : length(profile_Tomo))
        for(t = 1 : length(profile_Tomo{l}.Rho_ExpVariogram{1}))
            if(profile_Tomo{l}.Rho_ExpVariogram{1}(t) == Stride_ToUse(k))
                gH_Val(k,l) = profile_Tomo{l}.Rho_ExpVariogram{2}(t);
                break;
            end
        end
    end
    for(k = 1 : length(Prf_Pos))
        for(l = 1 : length(Prf_Pos))
            tmpDist = abs(Prf_Pos(k) - Prf_Pos(l));
            if(Check_In_Vec(listH, tmpDist, tol_Dist) == 0)
                listH = [listH tmpDist];
            end
        end
    end
end
    listH = sort(listH, 'ascend');
    disp(listH);
%The first step is to process the dichotomy and guess the trend with the
%FFT Algorithm
gH_Val_Trend = []; gH_Val_Residual = [];
Param_K = [60 5 5];
for(k = 1 : size(gH_Val,1))
    [outMdl, outBias, outVal] =  DFT_Modeling(k, Prf_Pos, gH_Val(k,:), Param_K, 'Trend Removal');

    gH_Val_Trend(k,:) = outVal ;   
    gH_Val_Residual(k,:) = abs(gH_Val_Trend(k,:) - gH_Val(k,:)) ;
end

%Let us for this study retreive the spectral value at p = 65, i.e at k = 7
Prf_Pos_Dmg = [];
for(k = 1 : length(Prf_Pos))
    if(k ~= 10 && k ~= 20 && k ~= 30 && k ~= 40 && k ~= 50 && k ~= 60)
        Prf_Pos_Dmg = [Prf_Pos_Dmg Prf_Pos(k)];
    end
end

gH_Val_Dmg = []; gH_Val_Trend_Dmg = []; gH_Val_Residual_Dmg = [];
for(k = 1 : size(gH_Val,1))
    tempVec = []; tempVec2 = []; tempVec3 = [];
    for(l = 1 : size(gH_Val, 2))
        
        if(l ~= 10 && l ~= 20 && l ~= 30 && l ~= 40 && l ~= 50 && l ~= 60)
            tempVec = [tempVec gH_Val(k,l)];
            tempVec2 = [tempVec2 gH_Val_Trend(k,l)];
            tempVec3 = [tempVec3 gH_Val_Residual(k,l)];
        end
    end
    gH_Val_Dmg = [gH_Val_Dmg ; tempVec];
    gH_Val_Trend_Dmg = [gH_Val_Trend_Dmg ; tempVec2];
    gH_Val_Residual_Dmg = [gH_Val_Residual_Dmg ; tempVec3];
end

% Model_FromProfiles1 = Stochastique_Modeling(listH, gH_Val_Residual, Prf_Pos);
% 
% %Let us compute the covariance matrices for each spectrum
% [Covariance_Data, Estimated_Data] = CovarianceMatrice_Comp(Model_FromProfiles1, Prf_Pos, gH_Val);
% for(k = 1 : size(gH_Val,1))
%     Mdl_fig = figure;
%     plot(Prf_Pos, gH_Val(k,:), 'r*');
%     hold on;
% 
%     x = linspace(min(Prf_Pos), max(Prf_Pos));
%     plot(x, Estimated_Data{1,k},'b-');
%     curent_Name = ['C:\Users\Kasys\Documents\MATLAB\MATLAB\13_07_2021_UNK_SYNTH_1\OutPuts\Reconstructed Singals\Var_Signals_G - ' num2str(k) '.png'];
%     print(Mdl_fig, '-dpng', curent_Name); 
%     
%     Mdl_fig2 = figure;
%     hist(gH_Val(k,:));
%     curent_Name = ['C:\Users\Kasys\Documents\MATLAB\MATLAB\13_07_2021_UNK_SYNTH_1\OutPuts\Original Histogram\Histo_Signals_G - ' num2str(k) '.png'];
%     print(Mdl_fig2, '-dpng', curent_Name); 
% end
% 
% 

Model_FromProfiles2 = Stochastique_Modeling(listH, gH_Val_Residual_Dmg, Prf_Pos_Dmg);

%Let us compute the covariance matrices for each spectrum
[Covariance_Data2, Estimated_Data2] = CovarianceMatrice_Comp(Model_FromProfiles2, Prf_Pos_Dmg, gH_Val_Residual_Dmg);
for(k = 1 : size(gH_Val_Residual_Dmg,1))
    Mdl_fig21 = figure;
    plot(Prf_Pos_Dmg, gH_Val_Residual_Dmg(k,:), 'r*');
    hold on;
    
    curent_Name = ['C:\Users\Administrator\Documents\MATLAB\UNK_SYNTH_1\OutPuts\Reconstructed Singals\Var_Signals_R - ' num2str(k) '.png'];
    print(Mdl_fig21, '-dpng', curent_Name); 

    x2 = linspace(min(Prf_Pos_Dmg), max(Prf_Pos_Dmg));
    plot(x2, Estimated_Data2{1,k},'b-');
    
    Mdl_fig22 = figure;
    hist(gH_Val_Residual_Dmg(k,:));
    curent_Name = ['C:\Users\Administrator\Documents\MATLAB\UNK_SYNTH_1\OutPuts\Original Histogram\Histo_Signals_R - ' num2str(k) '.png'];
    print(Mdl_fig22, '-dpng', curent_Name);
end

%Now the reconstruction on the Damaged position
Pos_2_Fix = [10 20 30 40 50 60]; gH_Val_Fixed = [];
[out1, out2, out3] = UNK_Reconst(gH_Val_Trend_Dmg, gH_Val_Residual_Dmg, Prf_Pos_Dmg, Pos_2_Fix, Covariance_Data2, Model_FromProfiles2);

gH_Val_Org = [gH_Val(:, 10) gH_Val(:, 20) gH_Val(:, 30) gH_Val(:, 40) gH_Val(:, 50) gH_Val(:, 60)]; gH_Val_Fixed = out1;
%Bias anaysis
mean_Values = []; std_Values = []; R_Values = [];
for(k = 1 : length(Pos_2_Fix))
    
    Y1 = gH_Val_Org(:, k) ; Y2 = gH_Val_Fixed(:, k);
    
     [m_val, std_Val, r_val] = Bias_Analysis(Stride_ToUse, Y1, Y2, num2str(k));
    mean_Values = [mean_Values m_val]; std_Values = [std_Values std_Val]; R_Values = [R_Values r_val];
end


Ending_Time = cputime - Strating_Time;
disp(Ending_Time); disp('TO EXECUTE THE CODE');