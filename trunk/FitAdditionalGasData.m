% Load the new gas data file
filename = 'AdditionalGasData.xlsx';
[~,sheets] = xlsfinfo(filename);
num = cell(length(sheets)-1,1);
%%
for ctFluid=1:length(sheets)-1
    [num{ctFluid},txt] = xlsread(filename,sheets{ctFluid});
    lstPhase(ctFluid) = txt(2,end);
end
swtVapor = strcmp(lstPhase,'vapor');

%% Fit polynomials & plot
R = 8.3144621;
betaV = zeros(length(sheets)-1,2);
betaT = zeros(length(sheets)-1,3);
tRange = zeros(length(sheets)-1,2);
% For gases He & H2 & N2 apply a different exponent in Sutherland's formula
% He:3.45, H2:3.35 N2:3.2 (for a better fit)
gammaV = ones(length(sheets)-1,1)*3;
gammaV(strcmp(sheets,'H2')) = 3.35;
gammaV(strcmp(sheets,'He')) = 3.45;
gammaV(strcmp(sheets,'N2')) = 3.2;
gammaV(~cellfun('isempty',strfind(sheets,'bLb'))) = 0;
for ctFluid=1:length(sheets)-1
    disp(ctFluid);
    %%
    if strfind(sheets{ctFluid},'bLb')
        T_V = (num{ctFluid}(1,1):5:num{ctFluid}(end,1))';
    else
        T_V = (200:20:2000)';
    end
    
    % Viscosity
    ctPar = 12;
    if swtVapor(ctFluid)
        % Gas - Sutherland's formula
        funV = @(b,T) b(1)*T.^(gammaV(ctFluid)/2)./(T+b(2)); 
        BetaVini = [1.57,118];
    else
        % Liquid - Arrhenius model
        funV = @(b,T) b(1)*exp(b(2)./(R*T)); 
        if strcmp(sheets{ctFluid},'H2ObLb')
            BetaVini = [1,17000]; % Water
        else
            BetaVini = [18,5000];
        end
    end
    betaV(ctFluid,:) = nlinfit(num{ctFluid}(:,1),num{ctFluid}(:,ctPar),funV,BetaVini);
    valV = funV(betaV(ctFluid,:),T_V);
    figure(1);
    plot(num{ctFluid}(:,1),num{ctFluid}(:,ctPar),'b+-',T_V,valV,'ro');
    title(sheets{ctFluid});
    ylabel(txt{1,ctPar});
    xlabel(txt{1,1});
    
    % Therm Cond
    ctPar = 13; 
    funT = @(b,T) b(1)+b(2)*sqrt(T)+b(3)*T;
    betaT(ctFluid,:) = nlinfit(num{ctFluid}(:,1),num{ctFluid}(:,ctPar),funT,[-0.007,0.0015,0.0002]);
    valT = funT(betaT(ctFluid,:),T_V);
    figure(2);
    plot(num{ctFluid}(:,1),num{ctFluid}(:,ctPar),'b+-',T_V,valT,'ro');
    title(sheets{ctFluid});
    ylabel(txt{1,ctPar});
    xlabel(txt{1,1});
    
    % Range
    tRange(ctFluid,:) = num{ctFluid}([1,end],1)';
    pause
end

%% Add to the MediumModel Data - Vapor only
load('@MediumModel\IdealGases.mat');
%%
% for ctFluid=1:length(sheets)-1
for ctFluid=find(swtVapor) % not liquids
    strMaster.(sheets{ctFluid}).c = [betaV(ctFluid,:),gammaV(ctFluid)];
    strMaster.(sheets{ctFluid}).d = betaT(ctFluid,:);
%     strMaster.(sheets{ctFluid}).tRange_cd = tRange(ctFluid,:);
end

save('@MediumModel\IdealGases.mat','strMaster');
