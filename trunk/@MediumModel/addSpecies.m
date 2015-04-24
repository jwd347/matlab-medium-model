function addSpecies(medium,cellGas)
persistent strMaster
if isempty(strMaster) %~exist('strMaster','var') Checks if database has been loaded already
    load IdealGases  %% Load database of relevant values
end
if ~iscell(cellGas)
    cellGas={cellGas};
end
cellGas=setdiff(cellGas,medium.names); % Only new names
if ~all(ismember(cellGas,fieldnames(strMaster)))
    error('MediumModel:GasNotFound','A gas name was not found.')
end
if isempty(cellGas), return; end;
ctGasOrig=length(medium.names);
medium.names=[medium.names,cellGas];
medium.notCondensed=[medium.notCondensed,nan(1,length(cellGas))];
for ctGas =1:length(cellGas)
    try
        medium.gas.(cellGas{ctGas})=strMaster.(cellGas{ctGas});
        medium.index.(cellGas{ctGas})=ctGas+ctGasOrig;
        medium.notCondensed(ctGas+ctGasOrig)=~strMaster.(cellGas{ctGas}).swtCondensed;
    catch ME
        error('MediumModel:GasNotSuppoted','Gas name not recognised : %s',cellGas{ctGas})
    end
end

medium.Z=[medium.Z,zeros(size(medium.Z,1),ctGas)];

medium.CheckShape;
if size(medium.Z,1)==1
    medium.Zeq=repmat(medium.Z,length(medium.T),1);
else
    medium.Zeq=medium.Z;
end
medium.props;    %% Get the properties of the gas/gas mixture
end
