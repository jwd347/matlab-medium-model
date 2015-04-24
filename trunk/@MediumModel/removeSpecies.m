function removeSpecies(medium,cellGas)
persistent strMaster
if isempty(strMaster) % Checks if database has been loaded already
    load IdealGases  % Load database of relevant values
end
if ~iscell(cellGas)
    cellGas={cellGas};
end
[cellGas,~,indMatch]=intersect(cellGas,medium.names); % Only existing names
if isempty(cellGas)
    error('MediumModel:GasNotFound','gas name was not found in current model.')
end
medium.names(indMatch)=[];
medium.notCondensed(indMatch)=[];
for ctGas =1:length(cellGas)
    try
        medium.gas = rmfield(medium.gas,cellGas{ctGas});
        medium.index = rmfield(medium.index,cellGas{ctGas});
    catch ME
        error('MediumModel:GasNotSuppoted','Gas name not recognised : %s',cellGas{ctGas})
    end
end
for ctGasLeft =1:length(medium.names)
    try
        medium.index.(medium.names{ctGasLeft})=ctGasLeft;
    catch ME
        error('MediumModel:GasNotSuppoted','Gas name not recognised : %s',cellGas{ctGas})
    end
end

newZ = medium.Z;
newZ(indMatch)=[];
medium.Z=newZ./sum(newZ);

medium.CheckShape;
if size(medium.Z,1)==1
    medium.Zeq=repmat(medium.Z,length(medium.T),1);
else
    medium.Zeq=medium.Z;
end
medium.props;    %% Get the properties of the gas/gas mixture

end
