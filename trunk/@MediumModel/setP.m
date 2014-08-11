function setP(me,P,varargin)


if isempty(varargin) 
    me.P=P;
else 
P0 = varargin{1};   
me.P0=P0;
me.P=P;
end

me.CheckShape;
me.props;

end
