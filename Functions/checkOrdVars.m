function varargout = checkOrdVars(n,varargin)

%Check the number of input data and expand it in a vector of n component if necessary
for k = 1:nargout
    if isnumeric(varargin{k}) && size(varargin{k},1)==1
        varargout{k}=varargin{k}.*ones(n,1);
    elseif ischar(varargin{k}) && size(varargin{k},1)==1
        varargout{k}=cell(n,1);
        for i=1:length(varargout{k})
            varargout{k}{i}=varargin{k};
        end
%     elseif isnumeric(varargin{k}) && size(varargin{k},3)>1
%         varargout{k}=varargin{k}'.*ones(n,1);
    else
        varargout{k}=varargin{k};
    end
end
end