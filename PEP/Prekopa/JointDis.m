function val = JointDis(state,param,varargin)
val = prod(normcdf(state,param,param/2));
end