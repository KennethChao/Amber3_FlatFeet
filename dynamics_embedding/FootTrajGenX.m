function x = FootTrajGenX(t)
%FOOTTRAJGENX
%    X = FOOTTRAJGENX(T)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    02-Oct-2014 15:09:29

t2 = t.^2;
x = t.*t2.*(t.*-1.5e1+t2.*6.0+1.0e1).*(1.0./1.0e1);
