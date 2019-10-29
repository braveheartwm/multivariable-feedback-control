% Model for the Room Heat application
function [ret,x0,str,ts,xts]=rhmodel(t,x,u,flag);
% RHMODEL is the M-file description of the SIMULINK system named RHMODEL.
%	The block-diagram can be displayed by typing: RHMODEL.
%
%	SYS=RHMODEL(T,X,U,FLAG) returns depending on FLAG certain
%	system values given time point, T, current state vector, X,
%	and input vector, U.
%	FLAG is used to indicate the type of output to be returned in SYS.
%
%	Setting FLAG=1 causes RHMODEL to return state derivatives, FLAG=2
%	discrete states, FLAG=3 system outputs and FLAG=4 next sample
%	time. For more information and other options see SFUNC.
%
%	Calling RHMODEL with a FLAG of zero:
%	[SIZES]=RHMODEL([],[],[],0),  returns a vector, SIZES, which
%	contains the sizes of the state vector and other parameters.
%		SIZES(1) number of states
%		SIZES(2) number of discrete states
%		SIZES(3) number of outputs
%		SIZES(4) number of inputs
%		SIZES(5) number of roots (currently unsupported)
%		SIZES(6) direct feedthrough flag
%		SIZES(7) number of sample times
%
%	For the definition of other parameters in SIZES, see SFUNC.
%	See also, TRIM, LINMOD, LINSIM, EULER, RK23, RK45, ADAMS, GEAR.

% Note: This M-file is only used for saving graphical information;
%       after the model is loaded into memory an internal model
%       representation is used.
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: rhmodel.m,v 1.4 2004/02/05 12:23:14 vidaral Exp $


% the system will take on the name of this mfile:
sys = mfilename;
new_system(sys)
%simver(1.3)
simver(3.0)
if (0 == (nargin + nargout))
     set_param(sys,'Location',[173,215,727,517])
     open_system(sys)
end;
set_param(sys,'Solver',		'ode45')
set_param(sys,'StopTime',	'1000')
set_param(sys,'StartTime', 	'0.0')
set_param(sys,'MaxStep',	'2')
set_param(sys,'MinStep',	'1')
set_param(sys,'RelTol','1e-3')
%set_param(sys,'Return vars',   '')

add_block('built-in/Transport Delay',[sys,'/','Measurement delay'])
set_param([sys,'/','Measurement delay'],...
		'orientation','left',...
		'Delay Time','thetam',...
		'position',[235,235,280,265])

add_block('built-in/Outport',[sys,'/','y'])
set_param([sys,'/','y'],...
		'position',[525,145,545,165])

add_block('built-in/Step Fcn',[sys,'/','d'])
set_param([sys,'/','d'],...
		'Time','0',...
		'SampleTime','0',...
		'After','d',...
		'position',[275,40,295,60])

add_block('built-in/Outport',[sys,'/','u'])
set_param([sys,'/','u'],...
		'orientation','left',...
		'Port','2',...
		'position',[145,85,165,105])

add_block('built-in/Sum',[sys,'/','Sum'])
set_param([sys,'/','Sum'],...
		'position',[460,145,480,165])

add_block('built-in/Transfer Fcn',[sys,'/','Gd'])
set_param([sys,'/','Gd'],...
		'Numerator','[10]',...
		'Denominator','[1000 1]',...
		'position',[340,32,410,68])

add_block('built-in/Transfer Fcn',[sys,'/','G'])
set_param([sys,'/','G'],...
		'Numerator','[20]',...
		'Denominator','[1000 1]',...
		'position',[350,140,435,180])

add_block('built-in/Transfer Fcn',[sys,'/','PID Controller'])
set_param([sys,'/','PID Controller'],...
		'Numerator','conv([taui 1],[taud 1])*kc',...
		'Denominator','taui*[0.1*taud 1 0]',...
		'position',[140,138,320,182])

add_block('built-in/Step Fcn',[sys,'/','r'])
set_param([sys,'/','r'],...
		'orientation','left',...
		'Time','0',...
		'SampleTime','0',...
		'After','r',...
		'position',[25,60,45,80])

add_block('built-in/Sum',[sys,'/','Sum1'])
set_param([sys,'/','Sum1'],...
		'inputs','+-',...
		'position',[105,150,125,170])

add_block('built-in/Transfer Fcn',[sys,'/','Kr'])
set_param([sys,'/','Kr'],...
		'Numerator','[3]',...
		'Denominator','[150 1]',...
		'position',[25,133,80,177])
add_line(sys,[440,160;455,160])
add_line(sys,[415,50;445,50;455,150])
add_line(sys,[325,160;345,160])
add_line(sys,[85,155;100,155])
add_line(sys,[300,50;335,50])
add_line(sys,[130,160;135,160])
add_line(sys,[485,155;520,155])
add_line(sys,[492,155;492,250;285,250])
add_line(sys,[230,250;87,250;87,165;100,165])
add_line(sys,[332,160;332,95;170,95])
add_line(sys,[20,70;7,70;7,155;20,155])

drawnow

% Return any arguments.
if (nargin | nargout)
	% Must use feval here to access system in memory
	if (nargin > 3)
		if (flag == 0)
			eval(['[ret,x0,str,ts,xts]=',sys,'(t,x,u,flag);'])
		else
			eval(['ret =', sys,'(t,x,u,flag);'])
		end
	else
		[ret,x0,str,ts,xts] = feval(sys);
	end
else
	drawnow % Flash up the model and execute load callback
end
