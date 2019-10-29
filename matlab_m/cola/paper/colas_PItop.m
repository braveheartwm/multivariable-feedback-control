function [ret,x0,str,ts,xts]=colas_sis(t,x,u,flag);
%COLAS_SIS	is the M-file description of the SIMULINK system named COLAS_SIS.
%	The block-diagram can be displayed by typing: COLAS_SIS.
%
%	SYS=COLAS_SIS(T,X,U,FLAG) returns depending on FLAG certain
%	system values given time point, T, current state vector, X,
%	and input vector, U.
%	FLAG is used to indicate the type of output to be returned in SYS.
%
%	Setting FLAG=1 causes COLAS_SIS to return state derivatives, FLAG=2
%	discrete states, FLAG=3 system outputs and FLAG=4 next sample
%	time. For more information and other options see SFUNC.
%
%	Calling COLAS_SIS with a FLAG of zero:
%	[SIZES]=COLAS_SIS([],[],[],0),  returns a vector, SIZES, which
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

% the system will take on the name of this mfile:
sys = mfilename;
new_system(sys)
simver(1.3)
if (0 == (nargin + nargout))
     set_param(sys,'Location',[124,80,893,602])
     open_system(sys)
end;
set_param(sys,'algorithm',     'Gear')
set_param(sys,'Start time',    '0.0')
set_param(sys,'Stop time',     '500')
set_param(sys,'Min step size', '0.01')
set_param(sys,'Max step size', '5')
set_param(sys,'Relative error','1e-6')
set_param(sys,'Return vars',   '')
set_param(sys,'AssignWideVectorLines','on');

add_block('built-in/Clock',[sys,'/','Clock'])
set_param([sys,'/','Clock'],...
		'position',[85,10,105,30])

add_block('built-in/To Workspace',[sys,'/','time'])
set_param([sys,'/','time'],...
		'mat-name','t',...
		'position',[485,12,535,28])

add_block('built-in/Constant',[sys,'/','F'])
set_param([sys,'/','F'],...
		'orientation',1,...
		'position',[170,65,200,85])

add_block('built-in/Constant',[sys,'/','zF'])
set_param([sys,'/','zF'],...
		'orientation',1,...
		'Value','0.5',...
		'position',[116,65,144,85])

add_block('built-in/Demux',[sys,'/','Demux'])
set_param([sys,'/','Demux'],...
		'outputs','[1,1,1,1,41]',...
		'position',[380,116,435,184])

add_block('built-in/To Workspace',[sys,'/','Comp.'])
set_param([sys,'/','Comp.'],...
		'mat-name','Comp',...
		'buffer','1200',...
		'position',[495,237,545,253])

add_block('built-in/To Workspace',[sys,'/','x_B'])
set_param([sys,'/','x_B'],...
		'mat-name','y2',...
		'buffer','1200',...
		'position',[495,142,545,158])

add_block('built-in/To Workspace',[sys,'/','M_D'])
set_param([sys,'/','M_D'],...
		'mat-name','y3',...
		'buffer','1200',...
		'position',[495,172,545,188])

add_block('built-in/To Workspace',[sys,'/','M_B'])
set_param([sys,'/','M_B'],...
		'mat-name','y4',...
		'buffer','1200',...
		'position',[495,202,545,218])

add_block('built-in/Constant',[sys,'/','rMD'])
set_param([sys,'/','rMD'],...
		'orientation',2,...
		'Value','0.5',...
		'position',[490,300,550,320])

add_block('built-in/Gain',[sys,'/','Gain'])
set_param([sys,'/','Gain'],...
		'orientation',2,...
		'Gain','-10',...
		'position',[330,290,360,320])

add_block('built-in/Sum',[sys,'/','Sum'])
set_param([sys,'/','Sum'],...
		'orientation',2,...
		'inputs','-+',...
		'position',[400,295,420,315])

add_block('built-in/Constant',[sys,'/','rMB'])
set_param([sys,'/','rMB'],...
		'orientation',2,...
		'Value','0.5',...
		'position',[490,350,550,370])

add_block('built-in/Gain',[sys,'/','Gain1'])
set_param([sys,'/','Gain1'],...
		'orientation',2,...
		'Gain','-10',...
		'position',[335,340,365,370])

add_block('built-in/S-Function',[sys,'/',['Distillation',13,'column',13,'(nonlinear)',13,'',13,'LV configuration']])
set_param([sys,'/',['Distillation',13,'column',13,'(nonlinear)',13,'',13,'LV configuration']],...
		'function name','colas',...
		'position',[285,117,355,183])

add_block('built-in/Mux',[sys,'/','Mux'])
set_param([sys,'/','Mux'],...
		'inputs','7',...
		'position',[225,116,265,184])

add_block('built-in/Constant',[sys,'/','qF'])
set_param([sys,'/','qF'],...
		'orientation',1,...
		'position',[220,65,250,85])

add_block('built-in/Sum',[sys,'/','Sum1'])
set_param([sys,'/','Sum1'],...
		'orientation',2,...
		'inputs','-+',...
		'position',[405,345,425,365])

add_block('built-in/Transfer Fcn',[sys,'/','Transfer Fcn'])
set_param([sys,'/','Transfer Fcn'],...
		'orientation',2,...
		'Numerator','[3.6 1]',...
		'Denominator','[3.6 0]',...
		'position',[340,412,375,448])

add_block('built-in/Sum',[sys,'/','Sum2'])
set_param([sys,'/','Sum2'],...
		'orientation',2,...
		'inputs','-+',...
		'position',[415,420,435,440])

add_block('built-in/Constant',[sys,'/','ryD'])
set_param([sys,'/','ryD'],...
		'orientation',2,...
		'Value','0.99',...
		'position',[490,425,550,445])

add_block('built-in/To Workspace',[sys,'/','L'])
set_param([sys,'/','L'],...
		'mat-name','u1',...
		'position',[500,477,550,493])

add_block('built-in/Sum',[sys,'/','Sum3'])
set_param([sys,'/','Sum3'],...
		'orientation',3,...
		'position',[85,185,105,205])

add_block('built-in/Gain',[sys,'/','Gain2'])
set_param([sys,'/','Gain2'],...
		'orientation',2,...
		'Gain','60',...
		'position',[255,415,285,445])

add_block('built-in/To Workspace',[sys,'/','y_D'])
set_param([sys,'/','y_D'],...
		'mat-name','y1',...
		'buffer','1200',...
		'position',[495,72,545,88])


%     Subsystem  'Graph'.

new_system([sys,'/','Graph'])
set_param([sys,'/','Graph'],'Location',[0,59,274,252])

add_block('built-in/Inport',[sys,'/','Graph/x'])
set_param([sys,'/','Graph/x'],...
		'position',[65,55,85,75])

add_block('built-in/S-Function',[sys,'/',['Graph/S-function',13,'M-file which plots',13,'lines',13,'']])
set_param([sys,'/',['Graph/S-function',13,'M-file which plots',13,'lines',13,'']],...
		'function name','sfuny',...
		'parameters','ax, color,dt',...
		'position',[130,55,180,75])
add_line([sys,'/','Graph'],[90,65;125,65])
set_param([sys,'/','Graph'],...
		'Mask Display','plot(0,0,100,100,[90,10,10,10,90,90,10],[65,65,90,40,40,90,90],[90,78,69,54,40,31,25,10],[77,60,48,46,56,75,81,84])',...
		'Mask Type','Graph scope.')
set_param([sys,'/','Graph'],...
		'Mask Dialogue','Graph scope using MATLAB graph window.\nEnter plotting ranges and line type.|Time range:|y-min:|y-max:|Line type (rgbw-:*). Seperate each plot by ''/'':')
set_param([sys,'/','Graph'],...
		'Mask Translate','color = @4; ax = [0, @1, @2, @3]; dt = -1;')
set_param([sys,'/','Graph'],...
		'Mask Help','This block plots to the MATLAB graph window and can be used as an improved version of the Scope block. Look at the m-file sfuny.m to see how it works. This block can take scalar or vector input signal.')
set_param([sys,'/','Graph'],...
		'Mask Entries','20\/0.989\/0.991\/''y-/g--/c-./w:/m*/ro/b+''\/')


%     Finished composite block 'Graph'.

set_param([sys,'/','Graph'],...
		'position',[635,101,665,139])

add_block('built-in/From Workspace',[sys,'/',['From',13,'Workspace']])
set_param([sys,'/',['From',13,'Workspace']],...
		'matl_expr','TL',...
		'position',[20,222,60,248])

add_block('built-in/From Workspace',[sys,'/',['From',13,'Workspace1']])
set_param([sys,'/',['From',13,'Workspace1']],...
		'matl_expr','TV',...
		'position',[20,302,60,328])
add_line(sys,[440,120;450,120;450,80;490,80])
add_line(sys,[440,150;460,150;460,180;490,180])
add_line(sys,[440,165;450,165;450,210;490,210])
add_line(sys,[440,180;439,180;439,245;490,245])
add_line(sys,[185,90;185,100;165,100;165,160;220,160])
add_line(sys,[130,90;130,100;155,100;155,170;220,170])
add_line(sys,[110,20;480,20])
add_line(sys,[360,150;375,150])
add_line(sys,[440,135;475,135;475,150;490,150])
add_line(sys,[270,150;280,150])
add_line(sys,[395,305;365,305])
add_line(sys,[485,310;425,310])
add_line(sys,[325,305;175,305;175,140;220,140])
add_line(sys,[440,150;460,150;460,300;425,300])
add_line(sys,[440,165;450,165;450,350;430,350])
add_line(sys,[485,360;430,360])
add_line(sys,[400,355;370,355])
add_line(sys,[330,355;185,355;185,150;220,150])
add_line(sys,[235,90;235,104;200,104;200,180;220,180])
add_line(sys,[410,430;380,430])
add_line(sys,[485,435;440,435])
add_line(sys,[65,235;90,235;90,210])
add_line(sys,[95,180;95,120;220,120])
add_line(sys,[335,430;290,430])
add_line(sys,[250,430;100,430;100,210])
add_line(sys,[250,430;213,430;213,485;495,485])
add_line(sys,[440,120;630,120])
add_line(sys,[65,315;135,315;135,130;220,130])
add_line(sys,[440,120;581,120;581,400;451,400;451,425;440,425])

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
