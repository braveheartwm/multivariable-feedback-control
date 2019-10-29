  	PROGRAM FENSKE
C To compile you typically use a command like: f77 -o fenske fenske.f
C
C	SIMULATES BINARY DISTILLATION COLUMN WITH CONSTANT
C	REALTIVE VOLATILITY (ALFA) AND CONSANT MOLAR FLOWS.
C  The program computes profiles and also steady-state gains etc. 
C
C	S. SKOGESTAD  01 MAY 85
C		      25 OCT 85
C                     12 feb 86

	PARAMETER (NMT=252)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	common /cdist1/ it,x,y,lb,vt,vb,b
	common /cdist2/ alfa,f,zf,yf,xf,qf,nt,nf
	DOUBLE PRECISION X(NMT),Y(NMT)
	DOUBLE PRECISION L,LB,LT,LT0,NMIN,lof
	REAL o1,o2,o3,o4
	real g11,g12,g21,g22,g11s,g12s,g21s,g22s,a,as,rga11,rga11s,cndno,cndnos
	CHARACTER*18 FILENM
	CHARACTER*1 ANSW
	LOGICAL USED
C
	TOL = 1.E-9
	ITMAX = 100
	ipr = -1
C
C
	WRITE(6,1)
1	FORMAT(1H1,//,'                  PROGRAM FENSKE (Feb. 86)'//
     *     '   Binary, Const. rel. volatility, Const. molar flows, ',
     *		 'Total condenser '///)
	WRITE(6,*) 'Enter Feed Data : F, zF, qF(frac.liq.) ?'
	READ(5,*) F,ZF,QF
	if(f.lt.0) then
		f = -f
		ipr = 1
		endif
	WRITE(6,*) 'Enter Product Compositions : yD, xB ?'
	READ(5,*) yD,XB
	WRITE(6,*) 'Enter Realtive Volatility,  ALFA ?'
	READ(5,*) ALFA
C
	call flash
	S = (yD/(1-yD)) / (XB/(1-XB))
	write(6,*) 'Separation factor S =',S
C	DISTILLATE AND BOTTOM PRODUCT RATES FROM TOTAL MASS BALANCE
C
	B = F * (yD-ZF)/(yD-XB)
	D = F * (ZF-XB)/(yD-XB)
	write(6,*) 'Product rates, B= ',b,'   D= ',d
C
C	MIN. NO. OF TRAYS BY FENSKE
C
	NMIN = LOG(S)/LOG(ALFA)
	o1 = nmin
	WRITE(6,*) 'Min no of theor. trays (at total reflux) = ',o1
C
50	WRITE(6,*) 'Enter total no of theor. trays in column,'
     *		,'including reboiler  ?'
	READ(5,*) NTH
	IF(NTH.GT.NMT .OR. NTH.LT.NMIN .or. nth.lt.2) GOTO 50
C   add total condenser...
	NT = NTH + 1
60	WRITE(6,*) 'Enter location of feed tray (reboiler=1), NF ?'
	READ(5,*) NF
	IF(NF.Gt.NT-1) GOTO 60
	IF(NF.Lt.2) GOTO 60
	rest = rrest(nth,yd,xb,d)
C
C
C	ITERATION LOOP FOR TRAY-BY-TRAY CALCULATIONS
C	START FROM TOP AND GUESS REFLUX RATE LT
C	PROCEED DOWNWARDS AND COMPUTE FEED PLATE COMPOSITION 
C       THEN START FROM BOTTOM AND COMPUTE COMPOSITION ON FEED PLATE
C          COMPARE.....
C
C   initial guess...
	LT = REST * D
	icol = 3
100	continue
	llnew = 0
	write(6,*)
	if(ipr.ge.1) then
     	write(6,*) 'it      L       D       yd      xb      xnf(top)   xnf(btm)'
	endif
C
C 	MAIN ITERATION LOOP START HERE
C
	DO 500 IT=1,ITMAX
	  call colit(d,yd,lt,xnftop,xnfbtm,ipr,xb)
C  NOW TRY TO MATCH FEED PLATE COMPOSITIONS...
	  FUNC = XNFTOP - XNFBTM
	  RERR = ABS(FUNC/XNFTOP)
	  IF(RERR.LT.TOL) GOTO 1000
c  update guess for d,yd or lt ......
	  if(icol.eq.1) then
  	     call nr(func,d)
	  else if(icol.eq.2) then
	     ydh = 1 - yd
	     call nr(func,ydh)
	     yd = 1 - ydh
	  else if(icol.eq.3) then
	     call nr(func,lt)
	  endif
500	CONTINUE
 
	write(6,*) 'Not converged....'
	goto 1050
C
1000	CONTINUE
	write(6,*) 'Converged.....  Reflux ratio R =',lt/d
	call out1(d,yd,xb,lt,xavg)
		if(llnew.eq.0) then
		  lt0 = lt
		  vb0 = vb
		  xb0 = xb
		  yd0 = yd
		  d0 = d
		  xavg0 = xavg
		  endif
		  S0 = (yD/(1-yD)) / (XB/(1-XB))
		llnew = 1

990 	print*,'Compute change in lnS with L and V ?'
	read(5,3) answ
	if(answ.eq.'Y' .or. answ.eq.'y') then
c
c  shortcut LV gain matrix...
c
	beta = b/(yd*(1.-yd)*f) + d/(xb*(1.-xb)*f)
	ty = (yd-xb)/(yd*(1.-yd))
	tx = (yd-xb)/(xb*(1.-xb))
	lof = lt/f
	dlnsdls = (nth/2.)*qf/(lof*(lof+qf))
	vof = vb/f
	dlnsdvs = (nth/2.)*(1-qf)/(vof*(vof+(1-qf)))
	g11s = ( tx + b*dlnsdls/f)/beta
	g12s = (-tx + b*dlnsdvs/f)/beta
	g21s = ( ty - d*dlnsdls/f)/beta
	g22s = (-ty - d*dlnsdvs/f)/beta
	as = g12s*g21s/(g11s*g22s)
	rga11s = 1/(1-as)
	cndnos = 2.*(1+abs(as))/abs(1-as)
	
		print*,'Enter size of change in xb (%)?'

		read*,delta
		delta = delta/100.
		xb = xb0*(1.+delta)
		dum = s0*xb/(1-xb)
		ydh = 1/(1+dum)
		d = f * (zf-xb)/((1-ydh)-xb)
		deld = d - d0
		lt = lt0 - deld
		print*,'--- Increase L keeping V constant : ---'
c  change lt and keep vb constant ....
	  DO 1010 IT=1,ITMAX
	  yd = 1 - ydh		
	  call colit(d,yd,lt,xnftop,xnfbtm,ipr,xb)
	  FUNC = XNFTOP - XNFBTM
	  RERR = ABS(FUNC/XNFTOP)
	  IF(RERR.LT.TOL) GOTO 1015
c  update guess for ydh=1-yd ......
	  call nr(func,ydh)
1010	  CONTINUE
	  goto 990
1015		continue	 
		print*
		call out1(d,yd,xb,lt,xavg)
		dlns1 = log((1-xb)*yd/(xb*(1-yd))) - log(s0)
		dl1 = lt - lt0
		dv1 = vb - vb0
		dlp = dl1*100./lt0
		dx1 = xb - xb0
		dy1 = yd - yd0
		print*,'  dl=',dl1,'(=',dlp,'%)'
		o1 = dx1
		o2 = dy1  
		o3 = dlns1
		print*,'  dxb,dyd,dlns=',o1,o2,o3
		lof = lt/f
		dlnsdls = (nth/2.)*qf/(lof*(lof+qf))
		dlnsdl = f*dlns1/dl1
		print*,'  dlnS/d(L/F) =',dlnsdl,'  (shortcut=',dlnsdls,')'
		dlnsdlb = dlnsdl/beta
		print*,'  (1/beta)* dlnS/d(L/F) =',dlnsdlb
		tau = f*(nt-2)*(xavg-xavg0)/(d*dy1+b*dx1)
		print*,'  Estm. time constant (excl.H/F & reb./cond.)=',tau
		WRITE(6,2)
2	        FORMAT(/' HIT RETURN TO CONTINUE..',$)
	        READ(5,3) ANSW
3	        FORMAT(1A1)
		
c   chane vb and keep lt constant ....
		print*,'--- Decrease V keeping L constant : ---'
		lt = lt0
		d = d0 + deld
	  DO 1020 IT=1,ITMAX
       	  yd = 1 - ydh
	  call colit(d,yd,lt,xnftop,xnfbtm,ipr,xb)
	  FUNC = XNFTOP - XNFBTM
	  RERR = ABS(FUNC/XNFTOP)
	  IF(RERR.LT.TOL .and. it.ne.1) GOTO 1025
c  update guess for ydh ......
	  call nr(func,ydh)
1020	  CONTINUE
	  goto 990
1025		continue
		print*
		call out1(d,yd,xb,lt,xavg)
		dlns2 = log((1-xb)*yd/(xb*(1-yd))) - log(s0)
		dl2 = lt - lt0
		dv2 = vb - vb0
		dvp = dv2*100./vb0
		dx2 = xb - xb0
		dy2 = yd - yd0
		print*,'  dv=',dv2,'(=',dvp,'%)'
		o1 = dx2
		o2 = dy2  
		o3 = dlns2
		print*,'  dxb,dyd,dlns=',o1,o2,o3
		vof = vb/f
		dlnsdvs = (nth/2.)*(1-qf)/(vof*(vof+(1-qf)))
		dlnsdv = f*dlns2/dv2
		print*,'  dlnS/d(V/F) =',dlnsdv,'  (shortcut=',dlnsdvs,')'
		dlnsdvb = dlnsdv/beta
		print*,'  (1/beta)* dlnS/d(V/F) =',dlnsdvb
		tau = f*(nt-2)*(xavg-xavg0)/(d*dy2+b*dx2)
		print*,'  Estm. time constant (excl.H/F & reb./cond.)=',tau
		WRITE(6,2)
		READ(5,3) ANSW
		print*,'  LV - configuration gain matrix (y1=yd,u1=L) :'
		print*,'        (shortcut in paranthesis)'
cut		g11 = f*dy1/dl1
cut		g12 = f*dy2/dv2
cut		g21 = f*dx1/dl1
cut		g22 = f*dx2/dv2
		g11 = ( tx + b*dlnsdl/f)/beta
		g12 = (-tx + b*dlnsdv/f)/beta
		g21 = ( ty - d*dlnsdl/f)/beta
		g22 = (-ty - d*dlnsdv/f)/beta
		print*,g11,g12,'    (',g11s,g12s,')'
		print*,g21,g22,'    (',g21s,g22s,')'		
		a = g12*g21/(g11*g22)
		rga11 = 1/(1-a)
		cndno = 2.*(1+abs(a))/abs(1-a)
		print*,'  Rijnsdorp A=',a,'    (',as,')'
		print*,'  RGA11=',rga11,'    (',rga11s,')'
     	 	print*,'  Min. cond. no. =',cndno,'    (',cndnos,')'
	
		WRITE(6,2)
		READ(5,3) ANSW
c   change vb and lt keeping D constant ....
		print*,'--- Increase L keeping D constant : ---'
		lt = lt0 - deld
		d = d0 
	  DO 1030 IT=1,ITMAX
       	  yd = 1 - ydh
	  call colit(d,yd,lt,xnftop,xnfbtm,ipr,xb)
	  FUNC = XNFTOP - XNFBTM
	  RERR = ABS(FUNC/XNFTOP)
	  IF(RERR.LT.TOL .and. it.ne.1) GOTO 1035
c  update guess for ydh ......
	  call nr(func,ydh)
1030	  CONTINUE
	  goto 990
1035		continue
		print*
		call out1(d,yd,xb,lt,xavg)
		dlns3 = log((1-xb)*yd/(xb*(1-yd))) - log(s0)
		dl3 = lt - lt0
		dv3 = vb - vb0
		dvp = dv3*100./vb0
		dx3 = xb - xb0
		dy3 = yd - yd0
		print*,'  dv=',dv3,'(=',dvp,'%)'
		o1 = dx3
		o2 = dy3  
		o3 = dlns3
		print*,'  dxb,dyd,dlns=',o1,o2,o3
		vof = vb/f
		dlnss = dlnsdls + dlnsdvs
		dlnsdv = f*dlns3/dv3
		print*,'  dlnS/d(V/F) =',dlnsdv,'  (shortcut=',dlnss,')'
		dlnsdvb = dlnsdv/beta
		print*,'  (1/beta)* dlnS/d(V/F) =',dlnsdvb
		goto 990
		endif

1050	write(6,*) 'Rerun (Y/N) ?'
	read(5,3) answ
	if(answ.eq.'Y' .or. answ.eq.'y') then
		llnew = 0  
		print*,'Vary d(1),yd(2) or lt(3) or change feed data(4) ?'
		read*,icol
		if(icol.eq.1) then
		 print*,'Enter yd and lt ?'
		 read*,yd,lt
		 d=d0
		else if(icol.eq.2) then
		 print*,'Enter d and lt ?'
	         read*,d,lt
		 yd = yd0
		 ydmx = zf*f/d
		 if(yd.gt.ydmx) yd = ydmx 
		else if(icol.eq.3) then
		 print*,'Enter d and yd ?'
		 read*,d,yd
		 lt = lt0
		else if(icol.eq.4) then
		 print*,'--- Will use unchanged d and lt !! ---'
		 print*,'Enter zf [',zf,']  ?'
		 read*,zf
		 print*,'Enter qf [',qf,']  ?'
		 read*,qf
		 print*,'Enter F  [',f, ']  ?'
		 read*,f
		 if(f.lt.0) then
		  f = -f
		  ipr = -ipr 
		  endif
		 call flash
		 icol = 2
		 d = d0
		 lt = lt0
	         yd = yd0
		 ydmx = zf*f/d
		 if(yd.gt.ydmx) yd = ydmx 
		endif		 
		goto 100
		endif
	WRITE(6,*)
C
C	PRINTOUT
C
	WRITE(6,*) 'Print-out of column profiles (Y/N) ?'
	READ(5,3) ANSW
	IF(ANSW.EQ.'Y' .OR. ANSW.EQ.'y') THEN
	WRITE(6,*) '       TRAY     L...  V...  X... Y...'  	
	DO 2000 N=NT,1,-1
	IF(N.GT.NF) THEN
		o1=VT
		o2=LT
	ELSE IF(N.LE.NF) THEN
		o1=VB
		o2=LB
	ENDIF
	IF(N.EQ.1) o2=f-d
	IF(N.EQ.NT) o1=D
	o3 = x(n)
	o4 = y(n)
	WRITE(6,*) N,o2,o1,o3,o4
2000	CONTINUE
	ENDIF
C
	WRITE(6,*) 'Print-out of column profiles to file (Y/N) ?'
	READ(5,3) ANSW
	IF(ANSW.EQ.'Y' .OR. ANSW.EQ.'y') THEN
	WRITE(6,*) 'Enter Filename ?'
	READ(5,4) FILENM
4	FORMAT(A18)
	INQUIRE(FILE=FILENM,exist=USED)
	IF(USED) THEN
	  OPEN(UNIT=1,FILE=FILENM,STATUS='OLD',FORM='FORMATTED')
	ELSE
          OPEN(UNIT=1,FILE=FILENM,STATUS='NEW',FORM='FORMATTED')
	ENDIF
C
	WRITE(1,*) 'FEED DATA : F,ZF,QF .... :'
	WRITE(1,*) F,ZF,QF
	WRITE(1,*) 'PRODUCT COMPOSITIONS : yD,XB .... :'
	WRITE(1,*) yD,XB
	WRITE(1,*) 'THEOR.TRAYS ,FEED LOC., REL.VOLAT. ; NT,NF,ALFA ....'
	WRITE(1,*)  NT,NF,ALFA
	WRITE(1,*) '       TRAY     L...  V...  X... Y...'  	
	DO 2700 N=NT,1,-1
	IF(N.GT.NF) THEN
		V=VT
		L=LT
	ELSE IF(N.LE.NF) THEN
		V=VB
		L=LB
	ENDIF
	IF(N.EQ.1) L=B
	IF(N.EQ.NT) V=D
	WRITE(1,*) N,L,V,X(N),y(n)
2700	CONTINUE
	REWIND 1
	ENDIF
C
	E N D

c ------------------------------------------------------------------
	subroutine colit(d,yd,lt,xnftop,xnfbtm,ipr,xb)
	PARAMETER (NMT=252)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	common /cdist1/ it,x,y,lb,vt,vb,b
	common /cdist2/ alfa,f,zf,yf,xf,qf,nt,nf
	DOUBLE PRECISION X(NMT),Y(NMT)
	DOUBLE PRECISION LB,LT
	real o1,o2,o3,o4,o5,o6
	EQUIL(XX) = (ALFA*XX) / (1 + (ALFA-1)*abs(XX))
	EQUIL2(YY) = YY / (YY+ALFA*(1-YY))
C   flow rates in top of column...
	VT = LT + D
C   flow rates in bottom part of column...
	FL = QF*F
	FV = F - FL
	LB = LT + FL
	VB = VT - FV
C
C	START CALCULATING FROM TOP......
C
C   total condenser...
	X(NT) = yD
	Y(NT-1) = yD
	X(NT-1) = EQUIL2(Y(NT-1))
C   top section of column...
	if (nf .eq. nt-1) goto 220
	if (nf .eq. nt-2) goto 210
	DO 200 N=NT-2,NF+1,-1
	Y(N) = Y(N+1) + (LT/VT)*(X(N+1)-X(N+2))
	X(N) = EQUIL2(Y(N))
200	CONTINUE
C   feed plate...
210	Y(NF) = (VT*Y(NF+1) + LT*(X(NF+1)-X(NF+2)) - FV*YF) / VB
	X(NF) = EQUIL2(Y(NF))
220	XNFTOP = X(NF)
C
C	NOW START FROM BOTTOM.....
C
C   reboiler (tray1) and tray2...
	b = f - d
	xb = f*zf/b - d*yd/b
	X(1) = XB
	Y(1) = EQUIL(X(1))
	X(2) = (VB*Y(1) + B*X(1)) / LB
C   bottom section...
	if (nf.eq.2) goto 310
	DO 300 N=3,NF
	Y(N-1) = EQUIL(X(N-1))
	X(N) = X(N-1) + (VB/LB)*(Y(N-1)-Y(N-2))
300	CONTINUE
310	XNFBTM = X(NF)
c
	if(ipr.eq.1) then
	o1 = lt
	o5 = d
	o4 = yd
	o6 = xb
	o2 = xnftop
	o3 = xnfbtm
	write(6,*) it,o1,o5,o4,o6,o2,o3
	endif
	return
	end

c-------------------------------------------------------------
	subroutine out1(d,yd,xb,lt,xavg)
	PARAMETER (NMT=252)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	common /cdist1/ it,x,y,lb,vt,vb,b
	common /cdist2/ alfa,f,zf,yf,xf,qf,nt,nf
	DOUBLE PRECISION X(NMT),Y(NMT)
	DOUBLE PRECISION LB,LT
	real o1,o2,o3,o4,o5,o6
	o1 = yd
	o2 = xb
	o3 = d
	o4 = lt
	o5 = vb
	print*,'yd,xb,d,lt,vb=',o1,o2,o3,o4,o5
c  compute avg. composition (excluding reboilerr and condenser)
	xsum = 0.
	do 100 i=2,nt-1
100	xsum = xsum + x(i)
	xavg = xsum/(nt-2)
	o1 = xavg
	s = (1.-xb)*yd/(xb*(1.-yd))
	o2 = log(yd/(1.-yd)) / log(s)
	zed = d*yd*(1.-yd)/f + b*xb*(1.-xb)/f
c  estimate time constant (w/o m/f) exluding reboiler and condenser
	o3 = (nt-2)/(zed*log(s))
	o4 = x(nf)
	print*,'xavg=',o1,'(s.cut=',o2,')  xnf=',o4,
     *          'tau-scut=',o3 
	o6 = log( s )
	o1 = 1. - yd
	o2 = (b/f)/(yd*(1.-yd)) + (d/f)/(xb*(1.-xb))
	print*,'ln S =',o6,'  1-yd =',o1,'  beta =',o2
	print*
	return
	end

c------------------------------------------------------------------

	subroutine nr(func,x)
	implicit double precision(a-h,o-z)
	common /cdist1/ it
	common /csave/ fold,xold
	IF(IT.EQ.1) THEN
		IF(FUNC.GT.0) XNEW = X*1.01
		IF(FUNC.LT.0) XNEW = X*0.99
		GOTO 400
	ENDIF
C
	DERIV = (X-XOLD)/(FUNC-FOLD)
	chg = - deriv*func
	XNEW = X + chg
	if(chg .gt. x) xNEW = x*2
	if(-chg .gt. x/2) xNEW = x*0.5
400	CONTINUE
	XOLD = X
	FOLD = FUNC
	X = XNEW
	return
	end

c-----------------------------------------------------------
	subroutine flash
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	common /cdist2/ alfa,f,zf,yf,xf,qf,nt,nf
C
C	INTERNAL FUNCTIONS FOR VLE
C
	EQUIL(XX) = (ALFA*XX) / (1 + (ALFA-1)*XX)
	EQUIL2(YY) = YY / (YY+ALFA*(1-YY))
C
C	FLASH FEED TO FIND XF AND ZF
C
	IF( ABS(QF) .LT. 1.E-4) THEN
		YF = ZF
		XF = EQUIL2(YF)
	ELSE
		A = QF*(ALFA-1)
		B = ALFA*(1-QF) + QF - (ALFA-1)*ZF
		C = -ZF
		XF = (-B + SQRT(B*B-4*A*C)) / (2*A)
		YF = EQUIL(XF)
	ENDIF
	write(6,*) 'Flash of feed: xF, yF =',xf,yf
	return
	end
C
C ------------------------------------------------------------------

	function rrest(nth,yd,xb,d)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	common /cdist2/ alfa,f,zf,yf,xf,qf,nt,nf
	S = (yD/(1-yD)) / (XB/(1-XB))
C	MIN. REFLUX RATIO BY FENSKE
C
	RMIN = ( yD/XF - ALFA*(1-yD)/(1-XF) ) / (ALFA-1)
	WRITE(6,*) 'Min. reflux ratio L/D = ',RMIN
C
C	ESTIMATED ACTUAL REFLUX RATIO BY JAFAREY  (I&EC PDD,2,200(1979)
C
	B1 = (S**(1.05/NTH)) / ALFA
	BT = 1 - B1*B1
	A = ZF*BT
	B = QF*BT + ZF*BT - 1
	C = -QF*(1-BT)
	REST = (-B + SQRT(B*B-4*A*C)) / (2*A)
	WRITE(6,*) 'Estimated actual reflux ratio (Jafarey,1979)',
     *	 	' L/D = ',REST
	IF(REST.LE.RMIN)  REST = RMIN
	if(rest.le.0) rest = 1.
C
C	ESTIMATED ACTUAL REFLUX RATIO BY Skogestad 
C
	B1 = (S**(1.0/NTH)) / ALFA
	BT = 1 - B1*B1
	A = (D/F)*BT
	B = QF*BT + (d/f)*BT - 1
	C = -QF*(1-BT)
	RSIS = (-B + SQRT(B*B-4*A*C)) / (2*A)
	WRITE(6,*) 'Estimated actual reflux ratio (Skogestad)',
     *	 	' L/D = ',RSIS
C
	WRITE(6,2)
2	FORMAT(/' HIT RETURN TO CONTINUE..',$)
	READ(5,3) ANSW
3	FORMAT(1A1)
	rrest = rest
	return
	end
