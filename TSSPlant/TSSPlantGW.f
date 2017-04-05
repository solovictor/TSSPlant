	program TSSPlantGW

c   May 13,  2015 
c   KAUST 
c
c   Plant Promoter Prediction Program (Genome-wide)  
c
c   Input File: Sequence(s) in FASTA format  
c
c   RUN: TSSPlant -i:<Input sequence(s) file> 
c             ][o:<Output file>][-o:<Output file>] 
c             [-t:<If in region of a given length both TATA and TATA-less promoters
c                         are predicted, exclude TATA-less promoter ("y" or "Y") 
c                    OR  
c                         take TSS closer to Gene start ("n" or "N")> 
c
c   Deault output file: TSSPlantGW.res 
c   Default value for <parT>: y  
c 
c   Data Files required: 
c     (1) ptss_path.dat   
c     (2) Data_TSSPlant subfolder with data diles required 
c         
c
c
c ___________________________________________________

        character*800 path 
         character*300 emp,s,sn,vs(2),ss
         character*300 infile,outfile,dataf,thrf,name 
         character*1 seq(1000),seqt(400000000) 
         character*1 seqc(1000),seqtc(400000000) 	 
  	 character*300 class(2) 
  	 character*300 pr(4),prr(4) 
  	 character*3 pp(4)	 
	 byte seqb(1000),seqtb(400000000),cont(1000000)
	 integer tata(1000000),tatax(1000000),ptype 
	 integer rex,tss(2,1000000),ntss(2),kcls(2)
	 integer prt,h,tatapos(1000000),tataposx 	
	 real frg(4),scortss(2,1000000),scor(16),scorx
	 integer tss1(1000000),tatapos1(1000000),tss2(1000000) 
	 real scortss1(1000000),scortata1(1000000),scortss2(1000000) 	 
	 real pi,scortata(1000000),scortatax,c(2),skx(11,2)
	 
	 integer tatam(1000000),tataxm(1000000) 
	 integer tssm(2,1000000), tataposm(1000000) 
	 real scortssm(2,1000000),scortata1m(1000000) 
	 integer tss1m(1000000),tatapos1m(1000000),tss2m(1000000) 
	 real scortss1m(1000000),scortss2m(1000000) 	 
	 real scortatam(1000000) 
	 integer nstm(1000000),nst2m(1000000)
	  
	 double precision th1,th2(2),dsr(2,3) 
	 double precision w(2,2,30,30),bias(2,2,30),scortssx	 
	 double precision min(2,30),max(2,30)
	 real fr1(4,30),fr1c(4,30) 
	 real fr2(2,4,30),fr2c(2,4,30)	
	 real fr3(2,4,30),fr4(2,4,30)		  
	 real fro3(2,64), fro4a(2,256),fro4b(2,256)
	 real fro6a(2,4096),fro6b(2,4096)
	 real frssk1(2,11),frssk2(2,11)
	 real frdis1(25),frdis2(25),frdis3(2,25)
	 real frre1(2,10000),frre2(2,10000)
	 integer postata(2),ltata(2),nst(1000000)
	 integer posinr(2,2),linr(2,2),nst2(1000000)
	 integer posyp(2,2),lyp(2) 
	 integer posdp(2,2),ldp(2) 
	 integer pos3(2,2)
	 integer pos4a(2,2),pos4b(2,2)
	 integer pos6a(2,2),pos6b(2,2)
	 integer posssk1(2,2),posssk2(2,2)
	 integer posdis1(2),posdis2(22),posdis3(2,2)
	 integer re1l(2,10000),re1m(2,10000)
	 integer re2l(2,10000),re2m(2,10000)
	 integer posre1(2,2),posre2(2,2)
	 integer nsite1(2),nsite2(2),nline(2)
	 character*1 re1seq(2,10000,80),re2seq(2,10000,80)
  	 character*300 hh(4)	 
	 
	 common/posssk1/posssk1
	 common/posssk2/posssk2
	 common/skx/skx
	 common/fr1/fr1
	 common/fr1c/fr1c 
	 common/fr2/fr2
	 common/fr2c/fr2c 
	 common/fr3/fr3
	 common/fr4/fr4 
	 common/nline/nline 
	 common/fro3/fro3
	 common/fro4a/fro4a
	 common/fro4b/fro4b 
	 common/fro6a/fro6a
	 common/fro6b/fro6b 
	 common/frssk1/frssk1
	 common/frssk2/frssk2 
	 common/frdis1/frdis1
	 common/nsite1/nsite1
	 common/nsite2/nsite2
	 common/frdis2/frdis2
	 common/frdis3/frdis3 
	 common/frre1/frre1
	 common/frre2/frre2 
	 common/postata/postata
	 common/ltata/ltata 
	 common/posinr/posinr
	 common/linr/linr 
	 common/posyp/posyp
	 common/lyp/lyp
	 common/posdp/posdp	 
	 common/ldp/ldp 
	 common/pos3/pos3
	 common/pos4a/pos4a
	 common/pos4b/pos4b 
	 common/pos6a/pos6a
	 common/pos6b/pos6b 
	 common/posdis1/posdis1
	 common/posdis2/posdis2
	 common/posdis3/posdis3 
	 common/re1l/re1l
	 common/re2l/re2l
	 common/re1m/re1m
	 common/re2m/re2m	 
 	 common/posre1/posre1
	 common/posre2/posre2
	 common/re1seq/re1seq 
	 common/re2seq/re2seq 
	 common/seq/seq
	 common/seqc/seqc	 
	 common/seqb/seqb
	 common/kpos/kpos
	 common/scorx/scorx
	 common/lseq/lseq	
	 common/frg/frg
	 common/emp/emp
	 common/dsr/dsr 
 	 common/minmax/min,max
	 common/iftx/iftx
	 common/th/th1,th2
	 common/kcls/kcls
         common/envp/path,lpath 
	 common/nf/nf
	 common/scor/scor
	 common/ipred/ipred
	 common/scortssx/scortssx
	 common/pi/pi	
	 common/c/c 
	 common/w/w
	 common/bias/bias
	 common/kint/kint
	 common/itss/itss	 
c  ___________________________________________________

       do i=1,300
          emp(i:i)=' '
       enddo
       do i=1,800
          path(i:i)=' '
       enddo
       	
        call getenv('TSSPlant_DATA',path)
        do i=800,1,-1
           if(path(i:i).ne.' ')goto 8182
        enddo
 8182     if(i.eq.0)then
        print *, ' '   
        print *, 'Please setup environmental variable ',
     *	         'TSSPlant_DATA with'   
        print *, 'the name of directory where data files ',
     *	         'are located'   
        print *,' For example:'
        print *,' setenv TSSPlant_DATA /...directory',
     *          ' with TSSPlant data files required'
        print *,' export TSSPlant_DATA="directory',
     *          ' with TSSPlant data files required"'	
        print *,' '
        stop
          endif

        if(path(i:i).ne.'/')then
	  	i=i+1
	  	path(i:i)='/'
	endif
	lpath=i
       
       
	data pp/'-i:','-o:','-t:','-l:'/
	lseq=251	!Length of Window to compute Features scores
	thrf=emp
	
	ntot1=0
	ntot2=0
	ntot3=0
	ntot4=0
	ntot5=0
	ltot=0

	ngene=0
	do i=1,100000
		nst(i)=-100000
	enddo	
       	
	itss=201
	lex=15		! Left  extention of INR search region 
	rex=15          ! Right extention of INR search region  
        pi=3.1415927
	kint=25
	kcls(1)=1	! TATA+
	kcls(2)=2	! TATA-
	 
	do i=1,2
		vs(i)=emp
	enddo
	
	skx(1,1)=-1.0
	skx(1,2)=0.0		
	skx(2,1)=0.0
	skx(2,2)=0.1		
	skx(3,1)=0.1
	skx(3,2)=0.2		
	skx(4,1)=0.2
	skx(4,2)=0.3		
	skx(5,1)=0.3
	skx(5,2)=0.4		
	skx(6,1)=0.4
	skx(6,2)=0.5		
	skx(7,1)=0.5
	skx(7,2)=0.6		
	skx(8,1)=0.6
	skx(8,2)=0.7		
	skx(9,1)=0.7
	skx(9,2)=0.8	
	skx(10,1)=0.8
	skx(10,2)=0.9		
	skx(11,1)=0.9
	skx(11,2)=1.0	
	  		    
 1	format(a) 
 2	format(2i5)      
 3	format(i5)
 4	format(i2) 
 5	format(f12.4)  
 6	format('Program TSSPlant',/,'Search for RNA II promoters (TSSs)',/,
     *         'Input file with query sequence(s): ',a)
 7	format(//,'Query: ',a,/,'Length of Query sequence: ',i6) 
 8	format(60a1) 
 9	format('Thresholds, for TATA      promoters: ',f6.2,/, 
     *         '                TATA-less promoters: ',f6.2)
 6699	format('Annotated Gene Start Position: ',i6) 
 6527	format(4i5)  
 6533   format(f12.6)  
 6554	format(30(f7.3)) 
 6555 	format(6f12.6) 
 6556	format(2i6)
 6557	format(f9.6)
 3500	format('Promoters not found')
 3501	format(i6,' promoter(s) predicted:')	
 3502	format('TSS position: ',i6,3x,'TSS score = ',f8.4,5x,
     *          'TATA-box position: ',i6,3x,'TATA-box score = ',f8.4)
 3503	format('(TATA-) TSS position: ',i6,3x,'TSS score = ',f8.4)  
 8778	format('If in interval -/+',i5,' nt around putative TSS ', 
     *   ' both TATA and TATA-less promoters ',/, 
     *   '   are predicted, TATA-less promoter is excluded') 	
 8779	format('If in interval -/+',i5,' nt around putative TSS ',
     *   ' both TATA and TATA-less promoters ',/,
     *   '   are predicted, promoter (TSS) closer to ',
     *   'Gene start is taken')  	

         
c ______ Get parameters ______ 
	
	do i=1,4
		pr(i)=emp
		prr(i)=emp
	enddo

	iargno=iargc()
	if((iargno.lt.1).or.(iargno.gt.6))then 
	    call help
            stop
	endif
	
	do i=1,iargno
		call getarg(i,prr(i)) 	
	enddo

			do 1045 i=1,iargno
		do j=1,4
	if(prr(i)(1:3).eq.pp(j))then
		pr(j)(1:297)=prr(i)(4:300)
		goto 1045
	endif
		enddo
 1045			continue

	if(pr(1).eq.emp)then 
	    call help
            stop
	endif
	infile=emp
	infile=pr(1)	
	
	dataf=emp
	dataf(1:12)='TSSPlant.dat'
	frg(1)=0.26915
	frg(2)=0.20302
	frg(3)=0.22121
	frg(4)=0.30661
	
	thrf(1:22)='thr_tata_inr.dat'		
	vs(1)(1:14)='visan_tata.dat'
 	vs(2)(1:18)='visan_tataless.dat'
			
 	ptype=1
 	if((pr(3)(1:1).eq.'n').or.(pr(3)(1:1).eq.'N'))ptype=0

 	outfile=emp
	if(pr(2).eq.emp)then
		outfile='TSSPlantGW.res'
	else
		outfile=pr(2)
	endif 
	
	if(pr(4).eq.emp)then
		inter=300
		goto 7891
	endif
	do i=1,300
		if(pr(4)(i:i).ne.' ')goto 7987
	enddo
 7987	k1=i
 	k2=len1(pr(4))		
	inter=intr(pr(4),k1,k2)
	if(inter.lt.100)inter=300	
	
c ______ Get Thresholds (for TATA and INR elements), Input files and Search Parameters 

 7891 	call openf(1,thrf) 
	read(1,1)
	read(1,*)th1
	read(1,1)
	read(1,*)th2(1)
	read(1,*)th2(2) 		
 3373 	close(1)
	
c ============================================ 

	call openf(2,dataf) 

c __________ Features for TATA+ promoters
	
 13	s=emp
	read(2,1,end=1899)s
	if(s(1:6).ne.'Class:')goto 13
	class(1)=emp	 
	class(1)(1:5)='TATA+'
	
	read(2,4)nf1 
	read(2,1)
	
c Feature: TATA  

	read(2,1)
	do i=1,2
		hh(i)=emp
		read(2,1,end=1899)hh(i)
	enddo
	read(2,2)postata(1),postata(2) 
	read(2,2)ltata(1),ltata(2) 

	call openf(31,hh(1))
	do i=1,4
		read(31,1)
	enddo	
	do i=1,4
		read(31,6554)(fr1(i,j),j=1,ltata(2))
	enddo
	close(31)
	
	call openf(31,hh(2))
	do i=1,4
		read(31,1)
	enddo	
	do i=1,4
		read(31,6554)(fr1c(i,j),j=1,ltata(1))
	enddo
	close(31)

c Feature: INR (TATA+)

	read(2,1)	
	do i=1,2
		hh(i)=emp
		read(2,1,end=1899)hh(i)
	enddo
	read(2,2)posinr(1,1),posinr(1,2) 
	read(2,2)linr(1,1),linr(1,2)
	posinr(1,1)=posinr(1,1)-lex
	posinr(1,2)=posinr(1,2)+rex 

	call openf(31,hh(1))
	do i=1,4
		read(31,1)
	enddo	
	do i=1,4
		read(31,6554)(fr2(1,i,j),j=1,linr(1,2)) 
	enddo
	close(31)
	
	call openf(31,hh(2))
	do i=1,4
		read(31,1)
	enddo	
	do i=1,4
		read(31,6554)(fr2c(1,i,j),j=1,linr(1,1)) 
	enddo
	close(31)	

c Feature: YP  (TATA+)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)
	read(2,2)posyp(1,1),posyp(1,2) 
	read(2,3)lyp(1)  

	call openf(31,hh(1))
	do i=1,4
		read(31,1)
	enddo	
	do i=1,4
		read(31,6554)(fr3(1,i,j),j=1,lyp(1)) 
	enddo
	close(31)

c Feature: DPE  (TATA+)

 	read(2,1)
 	hh(1)=emp
	read(2,1,end=1899)hh(1)
	read(2,2)posdp(1,1),posdp(1,2) 
	read(2,3)ldp(1)  

	call openf(31,hh(1))
	do i=1,4
		read(31,1)
	enddo	
	do i=1,4
		read(31,6554)(fr4(1,i,j),j=1,ldp(1)) 
	enddo
	close(31)
	
c Feature: Triplets (TATA+)

 	read(2,1)
 	hh(1)=emp
	read(2,1,end=1899)hh(1)	
	read(2,2)pos3(1,1),pos3(1,2) 
	
	call openf(31,hh(1))
	do i=1,2
		read(31,1)
	enddo		
	do i=1,64,6
		i1=i
		i2=i1+5
		if(i2.gt.64)i2=64
		read(31,6555)(fro3(1,ii),ii=i1,i2)
	enddo
	close(31)	
	
c Feature: Tetra-1 (TATA+)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)pos4a(1,1),pos4a(1,2) 
	
	call openf(31,hh(1))
	do i=1,2
		read(31,1)
	enddo		
	do i=1,256,6
		i1=i
		i2=i1+5
		if(i2.gt.256)i2=256
		read(31,6555)(fro4a(1,ii),ii=i1,i2)
	enddo
	close(31)

c Feature: Tetra-2 (TATA+)

	read(2,1) 
 	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)pos4b(1,1),pos4b(1,2) 
	
	call openf(31,hh(1))
	do i=1,2
		read(31,1)
	enddo		
	do i=1,256,6
		i1=i
		i2=i1+5
		if(i2.gt.256)i2=256
		read(31,6555)(fro4b(1,ii),ii=i1,i2)
	enddo
	close(31)	 

c Feature: Hexa-1 (TATA+)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)pos6a(1,1),pos6a(1,2) 
	
	call openf(31,hh(1))
	do i=1,2
		read(31,1)
	enddo		
	do i=1,4096,6
		i1=i
		i2=i1+5
		if(i2.gt.4096)i2=4096
		read(31,6555)(fro6a(1,ii),ii=i1,i2)
	enddo
	close(31)	 
	
c Feature: Hexa-2 (TATA+)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)pos6b(1,1),pos6b(1,2) 
	
	call openf(31,hh(1))
	do i=1,2
		read(31,1)
	enddo		
	do i=1,4096,6
		i1=i
		i2=i1+5
		if(i2.gt.4096)i2=4096
		read(31,6555)(fro6b(1,ii),ii=i1,i2)
	enddo
	close(31)	
	
c Feature: CG-scew (TATA+)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)posssk1(1,1),posssk1(1,2) 
	
	call openf(31,hh(1))
	do i=1,3
		read(31,1)
	enddo		
	do i=1,11,6
		i1=i
		i2=i1+5
		if(i2.gt.11)i2=11
		read(31,6555)(frssk1(1,ii),ii=i1,i2)
	enddo	 
	close(31)		 
 
c Feature: AC-skew (TATA+)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)posssk2(1,1),posssk2(1,2) 
	
	call openf(31,hh(1))
	do i=1,3
		read(31,1)
	enddo		
	do i=1,11,6
		i1=i
		i2=i1+5
		if(i2.gt.11)i2=11
		read(31,6555)(frssk2(1,ii),ii=i1,i2)
	enddo	 
	close(31) 

c Feature: TATA-TSS distance (TATA+)

	read(2,1)
	do i=1,2
		hh(i)=emp
		read(2,1,end=1899)hh(i)	
	enddo		 

	call openf(31,hh(2))
	do i=1,3
		read(31,1)
	enddo		
	do i=1,kint,6
		i1=i
		i2=i1+5
		if(i2.gt.kint)i2=kint
		read(31,6555)(frdis1(ii),ii=i1,i2)
	enddo	 
	close(31)
		
  	call openf(13,hh(1))
 	read(13,1)
	read(13,1)
	read(13,6533)dsr(1,1)
  	close(13)	
	
c Feature: TATA-INR distance (TATA+)

	read(2,1)
	do i=1,2
		hh(i)=emp
		read(2,1,end=1899)hh(i)	
	enddo		
	
	call openf(31,hh(2))
	do i=1,3
		read(31,1)
	enddo		
	do i=1,kint,6
		i1=i
		i2=i1+5
		if(i2.gt.kint)i2=kint
		read(31,6555)(frdis2(ii),ii=i1,i2)
	enddo	 
	close(31)
	
  	call openf(13,hh(1))	
	do j=1,4
		read(13,1)
	enddo	
	read(13,6533)dsr(1,2)
  	close(13) 
		
c Feature: INR-TSS distance (TATA+)

	read(2,1)
	do i=1,2
		hh(i)=emp
		read(2,1,end=1899)hh(i)	
	enddo		
	
	call openf(31,hh(2))
	do i=1,3
		read(31,1)
	enddo		
	do i=1,kint,6
		i1=i
		i2=i1+5
		if(i2.gt.kint)i2=kint
		read(31,6555)(frdis3(1,ii),ii=i1,i2)
	enddo	 
	close(31)
	
  	call openf(13,hh(1))	
	do j=1,6
		read(13,1)
	enddo	
	read(13,6533)dsr(1,3)
 	close(13) 
	
c Feature: SRED-1 (TATA+)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)	
	hh(2)=emp
	read(2,1,end=1899)hh(2)			
	read(2,2)posre1(1,1),posre1(1,2) 
	
	call openf(21,hh(1))
	read(21,1)
	read(21,3)nsite1(1)		
	do i=1,nsite1(1)
		read(21,6557)frre1(1,i) 
	enddo
	close(21)
	
	call openf(21,hh(2))
	read(21,1)
	read(21,1)
	do i=1,nsite1(1)
		read(21,1)
		s=emp
		read(21,1)s
		read(21,6556)re1l(1,i),re1m(1,i) 
		do j=1,re1l(1,i) 
			re1seq(1,i,j)=s(j:j)			
		enddo
	enddo	
	close(21)		

c Feature: SRED-2 (TATA+)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)	
	hh(2)=emp
	read(2,1,end=1899)hh(2)			
	read(2,2)posre2(1,1),posre2(1,2) 
	
	call openf(21,hh(1))
	read(21,1)
	read(21,3)nsite2(1)		
	do i=1,nsite2(1)
		read(21,6557)frre2(1,i) 
	enddo	
	close(21)
	
	call openf(21,hh(2))
	read(21,1)
	read(21,1)
	do i=1,nsite2(1)
		read(21,1)
		s=emp
		read(21,1)s
		read(21,6556)re2l(1,i),re2m(1,i) 
		do j=1,re2l(1,i) 
			re2seq(1,i,j)=s(j:j)
		enddo	 
	enddo	
	close(21)
	
c  Visan Learning parameters  
 
 	call openf(35,vs(1)) 
	read(35,4)nline(1) 
	read(35,*)c(1)	!Visan threshold
	read(35,*)(min(1,i),i=1,nf1)	!Weights Minimum 
	read(35,*)(max(1,i),i=1,nf1)	!Weights Maximum 
	h=1
	read(35,1) 
		do j=1,nline(1)
	read(35,*)(w(1,h,j,i),i=1,nf1)
		enddo
		
	h=2
	read(35,1)
		do j=1,2
	read(35,*)(w(1,h,j,i),i=1,nline(1))	
		enddo
			
	h=1
	read(35,1)
	read(35,1)
	do j=1,nline(1)			
		read(35,*)bias(1,h,j)
	enddo	

	h=2
	read(35,1)
	do j=1,2			
		read(35,*)bias(1,h,j)
	enddo	
 		
	close(35)
		 
c __________ Features for TATA- promoters
	
 113	s=emp
	read(2,1,end=1899)s
	if(s(1:6).ne.'Class:')goto 113
	class(2)=emp	 
	class(2)(1:5)='TATA-'
	
	read(2,4)nf2
	read(2,1) 
	
c Feature: INR (TATA-)

	read(2,1)
	do i=1,2
		hh(i)=emp
		read(2,1,end=1899)hh(i)
	enddo
	read(2,2)posinr(2,1),posinr(2,2) 
	read(2,2)linr(2,1),linr(2,2) 
	posinr(2,1)=posinr(2,1)-lex
	posinr(2,2)=posinr(2,2)+rex 
	

	call openf(31,hh(1))
	do i=1,4
		read(31,1)
	enddo	
	do i=1,4
		read(31,6554)(fr2(2,i,j),j=1,linr(2,2)) 
	enddo
	close(31)
	
	call openf(31,hh(2))
	do i=1,4
		read(31,1)
	enddo	
	do i=1,4
		read(31,6554)(fr2c(2,i,j),j=1,linr(2,1)) 
	enddo
	close(31)	

c Feature: YP  (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)
	read(2,2)posyp(2,1),posyp(2,2) 
	read(2,3)lyp(2)  

	call openf(31,hh(1))
	do i=1,4
		read(31,1)
	enddo	
	do i=1,4
		read(31,6554)(fr3(2,i,j),j=1,lyp(2)) 
	enddo
	close(31)

c Feature: DPE  (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)
	read(2,2)posdp(2,1),posdp(2,2) 
	read(2,3)ldp(2)  

	call openf(31,hh(1))
	do i=1,4
		read(31,1)
	enddo	
	do i=1,4
		read(31,6554)(fr4(2,i,j),j=1,ldp(2)) 
	enddo
	close(31)
	
c Feature: Triplets (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)	
	read(2,2)pos3(2,1),pos3(2,2) 
	
	call openf(31,hh(1))
	do i=1,2
		read(31,1)
	enddo		
	do i=1,64,6
		i1=i
		i2=i1+5
		if(i2.gt.64)i2=64
		read(31,6555)(fro3(2,ii),ii=i1,i2)
	enddo
	close(31)	 
	
c Feature: Tetra-1 (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)pos4a(2,1),pos4a(2,2) 
	
	call openf(31,hh(1))
	do i=1,2
		read(31,1)
	enddo		
	do i=1,256,6
		i1=i
		i2=i1+5
		if(i2.gt.256)i2=256
		read(31,6555)(fro4a(2,ii),ii=i1,i2)
	enddo
	close(31)	 
	
c Feature: Tetra-2 (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)pos4b(2,1),pos4b(2,2) 
	
	call openf(31,hh(1))
	do i=1,2
		read(31,1)
	enddo		
	do i=1,256,6
		i1=i
		i2=i1+5
		if(i2.gt.256)i2=256
		read(31,6555)(fro4b(2,ii),ii=i1,i2)
	enddo
	close(31)	 

c Feature: Hexa-1 (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)pos6a(2,1),pos6a(2,2) 
	
	call openf(31,hh(1))
	do i=1,2
		read(31,1)
	enddo		
	do i=1,4096,6
		i1=i
		i2=i1+5
		if(i2.gt.4096)i2=4096
		read(31,6555)(fro6a(2,ii),ii=i1,i2)
	enddo
	close(31)	 

c Feature: Hexa-2 (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)pos6b(2,1),pos6b(2,2) 
	
	call openf(31,hh(1))
	do i=1,2
		read(31,1)
	enddo		
	do i=1,4096,6
		i1=i
		i2=i1+5
		if(i2.gt.4096)i2=4096
		read(31,6555)(fro6b(2,ii),ii=i1,i2)
	enddo
	close(31)	
	
c Feature: CG-skew (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)posssk1(2,1),posssk1(2,2) 
	
	call openf(31,hh(1))
	do i=1,3
		read(31,1)
	enddo		
	do i=1,11,6
		i1=i
		i2=i1+5
		if(i2.gt.11)i2=11
		read(31,6555)(frssk1(2,ii),ii=i1,i2)
	enddo	 
	close(31)	 
 
c Feature: AC-skew (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)		
	read(2,2)posssk2(2,1),posssk2(2,2) 
	
	call openf(31,hh(1))
	do i=1,3
		read(31,1)
	enddo		
	do i=1,11,6
		i1=i
		i2=i1+5
		if(i2.gt.11)i2=11
		read(31,6555)(frssk2(2,ii),ii=i1,i2)
	enddo	 
	close(31)	 
		
c Feature: INR-TSS distance (TATA-)

	read(2,1)
	do i=1,2
		hh(i)=emp
		read(2,1,end=1899)hh(i)	
	enddo		

	call openf(31,hh(2))
	do i=1,3
		read(31,1)
	enddo		
	do i=1,kint,6
		i1=i
		i2=i1+5
		if(i2.gt.kint)i2=kint
		read(31,6555)(frdis3(2,ii),ii=i1,i2)
	enddo	 
	close(31)
	
  	call openf(13,hh(1))
	do j=1,8
		read(13,1)
	enddo	
	read(13,6533)dsr(2,3)
 	close(13) 
	
	
c Feature: SRED-1 (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)	
	hh(2)=emp
	read(2,1,end=1899)hh(2)			
	read(2,2)posre1(2,1),posre1(2,2) 
	
	call openf(21,hh(1))
	read(21,1)
	read(21,3)nsite1(2)		
	do i=1,nsite1(2)
		read(21,6557)frre1(2,i) 
	enddo	
	close(21)
	
	call openf(21,hh(2))
	read(21,1)
	read(21,1)
	do i=1,nsite1(2)
		read(21,1)
		s=emp
		read(21,1)s
		read(21,6556)re1l(2,i),re1m(2,i) 
		do j=1,re1l(2,i) 
			re1seq(2,i,j)=s(j:j)
		enddo	 
	enddo	
	close(21)	

c Feature: SRED-2 (TATA-)

	read(2,1)
	hh(1)=emp
	read(2,1,end=1899)hh(1)	
	hh(2)=emp
	read(2,1,end=1899)hh(2)			
	read(2,2)posre2(2,1),posre2(2,2) 
	
	call openf(21,hh(1))
	read(21,1)
	read(21,3)nsite2(2)		
	do i=1,nsite2(2)
		read(21,6557)frre2(2,i) 
	enddo	
	close(21)
	
	call openf(21,hh(2))
	read(21,1)
	read(21,1)
	do i=1,nsite2(2)
		read(21,1)
		s=emp
		read(21,1)s
		read(21,6556)re2l(2,i),re2m(2,i) 
		do j=1,re2l(2,i) 
			re2seq(2,i,j)=s(j:j)
		enddo	 
	enddo	
	close(21)
	
c  Visan Learning parameters  
 
 	call openf(35,vs(2)) 
	read(35,4)nline(2)
	read(35,*)c(2)	!Visan threshold
	read(35,*)(min(2,i),i=1,nf2)	!Weights Minimum 
	read(35,*)(max(2,i),i=1,nf2)	!Weights Maximum 
	h=1
	read(35,1) 
		do j=1,nline(2) 
	read(35,*)(w(2,h,j,i),i=1,nf2)
		enddo
		
	h=2
	read(35,1)
		do j=1,2
	read(35,*)(w(2,h,j,i),i=1,nline(2))
		enddo
			
	h=1
	read(35,1)
	read(35,1)
	do j=1,nline(2)			
		read(35,*)bias(2,h,j)
	enddo	

	h=2
	read(35,1)
	do j=1,2			
		read(35,*)bias(2,h,j)
	enddo	
 		
	close(35)
	
			
 1899	close(2)
 
	open(3,file=outfile,status='unknown')	
	write(3,6)infile(1:len1(infile))
	write(3,9)c(1),c(2)
	if(ptype.eq.1)write(3,8778)inter
	if(ptype.eq.0)write(3,8779)inter	
	
c ________________ Read Query sequence start _____________________

	open(1,file=infile,status='old')
	iend=0
 10	s=emp
	read(1,1,end=88)s
	if(s(1:1).ne.'>')goto 10
	ngene=ngene+1
	name=emp
	name=s
	lseqt=0
		
	do i=1,2
		ntss(i)=0
	enddo	
	
 12	s=emp
	read(1,1,end=8888)s
	if(s(1:1).eq.'>')then
		sn=emp
		sn=s
		goto 1112
	endif	
  	do i=1,len1(s)
		if(s(i:i).ne.' ')then
			lseqt=lseqt+1
			seqt(lseqt)=s(i:i)
		endif
	enddo
	goto 12
 8888	iend=1
 1112	if(lseqt.lt.251)then
 		if(iend.eq.1)goto 88
   		name=emp
 		name=sn 
		lseqt=0	
		ngene=ngene+1
		goto 12
	endif	
	
			do 1114 j=1,lseqt
			
	if((seqt(j).eq.'a').or.
     *	   (seqt(j).eq.'A'))then
     		seqtb(j)=1
     	        seqt(j)='a'
		goto 1114
	endif
	if((seqt(j).eq.'c').or.
     *	   (seqt(j).eq.'C'))then
     		seqtb(j)=2
     	        seqt(j)='c'
		goto 1114		
	endif
	if((seqt(j).eq.'g').or.
     *	   (seqt(j).eq.'G'))then
     		seqtb(j)=3
     	        seqt(j)='g'
		goto 1114		
	endif	
	if((seqt(j).eq.'t').or.
     *	   (seqt(j).eq.'T'))then
     		seqtb(j)=4
     	        seqt(j)='t'
		goto 1114		
	endif	
	seqtb(j)=2
	seqt(j)='c'					
		
 1114			continue
 
 	lq=0		
	do j=lseqt,1,-1
		lq=lq+1
		seqtc(lq)='g'
		if(seqt(j).eq.'a')seqtc(lq)='t'
		if(seqt(j).eq.'c')seqtc(lq)='g'
		if(seqt(j).eq.'g')seqtc(lq)='c'
		if(seqt(j).eq.'t')seqtc(lq)='a'
	enddo
			
	if(mtss.eq.0)mtss=lseqt
	if(mtss.gt.lseqt)mtss=lseqt			
	write(3,7)name(1:len1(name)),lseqt
	write(3,6699)mtss	
	if(prt.eq.1)then
		write(3,1)
		do i=1,lseqt,60
			i1=i
			i2=i1+59
			if(i2.gt.lseqt)i2=lseqt
			write(3,8)(seqt(j),j=i1,i2)
		enddo
		write(3,1)
	endif	
	

c ________ Compute scores and predict TSSs _______

	ltot=ltot+lseqt

	ntss(1)=0
	ntss(2)=0
	    
		do 1212 ir=1,lseqt-250
	i1=ir
	i2=i1+250
	ktss=i1+200
	lseq=0
	do i=i1,i2
		lseq=lseq+1
		seq(lseq)=seqt(i) 
		seqb(lseq)=seqtb(i)
	enddo

	j1=lseqt-i2+1
	j2=j1+250
	lq=0
	do j=j1,j2
		lq=lq+1
		seqc(lq)=seqtc(j) 
	enddo	

c ......... Search for TATA+ promoter 
 
 	call stata
	if(scorx.lt.th1)goto 1222
	scor(1)=scorx
	tataposx=ktss+kpos-200
	scortatax=scorx
  
 	call sinr(kcls(1)) 
	scor(2)=scorx

 	call syp(kcls(1)) 
	scor(3)=scorx
	
 	call sdp(kcls(1)) 
	scor(4)=scorx
	
 	call strip(kcls(1)) 
	scor(5)=scorx
	
 	call stetra1(kcls(1)) 
	scor(6)=scorx

 	call stetra2(kcls(1)) 
	scor(7)=scorx
 
 	call shexa1(kcls(1)) 
	scor(8)=scorx

 	call shexa2(kcls(1))
	scor(9)=scorx
	
 	call ssk1(kcls(1)) 
	scor(10)=scorx

 	call ssk2(kcls(1))
	scor(11)=scorx
	
 	call sdis1 
	scor(12)=scorx

 	call sdis2(kcls(1)) 
	scor(13)=scorx

 	call sdis3(kcls(1))
	scor(14)=scorx
	
 	call sred1(kcls(1))
	scor(15)=scorx

 	call sred2(kcls(1))
	scor(16)=scorx
	 
c ______ call VISAN ______ 
 
	call visanx(kcls(1),nf1)
	if(ipred.eq.1)then
 		ntss(1)=ntss(1)+1
 		tss(1,ntss(1))=ktss
 		scortss(1,ntss(1))=scortssx
		tatapos(ntss(1))=tataposx 
		scortata(ntss(1))=scortatax
	endif

	
c ......... Search for TATA- promoter 
 
 1222	call sinr(kcls(2)) 
	scor(1)=scorx

 	call syp(kcls(2)) 
	scor(2)=scorx
	
 	call sdp(kcls(2)) 
	scor(3)=scorx
	
 	call strip(kcls(2)) 
	scor(4)=scorx
	
 	call stetra1(kcls(2)) 
	scor(5)=scorx

 	call stetra2(kcls(2)) 
	scor(6)=scorx
 
 	call shexa1(kcls(2)) 
	scor(7)=scorx

 	call shexa2(kcls(2))
	scor(8)=scorx
	
 	call ssk1(kcls(2)) 
	scor(9)=scorx

 	call ssk2(kcls(2))
	scor(10)=scorx
	
 	call sdis3(kcls(2)) 
	scor(11)=scorx

 	call sred1(kcls(2))
	scor(12)=scorx

 	call sred2(kcls(2))
	scor(13)=scorx
	  
c ______ call VISAN ______ 
 
	call visanx(kcls(2),nf2)
	if(ipred.eq.1)then
 		ntss(2)=ntss(2)+1
 		tss(2,ntss(2))=ktss		
 		scortss(2,ntss(2))=scortssx
	endif

 1212		        continue
 
c Sorting TATA promoters	

	nt1=0			
	if(ntss(1).eq.0)goto 1001 
	do i=1,ntss(1) 
		cont(i)=1
	enddo
	if(ntss(1).eq.1)then
		nt1=1
		goto 5432
	endif	 

 	do i=1,ntss(1)
	do j=i+1,ntss(1)	
		if(scortss(1,i).lt.scortss(1,j))then 
			sk1=scortss(1,i)
			sk2=scortata(i)
			k1=tss(1,i)
			k2=tatapos(i)
			scortss(1,i)=scortss(1,j)
			scortata(i)=scortata(j)
			tss(1,i)=tss(1,j)
			tatapos(i)=tatapos(j)
			scortss(1,j)=sk1
			scortata(j)=sk2 
			tss(1,j)=k1 
			tatapos(j)=k2 
		endif
	enddo
	enddo
		
	do 261 i=1,ntss(1)
	do 2261 j=i+1,ntss(1)
		if(cont(j).eq.0)goto 2261
		if(abs(tss(1,i)-tss(1,j)).le.inter)cont(j)=0
 2261	continue 
 261	continue
 
 	do 963 i=1,ntss(1)
		if(cont(i).eq.0)goto 963
		nt1=nt1+1
		tss1(nt1)=tss(1,i)
		tatapos1(nt1)=tatapos(i)
		scortss1(nt1)=scortss(1,i)
		scortata1(nt1)=scortata(i) 
 963	continue
 
 5432	if(nt1.eq.1)then 
		tss1(nt1)=tss(1,1)
		tatapos1(nt1)=tatapos(1)
		scortss1(nt1)=scortss(1,1)
		scortata1(nt1)=scortata(1) 
		goto 1001
  	endif

 	do i=1,nt1 
	do j=i+1,nt1 
		if(tss1(i).lt.tss1(j))then 
			sk1=scortss1(i)
			sk2=scortata(i)
			k1=tss1(i)
			k2=tatapos(i)
			icontx=cont(i) 
			scortss1(i)=scortss1(j)
			scortata(i)=scortata(j)
			tss1(i)=tss1(j)
			tatapos(i)=tatapos(j)
			cont(i)=cont(j)
			scortss1(j)=sk1
			scortata(j)=sk2 
			tss1(j)=k1 
			tatapos(j)=k2 
			cont(j)=icontx
		endif
	enddo
	enddo
 	
c Sorting TATA-less promoters	

 1001	nt2=0
 	if(ntss(2).eq.0)goto 2002
  	do i=1,ntss(2) 
		cont(i)=1
	enddo
	if(ntss(2).eq.1)goto 3002
	if(ptype.eq.0)goto 3628
		
	do 262 i=1,nt1
	do 263 j=1,ntss(2) 
		if(cont(j).eq.0)goto 263
		if(abs(tss1(i)-tss(2,j)).le.inter)cont(j)=0
 263	continue
 262	continue
 
 3628	do i=1,ntss(2)
	do j=i+1,ntss(2)	
		if(scortss2(i).lt.scortss2(j))then 
			sk1=scortss(2,i)
			k1=tss(2,i)
			icontx=cont(i)
			scortss(2,i)=scortss(2,j)
			tss(2,i)=tss(2,j)
			cont(i)=cont(j)
			scortss(2,j)=sk1
			tss(2,j)=k1 
			cont(j)=icontx
		endif
	enddo
	enddo	
	
	do 264 i=1,ntss(2)
		if(cont(i).eq.0)goto 264	
	do 2264 j=i+1,ntss(2)
		if(cont(j).eq.0)goto 2264		
		if(abs(tss(2,i)-tss(2,j)).le.inter)cont(j)=0
 2264	continue 
 264	continue

 3002 	do 265 i=1,ntss(2)
		if(cont(i).eq.0)goto 265 
		nt2=nt2+1 
		tss2(nt2)=tss(2,i)
		scortss2(nt2)=scortss(2,i)
 265	continue
	if(ntss(2).eq.1)goto 2002 
 
 	do i=1,ntss(2)
	do j=i+1,ntss(2)	
		if(tss2(i).lt.tss2(j))then 
			sk1=scortss2(i)
			k1=tss2(i)
			icontx=cont(i)
			scortss2(i)=scortss2(j)
			tss2(i)=tss2(j)
			cont(i)=cont(j)
			scortss2(j)=sk1
			tss2(j)=k1 
			cont(j)=icontx
		endif
	enddo
	enddo	
 	
 2002	nt=nt1+nt2
 	if(nt.eq.0)then 
 		write(3,3500)
		goto 2000	
	endif 
 	if(nt2.eq.0)goto 5002	
	ii=nt1 
	do i=1,nt2
		ii=ii+1
		scortss1(ii)=scortss2(i)
		tss1(ii)=tss2(i)
		tatapos1(ii)=0
	enddo	
	
c ________ Arrangement of TSSs ________
	
 5002 	do i=1,nt
	do j=i+1,nt	
		if(tss1(i).lt.tss1(j))then 
			sk1=scortss1(i)
			sk2=scortata1(i)
			k1=tss1(i)
			k2=tatapos1(i)
			scortss1(i)=scortss1(j)
			scortata1(i)=scortata1(j)
			tss1(i)=tss1(j)
			tatapos1(i)=tatapos1(j)
			scortss1(j)=sk1
			scortata1(j)=sk2 
			tss1(j)=k1 
			tatapos1(j)=k2 
		endif
	enddo
	enddo
	
	do i=1,nt
		cont(i)=1
	enddo	
	if(nt.eq.1)goto 1002
	
  		do 1006 i=1,nt-1

	if(cont(i).eq.0)goto 1006

		do 1005 j=i+1,nt	

	if(cont(j).eq.0)goto 1005
	if(abs(tss1(i)-tss1(j)).le.inter)then
		if(ptype.eq.1)goto 2502
		kx1=abs(tss1(i)-mtss)
		kx2=abs(tss1(j)-mtss)
		if(kx1.le.kx2)then
			cont(j)=0
		else
			cont(i)=0
		endif
		goto 1005
		
 2502	if((tatapos1(i).eq.0).and.(tatapos1(j).ne.0))cont(i)=0 
	if((tatapos1(i).ne.0).and.(tatapos1(j).eq.0))cont(j)=0
	 	
	endif	

 1005		continue  
 1006		continue 
	
c ______ Statistics  

 1002	ntx=0
   	do i=1,nt
 		if(cont(i).eq.1)then
			ntot3=ntot3+1
			ntx=ntx+1
		endif	
	enddo	
	ks=-100000
	ik=0
	do 2015 i=1,nt
		if(cont(i).eq.0)goto 2015
	        kk=tss1(i)-mtss
		if(kk.gt.50)goto 2015
		if(kk.gt.ks)then
			ks=kk
			ik=i
		endif	
 2015	continue 
 	if(ik.eq.0)goto 3468
	nst(ngene)=tss1(ik)-mtss
	if(tatapos1(ik).ne.0)then
		ntot4=ntot4+1
		tata(ngene)=1
	else
		ntot5=ntot5+1
		tata(ngene)=0
	endif	 		 	
	
c _______ Print search results _______

 3468  	write(3,3501)ntx
			do 1007 i=1,nt
		if(cont(i).eq.0)goto 1007
		if(tatapos1(i).ne.0)then 
	write(3,3502)tss1(i),scortss1(i),tatapos1(i),scortata1(i)
	ntot1=ntot1+1
		else
	write(3,3503)tss1(i),scortss1(i) 	
	ntot2=ntot2+1	
		endif	
 1007			continue 

c _______ Go to analyze the next Query _______
 	
 2000	if(iend.eq.0)then
		name=emp
		name=sn
		lseqt=0
		ngene=ngene+1	
		do i=1,2
			ntss(i)=0
		enddo				
		goto 12
	endif	
 88	close(1)
 
c ______ Print Statistics 

	kgene=0
	do 1992 i=1,ngene
		if(nst(i).eq.-100000)goto 1992
		kgene=kgene+1
		nst2(kgene)=nst(i)
		tatax(kgene)=tata(i)
 1992	continue 

	k10=0
	kt10a=0
	kt10b=0
	
	k50=0
	kt50a=0
	kt50b=0
	
	k100=0
	kt100a=0
	kt100b=0
	
	k200=0
	kt200a=0
	kt200b=0
	
	k300=0
	kt300a=0
	kt300b=0
	
	k400=0
	kt400a=0
	kt400b=0
	
	k500=0	
	kt500a=0
	kt500b=0
	
	k600=0 
	kt600a=0
	kt600b=0
	
	k601=0
	kt601a=0
	kt601b=0
 			do 1993 i=1,kgene
			
	if((nst2(i).ge.-10).and.(nst2(i).le.10))then
		k10=k10+1
		if(tatax(i).eq.1)then
			kt10a=kt10a+1
		else
			kt10b=kt10b+1
		endif	
		goto 1993
	endif	
	if(nst2(i).ge.-50)then
		k50=k50+1
		if(tatax(i).eq.1)then
			kt50a=kt50a+1
		else
			kt50b=kt50b+1
		endif	
		goto 1993
	endif	
	if(nst2(i).ge.-100)then
		k100=k100+1
		if(tatax(i).eq.1)then
			kt100a=kt100a+1		
		else
			kt100b=kt100b+1
		endif			
		goto 1993
	endif				
	if(nst2(i).ge.-200)then
		k200=k200+1
		if(tatax(i).eq.1)then
			kt200a=kt200a+1
		else
			kt200b=kt200b+1
		endif			
		goto 1993
	endif				
	if(nst2(i).ge.-300)then
		k300=k300+1
		if(tatax(i).eq.1)then
			kt300a=kt300a+1
		else
			kt300b=kt300b+1
		endif			
		goto 1993
	endif				
	if(nst2(i).ge.-400)then
		k400=k400+1
		if(tatax(i).eq.1)then
			kt400a=kt400a+1
		else
			kt400b=kt400b+1
		endif			
		goto 1993
	endif				
	if(nst2(i).ge.-500)then
		k500=k500+1
		if(tatax(i).eq.1)then
			kt500a=kt500a+1
		else
			kt500b=kt500b+1
		endif			
		goto 1993
	endif				
	if(nst2(i).ge.-600)then
		k600=k600+1
		if(tatax(i).eq.1)then
			kt600a=kt600a+1
		else
			kt600b=kt600b+1
		endif			
		goto 1993
	endif				
	k601=k601+1
	if(tatax(i).eq.1)then
		kt601a=kt601a+1
	else
		kt601b=kt601b+1
	endif	
	
	 				
 1993			continue
 
 
 	if(kgene.eq.0)stop
	
 	p10=k10*100./kgene
 	p50=k50*100./kgene
 	p100=k100*100./kgene
 	p200=k200*100./kgene
 	p300=k300*100./kgene
 	p400=k400*100./kgene
 	p500=k500*100./kgene
 	p600=k600*100./kgene
 	p601=k601*100./kgene


 	if(k10.eq.0)then
		pt10a=0.
		pt10b=0.
		goto 3601
	endif	
	pt10a=kt10a*100./k10 
 	pt10b=kt10b*100./k10 

 3601 	if(k50.eq.0)then
		pt50a=0.
		pt50b=0.
		goto 3602
	endif	 	
        pt50a=kt50a*100./k50 	
 	pt50b=kt50b*100./k50 
	
 3602 	if(k100.eq.0)then
		pt100a=0.
		pt100b=0.
		goto 3603
	endif	
 	pt100a=kt100a*100./k100 
 	pt100b=kt100b*100./k100 
	
 3603 	if(k200.eq.0)then
		pt200a=0.
		pt200b=0.
		goto 3604
	endif	
 	pt200a=kt200a*100./k200 
 	pt200b=kt200b*100./k200 

 3604 	if(k300.eq.0)then
		pt300a=0.
		pt300b=0.
		goto 3605
	endif	
 	pt300a=kt300a*100./k300 
 	pt300b=kt300b*100./k300 

 3605 	if(k400.eq.0)then
		pt400a=0.
		pt400b=0.
		goto 3606
	endif	
 	pt400a=kt400a*100./k400 
 	pt400b=kt400b*100./k400
	
 3606 	if(k500.eq.0)then
		pt500a=0.
		pt500b=0.
		goto 3607
	endif		 
 	pt500a=kt500a*100./k500 
 	pt500b=kt500b*100./k500 
	
 3607 	if(k600.eq.0)then
		pt600a=0.
		pt600b=0.
		goto 3608
	endif	
 	pt600a=kt600a*100./k600 
 	pt600b=kt600b*100./k600 

 3608 	if(k601.eq.0)then
		pt601a=0.
		pt601b=0.
		goto 3609
	endif	
 	pt601a=kt601a*100./k601 
 	pt601b=kt601b*100./k601
	
 3609	ptot4=ntot4*100./kgene
	ptot5=ntot5*100./kgene

c _______________________________
 
	kdns=ltot/ntot3	
 	write(3,1995)
 1995	format(/,50('_'),/,'Summary:',/)
 	write(3,1996)ntot3,kgene,ngene,ntot1,ntot2,
     *	   kdns,ltot
	write(3,1997)kgene,k10,kt10a,kt10b,p10,pt10a,pt10b,
     *	     k50,kt50a,kt50b,p50,pt50a,pt50b,
     *	     k100,kt100a,kt100b,p100,pt100a,pt100b,
     *	     k200,kt200a,kt200b,p200,pt200a,pt200b,
     *	     k300,kt300a,kt300b,p300,pt300a,pt300b,
     *	     k400,kt400a,kt400b,p400,pt400a,pt400b,
     *	     k500,kt500a,kt500b,p500,pt500a,pt500b,
     *	     k600,kt600a,kt600b,p600,pt600a,pt600b,
     *	     k601,kt601a,kt601b,p601,pt601a,pt601b,
     *       ntot4,ptot4,ntot5,ptot5	     
 1996	format(i6,' TSSs in ',i6,' genes out of ',i6,' genes',/,
     *    10x,'TATA+ proms: ',i6,/,10x,'TATA- proms: ',i6,
     *    /,'Prom density:',i6,5x,'in ',i9,' nt') 
 1997	format(//,'Distribution of Closest ',
     *       i6,' TSSs relative to Gene Start '
     *       '(Around or Left):',/, 
     *       '  0-10   Left or Right   ... ',
     *    i6,3x,'[',i6,' TATA+ and ',i6,' TATA-] ... ',
     *    f7.3,'%   [',f7.3,'% .. ',f7.3,'%]',/, 
     *       ' 11-50   Left or Right   ... ',
     *    i6,3x,'[',i6,' TATA+ and ',i6,' TATA-] ... ',
     *    f7.3,'%   [',f7.3,'% .. ',f7.3,'%]',/,      
     *       ' 51-100  Left            ... ',  
     *    i6,3x,'[',i6,' TATA+ and ',i6,' TATA-] ... ',
     *    f7.3,'%   [',f7.3,'% .. ',f7.3,'%]',/,
     *       '101-200  Left            ... ',
     *    i6,3x,'[',i6,' TATA+ and ',i6,' TATA-] ... ',
     *    f7.3,'%   [',f7.3,'% .. ',f7.3,'%]',/,            
     *       '201-300  Left            ... ',
     *    i6,3x,'[',i6,' TATA+ and ',i6,' TATA-] ... ',
     *    f7.3,'%   [',f7.3,'% .. ',f7.3,'%]',/,        
     *       '301-400  Left            ... ',
     *    i6,3x,'[',i6,' TATA+ and ',i6,' TATA-] ... ',
     *    f7.3,'%   [',f7.3,'% .. ',f7.3,'%]',/,       
     *       '401-500  Left            ... ',
     *    i6,3x,'[',i6,' TATA+ and ',i6,' TATA-] ... ',
     *    f7.3,'%   [',f7.3,'% .. ',f7.3,'%]',/,            
     *       '501-600  Left            ... ',
     *    i6,3x,'[',i6,' TATA+ and ',i6,' TATA-] ... ',
     *    f7.3,'%   [',f7.3,'% .. ',f7.3,'%]',/,               
     *       '601/more Left            ... ',
     *    i6,3x,'[',i6,' TATA+ and ',i6,' TATA-] ... ',
     *    f7.3,'%   [',f7.3,'% .. ',f7.3,'%]',//,      
     *       'Totally: ',i6,' [',f7.3,'%] TATA+  and ',
     *    i6,' [',f7.3,'%] TATA- Proms') 
        
	                   
	close(3)
	close(4)

	stop
	end
	
c  _____________________________________________________________________

	function len1(s1)
 
        character*300 s1

        do i=300,1,-1 
                if(s1(i:i).ne.' ')goto 1
        enddo
 1      len1=i

	  return
 	  end

c  _____________________________________________________________________

	function intr(ss,n1,n2)
	character*300 ss,ssw
	integer numb(9)
	do i=1,9
		numb(i)=0
	enddo
	k=0
	do 300 i=n1,n2
		if(ss(i:i).eq.' ')goto 300
		if(ss(i:i).eq.',')goto 300		
		if(ss(i:i).eq.'.')goto 100
		k=k+1
		if(k.gt.9)then 
			intr=0
			return
		endif
		ssw(k:k)=ss(i:i)
 300	continue

 100	do i=1,k
		numb(i)=ichar(ssw(i:i))-48
	enddo
	j=0
	in=0
	do jj=k,1,-1
		j=j+1
		in=in+numb(j)*10**(jj-1)
	enddo
	intr=in
	return
	end

c  _____________________________________________________________________

        subroutine openf(in,xf)
       
c       open file with input name dataf(n)

        character*800 path 
        character*300 fname,xf,emp
        common/envp/path,lpath
	common/emp/emp

	k=len1(xf)
        fname=emp 
        fname(1:lpath)=path(1:lpath)
        fname(lpath+1:lpath+k)=xf(1:k)
        open(in,file=fname,status='old')
        return
        end

c  _____________________________________________________________________

	subroutine stata  

c  Compute TATA feature

         character*1 seq(1000)
	 byte seqb(1000)
	 integer tataposx
	 integer postata(2),ltata(2)	 
	 real frg(4),scorx 
	 real fr1(4,30),frv(4,30),frm(4,30)
	 double precision sv,smax 

	common/seq/seq
	common/seqb/seqb
	common/scorx/scorx
	common/lseq/lseq	
	common/frg/frg
        common/kpos/kpos
	common/postata/postata	
	common/fr1/fr1
	common/ltata/ltata
 
 	do i=1,4
	do j=1,ltata(2) 
		if(fr1(i,j).eq.0.)then
			frv(i,j)=0.00001
		else
			frv(i,j)=fr1(i,j)
		endif	
	enddo
	enddo
		
    	do i=1,4
	do j=1,ltata(2)
		frm(i,j)=log(frv(i,j)/frg(i))								
	enddo
	enddo
	
	smax=-1000000.
	kpos=0
	do i=postata(1),postata(2)-ltata(2)+1 
		sv=0.
		i1=i
		i2=i1+ltata(2)-1
		jj=0
		do j=i1,i2
			kp=seqb(j)
			jj=jj+1	
			sv=sv+frm(kp,jj)
		enddo
		if(smax.lt.sv)then
			smax=sv
			kpos=i
		endif			
	enddo
	scorx=smax

	return
	end	
c  _____________________________________________________________________

	subroutine sinr(k) 

c  Compute INR feature

         character*1 seq(1000)
	 byte seqb(1000)
	 integer posinr(2,2),linr(2,2) 
	 real frg(4),scorx 
	 real fr2(2,4,30),frv(4,30),frm(4,30)
	 double precision sv,smax 

	common/seq/seq
	common/seqb/seqb
	common/scorx/scorx
	common/lseq/lseq	
	common/frg/frg
        common/kpos/kpos
	common/posinr/posinr
	common/fr2/fr2
	common/linr/linr

	do i=1,4
	do j=1,linr(k,2) 
		if(fr2(k,i,j).eq.0.)then
			frv(i,j)=0.00001
		else
			frv(i,j)=fr2(k,i,j)
		endif	
	enddo
	enddo
		
   	do i=1,4
	do j=1,linr(k,2)
		frm(i,j)=log(frv(i,j)/frg(i))	
	enddo
	enddo
	
	smax=-1000000.
	kpos=0
	do i=posinr(k,1),posinr(k,2)-linr(k,2)+1 
		sv=0.
		i1=i
		i2=i1+linr(k,2)-1
		jj=0
		do j=i1,i2
			kp=seqb(j)
			jj=jj+1	
			sv=sv+frm(kp,jj)		
		enddo
		if(smax.lt.sv)then
			smax=sv
			kpos=i
		endif			
	enddo
	scorx=smax

	return
	end	
c  _____________________________________________________________________  

	subroutine syp(k) 

c  Compute YP feature

         character*1 seq(1000)
	 byte seqb(1000)
	 integer posyp(2,2),lyp(2) 
	 real frg(4),scorx 
	 real fr3(2,4,30),frv(4,30),frm(4,30)
	 double precision sv,smax 

	common/seq/seq
	common/seqb/seqb
	common/scorx/scorx
	common/lseq/lseq	
	common/frg/frg
        common/kpos/kpos
	common/posyp/posyp
	common/lyp/lyp 
	common/fr3/fr3
		
	do i=1,4
	do j=1,lyp(k) 
		if(fr3(k,i,j).eq.0.)then
			frv(i,j)=0.00001
		else
			frv(i,j)=fr3(k,i,j)
		endif	
	enddo
	enddo
	
  	do i=1,4
	do j=1,lyp(k)
		frm(i,j)=log(frv(i,j)/frg(i))
	enddo
	enddo
	
	smax=-1000000.
	kpos=0
	do i=posyp(k,1),posyp(k,2)-lyp(k)+1 
		sv=0.
		i1=i
		i2=i1+lyp(k)-1
		jj=0
		do j=i1,i2
			kp=seqb(j)
			jj=jj+1	
			sv=sv+frm(kp,jj)		
		enddo
		if(smax.lt.sv)then
			smax=sv
			kpos=i
		endif			
	enddo
	scorx=smax

	return
	end	
c  _____________________________________________________________________

	subroutine sdp(k) 

c  Compute DPE feature

         character*1 seq(1000)
	 byte seqb(1000)
	 integer posdp(2,2),ldp(2) 
	 real frg(4),scorx 
	 real fr4(2,4,30),frv(4,30),frm(4,30)
	 double precision sv,smax 

	common/seq/seq
	common/seqb/seqb
	common/scorx/scorx
	common/lseq/lseq	
	common/frg/frg
        common/kpos/kpos
	common/posdp/posdp
	common/ldp/ldp 
	common/fr4/fr4	

	do i=1,4
	do j=1,ldp(k) 
		if(fr4(k,i,j).eq.0.)then
			frv(i,j)=0.00001
		else
			frv(i,j)=fr4(k,i,j)
		endif	
	enddo
	enddo
	
  	do i=1,4
	do j=1,ldp(k)
		frm(i,j)=log(frv(i,j)/frg(i))
	enddo
	enddo
	
	smax=-1000000.
	kpos=0
	do i=posdp(k,1),posdp(k,2)-ldp(k)+1 
		sv=0.
		i1=i
		i2=i1+ldp(k)-1
		jj=0
		do j=i1,i2
			kp=seqb(j)
			jj=jj+1	
			sv=sv+frm(kp,jj)		
		enddo
		if(smax.lt.sv)then
			smax=sv
			kpos=i
		endif			
	enddo
	scorx=smax

	return
	end	
c  _____________________________________________________________________

	subroutine strip(k) 
	
c  Compute Triplets feature	
	
         character*1 seq(1000)
	 byte seqb(1000)
	 real scorx,fro3(2,64) 
	 integer pos3(2,2)
	 double precision sv,smax 
	 
	common/seq/seq
	common/seqb/seqb
	common/pos3/pos3 
	common/scorx/scorx
	common/lseq/lseq
	common/fro3/fro3
	  	
	kol1=pos3(k,1) 
	kol3=pos3(k,2)-2 
	lol=64

	sv=0.

 		do i=kol1,kol3 
	i1=seqb(i)
	i2=seqb(i+1)-1
	i3=seqb(i+2)-1
        ind=i1+i2*4+i3*16 
	sv=sv+fro3(k,ind)
        	enddo
c  	scorx=sv
	kz=pos3(k,2)-pos3(k,1)-1
	zk=kz/64.
  	scorx=sv/zk	
			
	return
	end
c  _____________________________________________________________________

	subroutine stetra1(k) 
	
c  Compute Tetra-1 feature	
	
         character*1 seq(1000)
	 byte seqb(1000)
	 real scorx,fro4a(2,256) 
	 integer pos4a(2,2)
	 double precision sv,smax 

	common/seq/seq
	common/seqb/seqb
	common/pos4a/pos4a 
	common/scorx/scorx
	common/lseq/lseq
	common/fro4a/fro4a		

	kol1=pos4a(k,1)	
	kol3=pos4a(k,2)-3
	kol=pos4a(k,2)-pos4a(k,1)-2		
	lol=256
	sv=0.	

		do i=kol1,kol3
	i1=seqb(i)
	i2=seqb(i+1)-1
	i3=seqb(i+2)-1
	i4=seqb(i+3)-1
        ind=i1+i2*4+i3*16+i4*64
	sv=sv+fro4a(k,ind) 
        	enddo
		
c  	scorx=sv

	kz=pos4a(k,2)-pos4a(k,1)-2
	zk=kz/256.
  	scorx=sv/zk	

			
	return
	end
c  _____________________________________________________________________

	subroutine stetra2(k) 
	
c  Compute Tetra-2 feature	
	
         character*1 seq(1000)
	 byte seqb(1000)
	 real scorx,fro4b(2,256) 
	 integer pos4b(2,2)
	 double precision sv,smax 

	common/seq/seq
	common/seqb/seqb
	common/pos4b/pos4b 
	common/scorx/scorx
	common/lseq/lseq
	common/fro4b/fro4b		

	kol1=pos4b(k,1)	
	kol3=pos4b(k,2)-3
	kol=pos4b(k,2)-pos4b(k,1)-2		
	lol=256
	sv=0.	

		do i=kol1,kol3
	i1=seqb(i)
	i2=seqb(i+1)-1
	i3=seqb(i+2)-1
	i4=seqb(i+3)-1
        ind=i1+i2*4+i3*16+i4*64
	sv=sv+fro4b(k,ind)
        	enddo
		
c  	scorx=sv

	kz=pos4b(k,2)-pos4b(k,1)-2
	zk=kz/256.
  	scorx=sv/zk	
			
	return
	end
c  _____________________________________________________________________

	subroutine shexa1(k) 
	
c  Compute Hexa-1 feature	
	
         character*1 seq(1000)
	 byte seqb(1000)
	 real scorx,fro6a(2,4096) 
	 integer pos6a(2,2)
	 double precision sv,smax 

	common/seq/seq
	common/seqb/seqb
	common/pos6a/pos6a
	common/scorx/scorx
	common/lseq/lseq
	common/fro6a/fro6a		

	kol1=pos6a(k,1)
	kol3=pos6a(k,2)-5
	kol=pos6a(k,2)-pos6a(k,1)-4 
	lol=4096
	sv=0.	

		do i=kol1,kol3
	i1=seqb(i)
	i2=seqb(i+1)-1
	i3=seqb(i+2)-1
	i4=seqb(i+3)-1
	i5=seqb(i+4)-1
 	i6=seqb(i+5)-1
        ind=i1+i2*4+i3*16+i4*64+i5*256+i6*1024
	sv=sv+fro6a(k,ind) 
        	enddo
		
c 	scorx=sv

	kz=pos6a(k,2)-pos6a(k,1)-4
	zk=kz/2096.
  	scorx=sv/zk	


			
	return
	end
c  _____________________________________________________________________

	subroutine shexa2(k) 
	
c  Compute Hexa-2 feature	
	
         character*1 seq(1000)
	 byte seqb(1000)
	 real scorx,fro6b(2,4096) 
	 integer pos6b(2,2)
	 double precision sv,smax 

	common/seq/seq
	common/seqb/seqb
	common/pos6b/pos6b
	common/scorx/scorx
	common/lseq/lseq
	common/fro6b/fro6b		

	kol1=pos6b(k,1)
	kol3=pos6b(k,2)-5
	kol=pos6b(k,2)-pos6b(k,1)-4 
	lol=4096
	sv=0.	

		do i=kol1,kol3
	i1=seqb(i)
	i2=seqb(i+1)-1
	i3=seqb(i+2)-1
	i4=seqb(i+3)-1
	i5=seqb(i+4)-1
 	i6=seqb(i+5)-1
        ind=i1+i2*4+i3*16+i4*64+i5*256+i6*1024
	sv=sv+fro6b(k,ind) 
        	enddo
		
c 	scorx=sv

	kz=pos6b(k,2)-pos6b(k,1)-4
	zk=kz/2096.
  	scorx=sv/zk	

			
	return
	end
c  _____________________________________________________________________

	subroutine ssk1(k)

c  Compute CG-skew feature
	
         character*300 emp,s 
         character*1 seq(1000)
	 byte seqb(1000)
	 integer posssk1(2,2)
	 real scorx,frssk1(2,11) 
	 real skx(11,2)
	 double precision skv 

	common/seq/seq
	common/seqb/seqb
        common/skx/skx	
	common/posssk1/posssk1
	common/scorx/scorx
	common/frssk1/frssk1
	common/lseq/lseq	
 		
	c=0.
	g=0.
	a=0.
	t=0.	
	do i=posssk1(k,1),posssk1(k,2)
		if(seqb(i).eq.2)c=c+1.
		if(seqb(i).eq.3)g=g+1.
		if(seqb(i).eq.1)a=a+1
		if(seqb(i).eq.4)t=t+1
	enddo	
	
	if(c.eq.0.)c=0.000001
	if(g.eq.0.)g=0.000001
	if(a.eq.0.)a=0.000001
	if(t.eq.0.)t=0.000001			

	skv=(c-g)/(c+g)
 	
 			do i=1,11
	if((skv.gt.skx(i,1)).and.(skv.le.skx(i,2)))then
		scorx=frssk1(k,i) 
		goto 13
	endif
			enddo
	if(skv.eq.skx(1,1))scorx=frssk1(k,1) 

 13	continue	
		
	return
	end
c  _____________________________________________________________________

	subroutine ssk2(k)

c  Compute AC-skew feature
	
         character*300 emp,s 
         character*1 seq(1000)
	 byte seqb(1000)
	 integer posssk2(2,2)
	 real scorx,frssk2(2,11) 
	 real skx(11,2)
	 double precision skv 
	 
	common/seq/seq
	common/seqb/seqb
        common/skx/skx	
	common/posssk2/posssk2
	common/scorx/scorx
	common/frssk2/frssk2
	common/lseq/lseq	
 		
	c=0.
	g=0.
	a=0.
	t=0.	
	do i=posssk2(k,1),posssk2(k,2)
		if(seqb(i).eq.2)c=c+1.
		if(seqb(i).eq.3)g=g+1.
		if(seqb(i).eq.1)a=a+1
		if(seqb(i).eq.4)t=t+1
	enddo	

	if(c.eq.0.)c=0.000001
	if(g.eq.0.)g=0.000001
	if(a.eq.0.)a=0.000001
	if(t.eq.0.)t=0.000001			

	skv=(a+c-g-t)/(a+c+g+t)
 	
 			do i=1,11
	if((skv.gt.skx(i,1)).and.(skv.le.skx(i,2)))then
		scorx=frssk2(k,i) 
		goto 13
	endif
			enddo
	if(skv.eq.skx(1,1))scorx=frssk2(k,1) 

 13	continue	
		
	return
	end
c  _____________________________________________________________________  

	subroutine sdis1 
	
c Compute TATA-TSS distance feature	
	
         character*300 emp,s
         character*1 seq(1000)
	 byte seqb(1000),icont
	 integer postata(2),ltata(2),d,inx(100,2) 
	 real scorx,frg(4),fr1c(4,30),frv(4,30),frm(4,30)
	 real frdis1(25) 	 		  
	 double precision sv,smax 
	 double precision th1,th2(2),dsr(2,3) 	 

	common/seq/seq
	common/seqb/seqb
	common/scorx/scorx
	common/lseq/lseq	
	common/dsr/dsr		
	common/frg/frg	
	common/th/th1,th2
	common/kint/kint
	common/postata/postata
	common/ltata/ltata 
	common/fr1c/fr1c 
	common/itss/itss
	common/frdis1/frdis1	
		
 	do i=1,kint
		inx(i,1)=(i-1)*2
		inx(i,2)=i*2
	enddo	

	do i=1,4
	do j=1,ltata(1)
		if(fr1c(i,j).eq.0.)then
			frv(i,j)=0.00001
		else
			frv(i,j)=fr1c(i,j)
		endif	
	enddo
	enddo
	
  	do i=1,4
	do j=1,ltata(1) 
		frm(i,j)=log(frv(i,j)/frg(i))
	enddo
	enddo
	
	smax=-1000000.
	do i=postata(1),postata(2)-ltata(1)+1 
		sv=0.
		i1=i
		i2=i1+ltata(1)-1
		jj=0
		do j=i1,i2
			kp=seqb(j)
			jj=jj+1	
			sv=sv+frm(kp,jj)		
		enddo
		if(smax.lt.sv)then
			smax=sv	
			k1=i1
			k2=i2
		endif						
	enddo		
	if(smax.lt.th1)then
 		scorx=dsr(1,1)
		return
	endif 
 	d=abs((2*itss-k1-k2-1)/2)
		do i=1,kint
	if((d.gt.inx(i,1)).and.(d.le.inx(i,2)))then
		scorx=frdis1(i)
		return
	endif	
		enddo
	if(d.eq.inx(1,1))scorx=frdis1(1)

	return
	end
c  _____________________________________________________________________

	subroutine sdis2(k)
	
c Compute TATA-INR distance feature	
	
         character*300 emp,s
         character*1 seq(1000)
	 byte seqb(1000),icont
	 integer postata(2),ltata(2),d,inx(100,2) 
	 real scorx,frg(4),fr1c(4,30),frv1(4,30) 	 		  
 	 integer posinr(2,2),linr(2,2) 
	 real fr2c(2,4,30),frv2(4,30),frm2(4,30)
	 real frdis2(25),frm1(4,30)	 		  
	 double precision sv,smax 
	 double precision th1,th2(2),dsr(2,3) 	 
 
	common/seq/seq
	common/seqb/seqb
	common/scorx/scorx
	common/lseq/lseq	
	common/dsr/dsr		
	common/frg/frg	
	common/th/th1,th2
	common/kint/kint
	common/postata/postata
	common/ltata/ltata 
	common/fr1c/fr1c 
	common/itss/itss
	common/frdis2/frdis2
 	common/posinr/posinr
	common/linr/linr
	common/fr2c/fr2c 
 	
 	do i=1,kint
		inx(i,1)=(i-1)*2
		inx(i,2)=i*2
	enddo	
	
c......... TATA

	do i=1,4
	do j=1,ltata(1)
		if(fr1c(i,j).eq.0.)then
			frv1(i,j)=0.00001
		else
			frv1(i,j)=fr1c(i,j)
		endif	
	enddo
	enddo
	
  	do i=1,4
	do j=1,ltata(1) 
		frm1(i,j)=log(frv1(i,j)/frg(i))
	enddo
	enddo
	
	smax=-1000000.
	do i=postata(1),postata(2)-ltata(1)+1 
		sv=0.
		i1=i
		i2=i1+ltata(1)-1
		jj=0
		do j=i1,i2
			kp=seqb(j)
			jj=jj+1	
			sv=sv+frm1(kp,jj)		
		enddo
		if(smax.lt.sv)then
			smax=sv	
			k1=i1
			k2=i2
		endif						
	enddo	
	
	if(smax.lt.th1)then
		scorx=dsr(1,2) 
		return
	endif 	
	
c......... INR

	do i=1,4
	do j=1,linr(k,1)
		if(fr2c(k,i,j).eq.0.)then
			frv2(i,j)=0.00001
		else
			frv2(i,j)=fr2c(k,i,j)
		endif	
	enddo
	enddo	
	
  	do i=1,4
	do j=1,linr(k,1)
		frm2(i,j)=log(frv2(i,j)/frg(i))	
	enddo
	enddo

	smax=-1000000.
	do i=posinr(k,1),posinr(k,2)-linr(k,1)+1 
		sv=0.
		i1=i
		i2=i1+linr(k,1)-1
		jj=0
		do j=i1,i2
			kp=seqb(j)
			jj=jj+1	
			sv=sv+frm2(kp,jj)		
		enddo			
		if(smax.lt.sv)then
			smax=sv	
			k3=i1
			k4=i2
		endif						
	enddo
			
	if(smax.lt.th2(1))then 
 		scorx=dsr(1,2)
		return
	endif 

 	d=abs((2*k3-k1-k2-1)/2)
		do i=1,kint
	if((d.gt.inx(i,1)).and.(d.le.inx(i,2)))then
		scorx=frdis2(i)
		return
	endif	
		enddo
	if(d.eq.inx(1,1))scorx=frdis2(1)
	
	return
	end	
c  _____________________________________________________________________

	subroutine sdis3(k)
	
c Compute INR-TSS distance feature	
	
         character*300 emp,s
         character*1 seq(1000)
	 byte seqb(1000),icont
	 integer posinr(2,2),linr(2,2),d,inx(100,2) 
	 real scorx,frg(4),fr2c(2,4,30),frv(4,30),frm(4,30)
	 real frdis3(2,25)	 		  
	 double precision sv,smax,th1 
	 double precision th2(2),dsr(2,3) 	 

	common/seq/seq
	common/seqb/seqb
	common/scorx/scorx
	common/lseq/lseq	
	common/dsr/dsr		
	common/frg/frg	
	common/th/th1,th2
	common/kint/kint
	common/posinr/posinr
	common/linr/linr
	common/fr2c/fr2c 
	common/itss/itss
	common/frdis3/frdis3
	
 	do i=1,kint
		inx(i,1)=(i-1)*2
		inx(i,2)=i*2
	enddo	
 
	do i=1,4
	do j=1,linr(k,1)
		if(fr2c(k,i,j).eq.0.)then
			frv(i,j)=0.00001
		else
			frv(i,j)=fr2c(k,i,j)
		endif	
	enddo
	enddo
	
  	do i=1,4
	do j=1,linr(k,1)
		frm(i,j)=log(frv(i,j)/frg(i))	
	enddo
	enddo
	
	smax=-1000000.
	do i=posinr(k,1),posinr(k,2)-linr(k,1)+1 
		sv=0.
		i1=i
		i2=i1+linr(k,1)-1
		jj=0
		do j=i1,i2
			kp=seqb(j)
			jj=jj+1	
			sv=sv+frm(kp,jj)		
		enddo
		if(smax.lt.sv)then
			smax=sv	
			k1=i1
			k2=i2
		endif						
	enddo	
 	
	if(smax.lt.th2(k))then 
 		scorx=dsr(k,3)
		return
	endif 
 	d=abs((2*itss-k1-k2-1)/2)
	
		do i=1,kint
	if((d.gt.inx(i,1)).and.(d.le.inx(i,2)))then
		scorx=frdis3(k,i)
		return
	endif	
		enddo
	if(d.eq.inx(1,1))scorx=frdis3(k,1)

	return
	end
c  _____________________________________________________________________

	subroutine sred1(k)
	
c Compute RE density 1 (+strand) feature	
	
         character*300 emp,s
         character*1 seq(1000),re1seq(2,10000,80)
	 integer re1l(2,10000),re1m(2,10000)
	 integer posre1(2,2),nsite1(2) 
	 real frre1(2,10000),scorx  

	common/seq/seq	
	common/scorx/scorx
	common/lseq/lseq	
        common/re1l/re1l
	common/re1m/re1m	 
 	common/posre1/posre1
	common/re1seq/re1seq 
	common/nsite1/nsite1
 	common/frre1/frre1
	 
	scorx=0.
	do 500 ir=1,nsite1(k)
	do 400 i=posre1(k,1),posre1(k,2)-re1l(k,ir)+1 
		i1=i
		i2=i1+re1l(k,ir)-1
		j=0
		ms=0
			
	do  300 ii=i1,i2
		j=j+1
		if(seq(ii).eq.re1seq(k,ir,j))goto 300
			ms=ms+1
 		if(ms.gt.re1m(k,ir))goto 400	
 300	continue 
  	scorx=scorx+frre1(k,ir) 
	goto 500
 400		continue
 500			continue

	return
	end
c  _____________________________________________________________________	

	subroutine sred2(k)
	
c Compute RE density 1 (+/- strands) feature	
	
         character*300 emp,s
         character*1 seq(1000),re2seq(2,10000,80)
         character*1 seqc(1000) 
	 integer re2l(2,10000),re2m(2,10000)
	 integer posre2(2,2),nsite2(2) 
	 real frre2(2,10000),scorx  

	common/seq/seq
	common/seqc/seqc	
	common/scorx/scorx
	common/lseq/lseq	
        common/re2l/re2l
	common/re2m/re2m	 
 	common/posre2/posre2
	common/re2seq/re2seq 
	common/nsite2/nsite2
 	common/frre2/frre2	
 	 
	kchain=2

	scorx=0.
	do 500 ir=1,nsite2(k)
	do 400 i=posre2(k,1),posre2(k,2)-re2l(k,ir)+1 
	i1=i
	i2=i1+re2l(k,ir)-1
	j=0
	ms=0
	do  300 ii=i1,i2
		j=j+1
		if(seq(ii).eq.re2seq(k,ir,j))goto 300
		ms=ms+1
		if(ms.gt.re2m(k,ir))goto 400
 300	continue 
  	scorx=scorx+frre2(k,ir) 
	goto 500
	
 400		continue
 
 	ibin1=lseq-posre2(k,2)+1 
	ibin2=lseq-posre2(k,1)+1 	
 
		do 1401 i=ibin1,ibin2-re2l(k,ir)+1 
	i1=i
	i2=i1+re2l(k,ir)-1
	j=0
	ms=0
	do  1301 ii=i1,i2
		j=j+1
		if(seqc(ii).eq.re2seq(k,ir,j))goto 1301
		ms=ms+1
		if(ms.gt.re2m(k,ir))goto 1401
 1301	continue 
 	scorx=scorx+frre2(k,ir)
	goto 500
	
 1401		continue
	
 500			continue

	return
	end
	
c  _____________________________________________________________________	

	 subroutine visanx(k,nf)
	 
	 real scor(16),pi,c(2)	  
	 integer nline(2),h,hh
	 double precision w(2,2,30,30),bias(2,2,30)
	 double precision neuron_layer(3,30),neuron_input(2,30)
	 double precision dotprod,scortssx,pix
	 double precision min(2,30),max(2,30)	
 
	 common/scor/scor
	 common/ipred/ipred
	 common/scortssx/scortssx
	 common/pi/pi
	 common/c/c 
	 common/nline/nline 
	 common/minmax/min,max
	 common/w/w
	 common/bias/bias	
	 
 1	 format(a) 	 

	pix=2./pi
	ipred=1

c ___ Compute Visan scores and estimate: "ktss" is TSS or not ___

		do i=1,nf
	neuron_layer(1,i)=2.*((scor(i)-min(k,i))/
     *  	(max(k,i)-min(k,i))-0.5)
		enddo
c _____
	 h=1
	 hh=2
	 		do j=1,nline(k) 
	 dotprod=0.
	 	do i=1,nf                             
	 dotprod=dotprod+neuron_layer(h,i)*w(k,h,j,i) 
	 	enddo
	 
	 neuron_input(h,j)=bias(k,h,j)+dotprod 
	 neuron_layer(hh,j)=pix*datan(neuron_input(h,j))	 
			 enddo

c _____			 
	 h=2
	 hh=3
	 		do j=1,2                           !  (for both TATA+ and TATA-)
	 dotprod=0.
	 	do i=1,nline(k)                            ! =16 (TATA+) or =13 (TATA-) 
	 dotprod=dotprod+neuron_layer(h,i)*w(k,h,j,i) 
	 	enddo 	
	 neuron_input(h,j)=bias(k,h,j)+dotprod 
	 neuron_layer(hh,j)=pix*datan(neuron_input(h,j)) 
			 enddo
c _____	 
 	 
	 scortssx=neuron_layer(3,2)-neuron_layer(3,1)
	 ipred=1 	  
	 if(scortssx.lt.c(k))ipred=0
	 
	 return
	 end

c ___________________________________________________________________

	subroutine help

	print*,' '
	print*,' '
	print*,' '			
        print*,'TSSPlant  -i:<parI> [-o:<parO>] ',
     *	       '[-p:<parP>] [-s:<parS>] [-t:<parT>]'  
	print*,' '
	print*,' '	
	print*,'Options/Arguments:'	
	print*,' '
 
    	print*,'<parI>  Input File  with Query DNA ',
     *	       'sequence(s) in FASTA format (max. length: '
        print*,'        100 000 nt). If this option/argument is ',
     *	       'absent or number of arguments is more than '
        print*,'        6 or no any argument is given, Program ',
     *         'display HELP information and ends.' 
        print*,'        In Query sequences symbols besides ',
     * 	       'of "a/A", "g/G", "c/C" and "t/T" are accepted as "n".'  
	print*,' '	    
	print*,'<parL>  Left and Right boundaries of an ',
     *	      'interval applied for selection' 
        print*,'        of a single TSS (for more information ',    
     *         'see below: <parT>).' 
	print*,'        Lminimum = 100' 
        print*,'        Default: 300'
	print*,' '
    	print*,'<parO>  Output File'  
 	print*,'        Default: TSSPlant.res' 	       
	print*,' '  
  	print*,'<parP>  Print (y) or not (n) Query sequence' 
	print*,'        Default: n'
	print*,' '    
	print*,'<parS>  Position of the annotated gene start' 
	print*,'        Default: Length of Query sequence (as ',
     *	       'annotated gene start)'
	print*,' '   
	print*,'<parT>  After identification of all TSSs ',
     *	       'for TATA and TATA-less promoters'
	print*,'        with score higher the threshold ', 
     *	       'determined during Learning procedure,' 
	print*,'        for every interval (+L <= TSS <= -L),',
     *	       ' 2 TSSs (for TATA and TATA-less' 
	print*,'        promoters) with highest score are selected.' 
	print*,'        Parameter "L" is given by <parL> option.'
	print*,'        If parT="y", only TSS for TATA promoter ',
     *	       'is selected as predicted TSS;'  
	print*,'        parT="n", only TSS closer to the annotated, ',
     *	       'gene start is selected as predicted TSS.'  
	print*,'        Default: "n"'   
	print*,' '   
	
	return
	end
	
	
