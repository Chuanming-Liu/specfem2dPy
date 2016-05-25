#!/usr/bin/env python
import obspy
import pyaftan as ftan
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.pylab as plb
import os
from functools import partial
import multiprocessing 
import pyasdf
import stations

class ftanParam(object):
    """ An object to handle ftan output parameters
    -----------------------------------------------------------------------------------------------------
    Basic FTAN parameters:
    nfout1_1 - output number of frequencies for arr1, (integer*4)
    arr1_1   - preliminary results.
              Description: real*8 arr1(8,n), n >= nfin)
              arr1(1,:) -  central periods, s
              arr1(2,:) -  observed periods, s
              arr1(3,:) -  group velocities, km/s
              arr1(4,:) -  phase velocities, km/s or phase if nphpr=0, rad
              arr1(5,:) -  amplitudes, Db
              arr1(6,:) -  discrimination function
              arr1(7,:) -  signal/noise ratio, Db
              arr1(8,:) -  maximum half width, s
              arr1(9,:) -  amplitudes
    arr2_1   - final results
    nfout2_1 - output number of frequencies for arr2, (integer*4)
              Description: real*8 arr2(7,n), n >= nfin)
              If nfout2 == 0, no final result.
              arr2(1,:) -  central periods, s
              arr2(2,:) -  observed periods, s
              arr2(3,:) -  group velocities, km/sor phase if nphpr=0, rad
              arr2(4,:) -  phase velocities, km/s
              arr2(5,:) -  amplitudes, Db
              arr2(6,:) -  signal/noise ratio, Db
              arr2(7,:) -  maximum half width, s
              arr2(8,:) -  amplitudes
    tamp_1      -  time to the beginning of ampo table, s (real*8)
    nrow_1      -  number of rows in array ampo, (integer*4)
    ncol_1      -  number of columns in array ampo, (integer*4)
    amp_1       -  Ftan amplitude array, Db, (real*8)
    ierr_1   - completion status, =0 - O.K.,           (integer*4)
                                 =1 - some problems occures
                                 =2 - no final results
    ==========================================================
    Phase-Matched-Filtered FTAN parameters:
    nfout1_2 - output number of frequencies for arr1, (integer*4)
    arr1_2   - preliminary results.
             Description: real*8 arr1(8,n), n >= nfin)
             arr1(1,:) -  central periods, s (real*8)
             arr1(2,:) -  apparent periods, s (real*8)
             arr1(3,:) -  group velocities, km/s (real*8)
             arr1(4,:) -  phase velocities, km/s (real*8)
             arr1(5,:) -  amplitudes, Db (real*8)
             arr1(6,:) -  discrimination function, (real*8)
             arr1(7,:) -  signal/noise ratio, Db (real*8)
             arr1(8,:) -  maximum half width, s (real*8)
             arr1(9,:) -  amplitudes 
    arr2_2   - final results
    nfout2_2 - output number of frequencies for arr2, (integer*4)
             Description: real*8 arr2(7,n), n >= nfin)
             If nfout2 == 0, no final results.
             arr2(1,:) -  central periods, s (real*8)
             arr2(2,:) -  apparent periods, s (real*8)
             arr2(3,:) -  group velocities, km/s (real*8)
             arr1(4,:) -  phase velocities, km/s (real*8)
             arr2(5,:) -  amplitudes, Db (real*8)
             arr2(6,:) -  signal/noise ratio, Db (real*8)
             arr2(7,:) -  maximum half width, s (real*8)
             arr2(8,:) -  amplitudes
    tamp_2      -  time to the beginning of ampo table, s (real*8)
    nrow_2      -  number of rows in array ampo, (integer*4)
    ncol_2      -  number of columns in array ampo, (integer*4)
    amp_2       -  Ftan amplitude array, Db, (real*8)
    ierr_2   - completion status, =0 - O.K.,           (integer*4)
                                =1 - some problems occures
                                =2 - no final results
    -----------------------------------------------------------------------------------------------------
    """
    def __init__(self):
        # Parameters for first iteration
        self.nfout1_1=0
        self.arr1_1=np.array([])
        self.nfout2_1=0
        self.arr2_1=np.array([])
        self.tamp_1=0.
        self.nrow_1=0
        self.ncol_1=0
        self.ampo_1=np.array([],dtype='float32')
        self.ierr_1=0
        # Parameters for second iteration
        self.nfout1_2=0
        self.arr1_2=np.array([])
        self.nfout2_2=0
        self.arr2_2=np.array([])
        self.tamp_2=0.
        self.nrow_2=0
        self.ncol_2=0
        self.ampo_2=np.array([])
        self.ierr_2=0
        # Flag for existence of predicted phase dispersion curve
        self.preflag=False

    def writeDISP(self, fnamePR):
        """
        Write FTAN parameters to DISP files given a prefix.
        fnamePR: file name prefix
        _1_DISP.0: arr1_1
        _1_DISP.1: arr2_1
        _2_DISP.0: arr1_2
        _2_DISP.1: arr2_2
        """
        if self.nfout1_1!=0:
            f10=fnamePR+'_1_DISP.0'
            f=open(f10,'w')
            for i in np.arange(self.nfout1_1):
                tempstr='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf  \n' %( i, self.arr1_1[0,i] , self.arr1_1[1,i] , self.arr1_1[2,i] , self.arr1_1[3,i]  \
                    , self.arr1_1[4,i] , self.arr1_1[5,i] , self.arr1_1[6,i] )
                f.writelines(tempstr)
            f.close()
        if self.nfout2_1!=0:
            f11=fnamePR+'_1_DISP.1'
            f=open(f11,'w')
            for i in np.arange(self.nfout2_1):
                tempstr='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %8.3lf  \n' %( i, self.arr2_1[0,i], self.arr2_1[1,i] , self.arr2_1[2,i] , self.arr2_1[3,i]  \
                    , self.arr2_1[4,i] , self.arr2_1[5,i]  )
                f.writelines(tempstr)
            f.close()
        if self.nfout1_2!=0:
            f20=fnamePR+'_2_DISP.0'
            f=open(f20,'w')
            for i in np.arange(self.nfout1_2):
                tempstr='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %12.4lf %8.3lf \n' %( i, self.arr1_2[0,i], self.arr1_2[1,i] , self.arr1_2[2,i] , self.arr1_2[3,i]  \
                    , self.arr1_2[4,i] , self.arr1_2[5,i] , self.arr1_2[6,i] )
                f.writelines(tempstr)
            f.close()
        if self.nfout2_2!=0:
            f21=fnamePR+'_2_DISP.1'
            f=open(f21,'w')
            for i in np.arange(self.nfout2_2):
                tempstr='%4d %10.4lf %10.4lf %12.4lf %12.4lf %12.4lf %8.3lf  \n' %( i, self.arr2_2[0,i], self.arr2_2[1,i] , self.arr2_2[2,i] , self.arr2_2[3,i]  \
                    , self.arr2_2[4,i] , self.arr2_2[5,i]  )
                f.writelines(tempstr)
            f.close()
        return

    def FTANcomp(self, inftanparam, compflag=4):
        """
        Compare aftan results for two ftanParam objects.
        """
        fparam1=self
        fparam2=inftanparam
        if compflag==1:
            obper1=fparam1.arr1_1[1,:fparam1.nfout1_1]
            gvel1=fparam1.arr1_1[2,:fparam1.nfout1_1]
            phvel1=fparam1.arr1_1[3,:fparam1.nfout1_1]
            obper2=fparam2.arr1_1[1,:fparam2.nfout1_1]
            gvel2=fparam2.arr1_1[2,:fparam2.nfout1_1]
            phvel2=fparam2.arr1_1[3,:fparam2.nfout1_1]
        elif compflag==2:
            obper1=fparam1.arr2_1[1,:fparam1.nfout2_1]
            gvel1=fparam1.arr2_1[2,:fparam1.nfout2_1]
            phvel1=fparam1.arr2_1[3,:fparam1.nfout2_1]
            obper2=fparam2.arr2_1[1,:fparam2.nfout2_1]
            gvel2=fparam2.arr2_1[2,:fparam2.nfout2_1]
            phvel2=fparam2.arr2_1[3,:fparam2.nfout2_1]
        elif compflag==3:
            obper1=fparam1.arr1_2[1,:fparam1.nfout1_2]
            gvel1=fparam1.arr1_2[2,:fparam1.nfout1_2]
            phvel1=fparam1.arr1_2[3,:fparam1.nfout1_2]
            obper2=fparam2.arr1_2[1,:fparam2.nfout1_2]
            gvel2=fparam2.arr1_2[2,:fparam2.nfout1_2]
            phvel2=fparam2.arr1_2[3,:fparam2.nfout1_2]
        else:
            obper1=fparam1.arr2_2[1,:fparam1.nfout2_2]
            gvel1=fparam1.arr2_2[2,:fparam1.nfout2_2]
            phvel1=fparam1.arr2_2[3,:fparam1.nfout2_2]
            obper2=fparam2.arr2_2[1,:fparam2.nfout2_2]
            gvel2=fparam2.arr2_2[2,:fparam2.nfout2_2]
            phvel2=fparam2.arr2_2[3,:fparam2.nfout2_2]
        plb.figure()
        ax = plt.subplot()
        ax.plot(obper1, gvel1, '--k', lw=3) #
        ax.plot(obper2, gvel2, '-.b', lw=3)
        plt.xlabel('Period(s)')
        plt.ylabel('Velocity(km/s)')
        plt.title('Group Velocity Comparison')
        if (fparam1.preflag==True and fparam2.preflag==True):
            plb.figure()
            ax = plt.subplot()
            ax.plot(obper1, phvel1, '--k', lw=3) #
            ax.plot(obper2, phvel2, '-.b', lw=3)
            plt.xlabel('Period(s)')
            plt.ylabel('Velocity(km/s)')
            plt.title('Phase Velocity Comparison')
        return

class snrParam(object):
    """
    SNR parameters:
        suffix: p=positve lag  n=negative lag s=symmetric lag
        amp: largest amplitude measurement for each period
        snr: SNR measurement for each period
        nrms: noise rms measurement for each period
        oper: observed period
    """
    def __init__(self):
        self.amp_p=np.array([])
        self.snr_p=np.array([])
        self.nrms_p=np.array([])
        self.oper_p=np.array([])
        self.amp_n=np.array([])
        self.snr_n=np.array([])
        self.nrms_n=np.array([])
        self.oper_n=np.array([])
        self.amp_s=np.array([])
        self.snr_s=np.array([])
        self.nrms_s=np.array([])
        self.oper_s=np.array([])
        return

    def writeAMPSNR(self, fnamePR):
        """
        writeAMPSNR:
        Write output SNR parameters to text files
        _pos_amp_snr - positive lag
        _neg_amp_snr - negative lag
        _amp_snr     - symmetric lag
        """
        len_p=len(self.amp_p)
        len_n=len(self.amp_n)
        len_s=len(self.amp_s)
        if len_p!=0:
            fpos=fnamePR+'_pos_amp_snr'
            f=open(fpos,'w')
            for i in np.arange(len_p):
                tempstr='%8.4f   %.5g  %8.4f  \n' %(  self.oper_p[i] , self.amp_p[i],  self.snr_p[i] )
                f.writelines(tempstr)
            f.close()
        if len_n!=0:
            fneg=fnamePR+'_neg_amp_snr'
            f=open(fneg,'w')
            for i in np.arange(len_n):
                tempstr='%8.4f   %.5g  %8.4f  \n' %(   self.oper_n[i] , self.amp_n[i],  self.snr_n[i] )
                f.writelines(tempstr)
            f.close()
        if len_s!=0:
            fsym=fnamePR+'_amp_snr'
            f=open(fsym,'w')
            for i in np.arange(len_s):
                tempstr='%8.4f   %.5g  %8.4f  \n' %(   self.oper_s[i] , self.amp_s[i],  self.snr_s[i] )
                f.writelines(tempstr)
            f.close()
        return

class specfem2dtrace(obspy.core.trace.Trace):
    """
    specfem2dtrace:
    A derived class inherited from obspy.core.trace.Trace. This derived class have a variety of new member functions
    """
    def init_ftanParam(self):
        """
        Initialize ftan parameters
        """
        self.ftanparam=ftanParam()
        
    def init_snrParam(self):
        """
        Initialize SNR parameters
        """
        self.SNRParam=snrParam()
        
    def reverse(self):
        """
        Reverse the trace
        """
        self.data=self.data[::-1]
        return
    
    def aftan(self, pmf=True, piover4=-1.0, vmin=1.5, vmax=5.0, tmin=4.0, \
        tmax=30.0, tresh=20.0, ffact=1.0, taperl=1.0, snr=0.2, fmatch=1.0, phvelname='', predV=np.array([]) ):

        """ (Automatic Frequency-Time ANalysis) aftan analysis:
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        pmf         - flag for Phase-Matched-Filtered output (default: True)
        piover4     - phase shift = pi/4*piover4, for cross-correlation piover4 should be -1.0
        vmin        - minimal group velocity, km/s
        vmax        - maximal group velocity, km/s
        tmin        - minimal period, s
        tmax        - maximal period, s
        tresh       - treshold for jump detection, usualy = 10, need modifications
        ffact       - factor to automatic filter parameter, usualy =1
        taperl      - factor for the left end seismogram tapering, taper = taperl*tmax,    (real*8)
        snr         - phase match filter parameter, spectra ratio to determine cutting point for phase matched filter
        fmatch      - factor to length of phase matching window
        fname       - SAC file name
        phvelname   - predicted phase velocity file name
        predV          - predicted phase velocity curve
        
        Output:
        self.ftanparam, a object of ftanParam class, to store output aftan results
        -----------------------------------------------------------------------------------------------------
        References:
        Levshin, A. L., and M. H. Ritzwoller. Automated detection, extraction, and measurement of regional surface waves.
             Monitoring the Comprehensive Nuclear-Test-Ban Treaty: Surface Waves. Birkh?user Basel, 2001. 1531-1545.
        Bensen, G. D., et al. Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements.
             Geophysical Journal International 169.3 (2007): 1239-1260.
        """
        try:
            self.ftanparam
        except:
            self.init_ftanParam()
        try:
            dist=self.stats.sac.dist
        except:
            dist = np.sqrt( (self.stats.sac.evlo - self.stats.sac.stlo)**2 + (self.stats.sac.evla - self.stats.sac.stla)**2 )
            self.stats.sac.dist=dist
        if (phvelname==''):
            phvelname='./ak135.disp'
        nprpv = 0
        phprper=np.zeros(300)
        phprvel=np.zeros(300)
        if predV.size != 0:
            phprper=predV[:,0]
            phprvel=predV[:,1]
            nprpv = predV[:,0].size
            phprper=np.append( phprper, np.zeros(300-phprper.size) )
            phprvel=np.append( phprvel, np.zeros(300-phprvel.size) )
            self.ftanparam.preflag=True
        elif os.path.isfile(phvelname):
            php=np.loadtxt(phvelname)
            phprper=php[:,0]
            phprvel=php[:,1]
            nprpv = php[:,0].size
            phprper=np.append( phprper, np.zeros(300-phprper.size) )
            phprvel=np.append( phprvel, np.zeros(300-phprvel.size) )
            self.ftanparam.preflag=True
        nfin = 64
        npoints = 5  #  only 3 points in jump
        perc    = 50.0 # 50 % for output segment
        tempsac=self.copy()
        tb=self.stats.sac.b
        length=len(tempsac.data)
        if length>32768:
            print "Warning: length of seismogram is larger than 32768!"
            nsam=32768
            tempsac.data=tempsac.data[:nsam]
            tempsac.stats.e=(nsam-1)*tempsac.stats.delta+tb
            sig=tempsac.data
        else:
            sig=np.append(tempsac.data, np.zeros( float(32768-tempsac.data.size) ) )
            nsam=int( float (tempsac.stats.npts) )### for unknown reasons, this has to be done, nsam=int(tempsac.stats.npts)  won't work as an input for aftan
        dt=tempsac.stats.delta
        # Start to do aftan utilizing pyaftan
        self.ftanparam.nfout1_1,self.ftanparam.arr1_1,self.ftanparam.nfout2_1,self.ftanparam.arr2_1,self.ftanparam.tamp_1, \
                self.ftanparam.nrow_1,self.ftanparam.ncol_1,self.ftanparam.ampo_1, self.ftanparam.ierr_1= ftan.aftanpg(piover4, nsam, \
                    sig, tb, dt, dist, vmin, vmax, tmin, tmax, tresh, ffact, perc, npoints, taperl, nfin, snr, nprpv, phprper, phprvel)
        if pmf==True:
            if self.ftanparam.nfout2_1<3:
                return
            npred = self.ftanparam.nfout2_1
            tmin2 = self.ftanparam.arr2_1[1,0]
            tmax2 = self.ftanparam.arr2_1[1,self.ftanparam.nfout2_1-1]
            pred=np.zeros((2,300))
            pred[:,0:100]=self.ftanparam.arr2_1[1:3,:]
            pred=pred.T
            self.ftanparam.nfout1_2,self.ftanparam.arr1_2,self.ftanparam.nfout2_2,self.ftanparam.arr2_2,self.ftanparam.tamp_2, \
                    self.ftanparam.nrow_2,self.ftanparam.ncol_2,self.ftanparam.ampo_2, self.ftanparam.ierr_2 = ftan.aftanipg(piover4,nsam, \
                        sig,tb,dt,dist,vmin,vmax,tmin2,tmax2,tresh,ffact,perc,npoints,taperl,nfin,snr,fmatch,npred,pred,nprpv,phprper,phprvel)
        self.ftanparam.station_id=self.stats.network+'.'+self.stats.station 
        return

    def plotftan(self, plotflag=3, sacname=''):
        """
        Plot ftan diagram:
        This function plot ftan diagram.
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        plotflag -
            0: only Basic FTAN
            1: only Phase Matched Filtered FTAN
            2: both
            3: both in one figure
        sacname - sac file name than can be used as the title of the figure
        -----------------------------------------------------------------------------------------------------
        """
        try:
            fparam=self.ftanparam
            if fparam.nfout1_1==0:
                return "Error: No Basic FTAN parameters!"
            dt=self.stats.delta
            dist=self.stats.sac.dist
            if (plotflag!=1 and plotflag!=3):
                v1=dist/(fparam.tamp_1+np.arange(fparam.ncol_1)*dt)
                ampo_1=fparam.ampo_1[:fparam.ncol_1,:fparam.nrow_1]
                obper1_1=fparam.arr1_1[1,:fparam.nfout1_1]
                gvel1_1=fparam.arr1_1[2,:fparam.nfout1_1]
                phvel1_1=fparam.arr1_1[3,:fparam.nfout1_1]
                plb.figure()
                ax = plt.subplot()
                p=plt.pcolormesh(obper1_1, v1, ampo_1, cmap='gist_rainbow',shading='gouraud')
                ax.plot(obper1_1, gvel1_1, '--k', lw=3) #
                if (fparam.preflag==True):
                    ax.plot(obper1_1, phvel1_1, '--w', lw=3) #

                if (fparam.nfout2_1!=0):
                    obper2_1=fparam.arr2_1[1,:fparam.nfout2_1]
                    gvel2_1=fparam.arr2_1[2,:fparam.nfout2_1]
                    phvel2_1=fparam.arr2_1[3,:fparam.nfout2_1]
                    ax.plot(obper2_1, gvel2_1, '-k', lw=3) #
                    if (fparam.preflag==True):
                        ax.plot(obper2_1, phvel2_1, '-w', lw=3) #
                cb = plt.colorbar(p, ax=ax)
                Tmin1=obper1_1[0]
                Tmax1=obper1_1[fparam.nfout1_1-1]
                vmin1= v1[fparam.ncol_1-1]
                vmax1=v1[0]
                plt.axis([Tmin1, Tmax1, vmin1, vmax1])
                plt.xlabel('Period(s)')
                plt.ylabel('Velocity(km/s)')
                plt.title('Basic FTAN Diagram '+sacname,fontsize=15)

            if fparam.nfout1_2==0 and plotflag!=0:
                return "Error: No PMF FTAN parameters!"
            if (plotflag!=0 and plotflag!=3):
                v2=dist/(fparam.tamp_2+np.arange(fparam.ncol_2)*dt)
                ampo_2=fparam.ampo_2[:fparam.ncol_2,:fparam.nrow_2]
                obper1_2=fparam.arr1_2[1,:fparam.nfout1_2]
                gvel1_2=fparam.arr1_2[2,:fparam.nfout1_2]
                phvel1_2=fparam.arr1_2[3,:fparam.nfout1_2]
                plb.figure()
                ax = plt.subplot()
                p=plt.pcolormesh(obper1_2, v2, ampo_2, cmap='gist_rainbow',shading='gouraud')
                ax.plot(obper1_2, gvel1_2, '--k', lw=3) #
                if (fparam.preflag==True):
                    ax.plot(obper1_2, phvel1_2, '--w', lw=3) #

                if (fparam.nfout2_2!=0):
                    obper2_2=fparam.arr2_2[1,:fparam.nfout2_2]
                    gvel2_2=fparam.arr2_2[2,:fparam.nfout2_2]
                    phvel2_2=fparam.arr2_2[3,:fparam.nfout2_2]
                    ax.plot(obper2_2, gvel2_2, '-k', lw=3) #
                    if (fparam.preflag==True):
                        ax.plot(obper2_2, phvel2_2, '-w', lw=3) #
                cb = plt.colorbar(p, ax=ax)
                Tmin2=obper1_2[0]
                Tmax2=obper1_2[fparam.nfout1_2-1]
                vmin2= v2[fparam.ncol_2-1]
                vmax2=v2[0]
                plt.axis([Tmin2, Tmax2, vmin2, vmax2])
                plt.xlabel('Period(s)')
                plt.ylabel('Velocity(km/s)')
                plt.title('PMF FTAN Diagram '+sacname,fontsize=15)

            if ( plotflag==3 ):
                v1=dist/(fparam.tamp_1+np.arange(fparam.ncol_1)*dt)
                ampo_1=fparam.ampo_1[:fparam.ncol_1,:fparam.nrow_1]
                obper1_1=fparam.arr1_1[1,:fparam.nfout1_1]
                gvel1_1=fparam.arr1_1[2,:fparam.nfout1_1]
                phvel1_1=fparam.arr1_1[3,:fparam.nfout1_1]
                plb.figure(num=None, figsize=(18, 16), dpi=80, facecolor='w', edgecolor='k')
                ax = plt.subplot(2,1,1)
                p=plt.pcolormesh(obper1_1, v1, ampo_1, cmap='gist_rainbow',shading='gouraud')
                ax.plot(obper1_1, gvel1_1, '--k', lw=3) #
                if (fparam.preflag==True):
                    ax.plot(obper1_1, phvel1_1, '--w', lw=3) #
                if (fparam.nfout2_1!=0):
                    obper2_1=fparam.arr2_1[1,:fparam.nfout2_1]
                    gvel2_1=fparam.arr2_1[2,:fparam.nfout2_1]
                    phvel2_1=fparam.arr2_1[3,:fparam.nfout2_1]
                    ax.plot(obper2_1, gvel2_1, '-k', lw=3) #
                    if (fparam.preflag==True):
                        ax.plot(obper2_1, phvel2_1, '-w', lw=3) #
                cb = plt.colorbar(p, ax=ax)
                Tmin1=obper1_1[0]
                Tmax1=obper1_1[fparam.nfout1_1-1]
                vmin1= v1[fparam.ncol_1-1]
                vmax1=v1[0]
                plt.axis([Tmin1, Tmax1, vmin1, vmax1])
                plt.xlabel('Period(s)')
                plt.ylabel('Velocity(km/s)')
                plt.title('Basic FTAN Diagram '+sacname)

                v2=dist/(fparam.tamp_2+np.arange(fparam.ncol_2)*dt)
                ampo_2=fparam.ampo_2[:fparam.ncol_2,:fparam.nrow_2]
                obper1_2=fparam.arr1_2[1,:fparam.nfout1_2]
                gvel1_2=fparam.arr1_2[2,:fparam.nfout1_2]
                phvel1_2=fparam.arr1_2[3,:fparam.nfout1_2]

                ax = plt.subplot(2,1,2)
                p=plt.pcolormesh(obper1_2, v2, ampo_2, cmap='gist_rainbow',shading='gouraud')
                ax.plot(obper1_2, gvel1_2, '--k', lw=3) #
                if (fparam.preflag==True):
                    ax.plot(obper1_2, phvel1_2, '--w', lw=3) #

                if (fparam.nfout2_2!=0):
                    obper2_2=fparam.arr2_2[1,:fparam.nfout2_2]
                    gvel2_2=fparam.arr2_2[2,:fparam.nfout2_2]
                    phvel2_2=fparam.arr2_2[3,:fparam.nfout2_2]
                    ax.plot(obper2_2, gvel2_2, '-k', lw=3) #
                    if (fparam.preflag==True):
                        ax.plot(obper2_2, phvel2_2, '-w', lw=3) #
                cb = plt.colorbar(p, ax=ax)
                Tmin2=obper1_2[0]
                Tmax2=obper1_2[fparam.nfout1_2-1]
                vmin2= v2[fparam.ncol_2-1]
                vmax2=v2[0]
                plt.axis([Tmin2, Tmax2, vmin2, vmax2])
                plt.xlabel('Period(s)')
                plt.ylabel('Velocity(km/s)')
                plt.title('PMF FTAN Diagram '+sacname)
        except AttributeError:
            print 'Error: FTAN Parameters are not available!'
        return


class InputFtanParam(object): ###
    """
    A subclass to store input parameters for aftan analysis and SNR Analysis
    -----------------------------------------------------------------------------------------------------
    Parameters:
    pmf         - flag for Phase-Matched-Filtered output (default: Fasle)
    piover4     - phase shift = pi/4*piover4, for cross-correlation piover4 should be -1.0
    vmin        - minimal group velocity, km/s
    vmax        - maximal group velocity, km/s
    tmin        - minimal period, s
    tmax        - maximal period, s
    tresh       - treshold for jump detection, usualy = 10, need modifications
    ffact       - factor to automatic filter parameter, usualy =1
    taperl      - factor for the left end seismogram tapering, taper = taperl*tmax,    (real*8)
    snr         - phase match filter parameter, spectra ratio to determine cutting point for phase matched filter
    fmatch      - factor to length of phase matching window
    fhlen       - half length of Gaussian width
    dosnrflag   - whether to do SNR analysis or not
    predV          - predicted phase velocity curve
    -----------------------------------------------------------------------------------------------------
    """
    def __init__(self):
        self.pmf=False
        self.piover4=-1.0
        self.vmin=1.5
        self.vmax=5.0
        self.tmin=4.0
        self.tmax=30.0
        self.tresh=20.0
        self.ffact=10.0
        self.taperl=1.0
        self.snr=0.2
        self.fmatch=1.0
        self.fhlen=0.008
        self.dosnrflag=False
        self.predV=np.array([])


class specfem2dASDF(pyasdf.ASDFDataSet):

    def readtxt(self, stafile, datadir, verbose=True):
        """ Read txt seismogram files into ASDF dataset according to given station list
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        stafile        - station list file name
        datadir      - data directory
        Output:
        self.waveforms
        -----------------------------------------------------------------------------------------------------
        """
        print 'Reading stationxml!'
        SLst=stations.StaLst()
        SLst.ReadStaList(stafile=stafile)
        StaInv=SLst.GetInventory() 
        self.add_stationxml(StaInv)
        print 'Start reading txt files!'
        for sta in SLst.stations:
            txtfname = datadir+'/'+sta.network+'.'+sta.stacode+'.'+'BXY.semd'
            if verbose==True:
                print 'Reading',txtfname
            inArr=np.loadtxt(txtfname)
            time=inArr[:, 0]
            data=inArr[:, 1]
            tr=obspy.core.Trace()
            tr.data=data
            tr.stats.station=sta.stacode
            tr.stats.delta=time[1]-time[0]
            tr.stats.network=sta.network
            self.add_waveforms(tr, tag='mem2d_raw')
        print 'End reading txt files!'
        return
    
    def AddEvent(self, x, z):
        """ Add event information to ASDF dataset
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        x,z       - event location, unit is km
        Output:
        self.events
        -----------------------------------------------------------------------------------------------------
        """
        print 'Attention: Event Location unit is km!'
        origin=obspy.core.event.origin.Origin(longitude=x, latitude=z, depth=0.)
        event=obspy.core.event.event.Event(origins=[origin])
        catalog=obspy.core.event.catalog.Catalog(events=[event])
        self.add_quakeml(catalog)
        return
    
    def aftan(self, tb=0.0, outdir=None, inftan=InputFtanParam(), vph=3.0, basic1=True, basic2=False,
            pmf1=False, pmf2=False):
        """ aftan analysis for ASDF Dataset
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        compindex  - component index in waveforms path (default = 0)
        tb                 -  begin time (default = 0)
        outdir          - directory for output disp txt files (default = None, no txt output)
        inftan          - input aftan parameters
        phvelname  - predicted phase velocity file (default = './ak135.disp' )
        basic1          - save basic aftan results or not
        basic2          - save basic aftan results(with jump correction) or not
        pmf1            - save pmf aftan results or not
        pmf2            - save pmf aftan results(with jump correction) or not
        
        Output:
        self.auxiliary_data.DISPbasic1, self.auxiliary_data.DISPbasic2, self.auxiliary_data.DISPpmf1, self.auxiliary_data.DISPpmf2
        -----------------------------------------------------------------------------------------------------
        """
        print 'Start aftan analysis!'
        try:
            evlo=self.events.events[0].origins[0].longitude
            evla=self.events.events[0].origins[0].latitude
        except:
            raise ValueError('No event specified to the datasets!')
        ###
        per=np.arange(10)*5.+5.
        v=np.ones(10)*vph
        predV=np.append(per, v)
        predV=predV.reshape(2,10)
        predV=predV.T
        ###
        for station_id in self.waveforms.list():
            # Get data from ASDF dataset
            tr=self.waveforms[station_id].mem2d_raw[0]
            tr.stats.sac={}
            tr.stats.sac.evlo=evlo
            tr.stats.sac.evla=evla
            tr.stats.sac.b=tb
            stlo=self.waveforms[station_id].coordinates['longitude']
            stla=self.waveforms[station_id].coordinates['latitude']
            tr.stats.sac.stlo=stlo*100. # see stations.StaLst.GetInventory
            tr.stats.sac.stla=stla*100.
            # aftan analysis
            ntrace=specfem2dtrace(tr.data, tr.stats)
            ntrace.aftan(pmf=inftan.pmf, piover4=inftan.piover4, vmin=inftan.vmin,
                vmax=inftan.vmax, tmin=inftan.tmin, tmax=inftan.tmax, tresh=inftan.tresh,
                ffact=inftan.ffact, taperl=inftan.taperl, snr=inftan.snr, fmatch=inftan.fmatch, predV=predV)
            print 'aftan analysis for', station_id#, ntrace.stats.sac.dist
            station_id_aux=tr.stats.network+tr.stats.station # station_id for auxiliary data("SW4AAA"), not the diference with station_id "SW4.AAA"
            # save aftan results to ASDF dataset
            if basic1==True:
                parameters={'Tc': 0, 'To': 1, 'Vgr': 2, 'Vph': 3, 'ampdb': 4, 'dis': 5, 'snrdb': 6, 'mhw': 7, 'amp': 8, 'Np': ntrace.ftanparam.nfout1_1,
                        'knetwk': tr.stats.network, 'kstnm': tr.stats.station}
                self.add_auxiliary_data(data=ntrace.ftanparam.arr1_1, data_type='DISPbasic1', path=station_id_aux, parameters=parameters)
            if basic2==True:
                parameters={'Tc': 0, 'To': 1, 'Vgr': 2, 'Vph': 3, 'ampdb': 4, 'snrdb': 5, 'mhw': 6, 'amp': 7, 'Np': ntrace.ftanparam.nfout2_1,
                        'knetwk': tr.stats.network, 'kstnm': tr.stats.station}
                self.add_auxiliary_data(data=ntrace.ftanparam.arr2_1, data_type='DISPbasic2', path=station_id_aux, parameters=parameters)
            if inftan.pmf==True:
                if pmf1==True:
                    parameters={'Tc': 0, 'To': 1, 'Vgr': 2, 'Vph': 3, 'ampdb': 4, 'dis': 5, 'snrdb': 6, 'mhw': 7, 'amp': 8, 'Np': ntrace.ftanparam.nfout1_2,
                        'knetwk': tr.stats.network, 'kstnm': tr.stats.station}
                    self.add_auxiliary_data(data=ntrace.ftanparam.arr1_2, data_type='DISPpmf1', path=station_id_aux, parameters=parameters)
                if pmf2==True:
                    parameters={'Tc': 0, 'To': 1, 'Vgr': 2, 'Vph': 3, 'ampdb': 4, 'snrdb': 5, 'mhw': 6, 'amp': 7, 'Np': ntrace.ftanparam.nfout2_2,
                        'knetwk': tr.stats.network, 'kstnm': tr.stats.station}
                    self.add_auxiliary_data(data=ntrace.ftanparam.arr2_2, data_type='DISPpmf2', path=station_id_aux, parameters=parameters)
            if outdir != None:
                foutPR=outdir+"/"+station_id
                tr.ftanparam.writeDISP(foutPR)
        print 'End aftan analysis!'
        return
    
    def SelectData(self, outfname, stafile, sacflag=True, compindex=np.array([0]), data_type='DISPbasic1' ):
        """ Select data from ASDF Dataset
        Code need checking
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        outfname    - output ASDF file name
        stafile          -  station list file name
        sacflag         - select sac data or not
        compindex  - component index in waveforms path (default = np.array([0]))
        data_type       - dispersion data type (default = DISPbasic1, basic aftan results)
        
        Output:
        Ndbase
        -----------------------------------------------------------------------------------------------------
        """
        SLst = stations.StaLst()
        SLst.ReadStaList(stafile=stafile)
        StaInv=SLst.GetInventory()
        Ndbase=sw4ASDF(outfname)
        Ndbase.add_stationxml(StaInv)
        Ndbase.add_quakeml(self.events)
        disptypelst=['DISPbasic1', 'DISPbasic2', 'DISPpmf1', 'DISPpmf2']
        for sta in SLst.stations:
            station_id=sta.network+'.'+sta.stacode
            station_id_aux=sta.network+sta.stacode
            if sacflag==True:
                try:
                    for cindex in compindex:
                        tr=self.waveforms[station_id].mem2d_raw[cindex]
                        Ndbase.add_waveforms(tr, tag='mem2d_raw')
                except:
                    print 'No sac data for:',station_id,'!'
            if data_type!='All' and data_type !='all':
                try:
                    data=self.auxiliary_data[data_type][station_id_aux].data.value
                    parameters=self.auxiliary_data[data_type][station_id_aux].parameters
                    # data=self.auxiliary_data[disptype][sta.stacode].data.value
                    # parameters=self.auxiliary_data[disptype][sta.stacode].parameters
                    Ndbase.add_auxiliary_data(data=data, data_type=data_type, path=station_id_aux, parameters=parameters)
                except:
                    print 'No', data_type, 'data for:', station_id, '!'
            else:
                for dispindex in disptypelst:
                    try:
                        data=self.auxiliary_data[dispindex][station_id_aux].data.value
                        parameters=self.auxiliary_data[dispindex][station_id_aux].parameters
                        # data=self.auxiliary_data[disptype][sta.stacode].data.value
                        # parameters=self.auxiliary_data[disptype][sta.stacode].parameters
                        Ndbase.add_auxiliary_data(data=data, data_type=dispindex, path=station_id_aux, parameters=parameters)
                    except:
                        print 'No', dispindex, 'data for:', station_id, '!'
        return Ndbase
    
    def InterpDisp(self, data_type='DISPbasic1', pers=np.array([10., 15., 20., 25.]), verbose=True):
        # outindex={'To': 0, 'Vgr': 1, 'Vph': 2,  'amp': 3, 'Np': pers.size}
        staidLst=self.auxiliary_data[data_type].list()
        for staid in staidLst:
            knetwk=str(self.auxiliary_data[data_type][staid].parameters['knetwk'])
            kstnm=str(self.auxiliary_data[data_type][staid].parameters['kstnm'])
            if verbose==True:
                print 'Interpolating dispersion curve for '+ knetwk + kstnm
            outindex={ 'To': 0, 'Vgr': 1, 'Vph': 2,  'amp': 3, 'inbound': 4, 'Np': pers.size, 'knetwk': knetwk, 'kstnm': kstnm }
            data=self.auxiliary_data[data_type][staid].data.value
            index=self.auxiliary_data[data_type][staid].parameters
            Np=index['Np']
            if Np < 5:
                print 'Not enough datapoints for: '+ knetwk+'.'+kstnm
            obsT=data[index['To']][:Np]
            Vgr=np.interp(pers, obsT, data[index['Vgr']][:Np] )
            Vph=np.interp(pers, obsT, data[index['Vph']][:Np] )
            amp=np.interp(pers, obsT, data[index['amp']][:Np] )
            inbound=(pers > obsT[0])*(pers < obsT[-1])*1
            interpdata=np.append(pers, Vgr)
            interpdata=np.append(interpdata, Vph)
            interpdata=np.append(interpdata, amp)
            interpdata=np.append(interpdata, inbound)
            interpdata=interpdata.reshape(5, pers.size)
            self.add_auxiliary_data(data=interpdata, data_type=data_type+'interp', path=staid, parameters=outindex)
        return
    
    def GetField(self, data_type='DISPbasic1', fieldtype='Vgr', pers=np.array([10.]), outdir=None, distflag=True, verbose=False ):
        ### Need Check
        print 'Getting field data!'
        data_type=data_type+'interp'
        tempdict={'Vgr': 'Tgr', 'Vph': 'Tph', 'amp': 'Amp'}
        if distflag==True:
            outindex={ 'x': 0, 'y': 1, tempdict[fieldtype]: 2,  'dist': 3 }
        else:
            outindex={ 'x': 0, 'y': 1, tempdict[fieldtype]: 2 }
        staidLst=self.auxiliary_data[data_type].list()
        evlo=self.events.events[0].origins[0].longitude
        evla=self.events.events[0].origins[0].latitude
        for per in pers:
            FieldArr=np.array([])
            Nfp=0
            for staid in staidLst:
                data=self.auxiliary_data[data_type][staid].data.value # Get interpolated aftan data
                index=self.auxiliary_data[data_type][staid].parameters # Get index
                knetwk=str(self.auxiliary_data[data_type][staid].parameters['knetwk'])
                kstnm=str(self.auxiliary_data[data_type][staid].parameters['kstnm'])
                if verbose==True:
                    print 'Getting field data', knetwk, kstnm
                station_id=knetwk+'.'+kstnm
                obsT=data[index['To']]
                outdata=data[index[fieldtype]]
                inbound=data[index['inbound']]
                fieldpoint=outdata[obsT==per]
                if fieldpoint == np.nan or fieldpoint==0:
                    print station_id+' has nan/zero value'+' T='+str(per)+'s'
                    continue
                # print fieldpoint
                inflag=inbound[obsT==per]
                if fieldpoint.size==0:
                    print 'No datapoint for'+ station_id+' T='+per+'s in interpolated disp dataset!'
                    continue
                if inflag == 0:
                    print 'Datapoint out of bound: '+ knetwk+'.'+kstnm+' T='+str(per)+'s!'
                    continue
                stlo=self.waveforms[station_id].coordinates['longitude']*100.
                stla=self.waveforms[station_id].coordinates['latitude']*100.
                distance=np.sqrt( (stlo-evlo)**2 + (stla-evla)**2 )
                if distance == 0:
                    continue
                FieldArr=np.append(FieldArr, stlo)
                FieldArr=np.append(FieldArr, stla)
                if fieldtype=='Vgr' or fieldtype=='Vph':
                    fieldpoint=distance/fieldpoint
                FieldArr=np.append(FieldArr, fieldpoint)
                if distflag==True:
                    FieldArr=np.append(FieldArr, distance)
                Nfp+=1
            if distflag==True:
                FieldArr=FieldArr.reshape( Nfp, 4)
            else:
                FieldArr=FieldArr.reshape( Nfp, 3)
            if outdir!=None:
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)
                txtfname=outdir+'/'+tempdict[fieldtype]+'_'+str(per)+'.txt'
                np.savetxt(txtfname, FieldArr, fmt='%g')
            self.add_auxiliary_data(data=FieldArr, data_type='Field'+data_type, path=tempdict[fieldtype]+str(int(per)), parameters=outindex)
        return
    
    def aftanParallel(self, tb=0.0, outdir=None, inftan=InputFtanParam(), vph=3.0, basic1=True, basic2=False,
            pmf1=False, pmf2=False):
        """
        Code need checking
        """
        print 'Start aftan analysis!'
        try:
            evlo=self.events.events[0].origins[0].longitude
            evla=self.events.events[0].origins[0].latitude
        except:
            raise ValueError('No event specified to the datasets!')
        ###
        per=np.arange(10)*5.+5.
        v=np.ones(10)*vph
        predV=np.append(per, v)
        predV=predV.reshape(2,10)
        predV=predV.T
        ###
        noiseStream=[]
        for station_id in self.waveforms.list():
            # Get data from ASDF dataset
            tr=self.waveforms[station_id].mem2d_raw[0]
            tr.stats.sac={}
            tr.stats.sac.evlo=evlo
            tr.stats.sac.evla=evla
            tr.stats.sac.b=tb
            stlo=self.waveforms[station_id].coordinates['longitude']
            stla=self.waveforms[station_id].coordinates['latitude']
            tr.stats.sac.stlo=stlo*100. # see stations.StaLst.GetInventory
            tr.stats.sac.stla=stla*100.
            # aftan analysis
            ntrace=specfem2dtrace(tr.data, tr.stats)
            noiseStream.append(ntrace)
        AFTAN = partial(aftan4mp, inftan=inftan)
        pool =multiprocessing.Pool(processes=2)
        FLst=pool.map(AFTAN, noiseStream, chunksize=2) #make our results with a map call
        pool.close() #we are not adding any more processes
        pool.join() #tell it to wait until all threads are done before going on
        return
    
def aftanManager():
    m = BaseManager()
    m.start()
    return m

def aftan4mp(nTr, inftan):
    
    nTr.aftan(pmf=inftan.pmf, piover4=inftan.piover4, vmin=inftan.vmin,
                vmax=inftan.vmax, tmin=inftan.tmin, tmax=inftan.tmax, tresh=inftan.tresh,
                ffact=inftan.ffact, taperl=inftan.taperl, snr=inftan.snr, fmatch=inftan.fmatch, predV=inftan.predV)
    # mylock.acquire() # will block if lock is already held
    print 'aftan analysis for', nTr.stats.station#, ntrace.stats.sac.dist 
    # flst.append(nTr.ftanparam)
    # mylock.release()
    return nTr.ftanparam


    