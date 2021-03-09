
% C************************************
% C The interface of the following original fuction is changed to match
% C the interface of the ETL routines.    Yong Han, 1999.
% 
% C      FUNCTION O2ABS(TEMP,PRES,VAPDEN,FREQ)
% C
% C     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
% C              IN NEPERS/KM
% C
% C      5/1/95  P. Rosenkranz 
% C      11/5/97  P. Rosenkranz - 1- line modification.
% c      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
% C
% C     ARGUMENTS:
% C      REAL TEMP,PRES,VAPDEN,FREQ
% C
% C     NAME    UNITS    DESCRIPTION        VALID RANGE
% C
% C     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
% c                                          valid for atmosphere
% C     PRES   MILLIBARS PRESSURE           3 TO 1000
% C     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
% C                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
% C     FREQ    GHZ      FREQUENCY          0 TO 900
% C
% C     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
% C     P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
% C      BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
% C     H.J. Liebe et al, JQSRT V.48, PP.629-643 (1992).
% c     M.J. Schwartz, Ph.D. thesis, M.I.T. (1997).
% C     SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
% c     This version differs from Liebe's MPM92 in two significant respects:
% c     1. It uses the modification of the 1- line width temperature dependence
% c     recommended by Schwartz: (1/T).
% c     2. It uses the same temperature dependence (X) for submillimeter 
% c     line widths as in the 60 GHz band: (1/T)**0.8 
% C************************************
%       subroutine O2N2_xxx (pdrykpa,vx,ekpa,frq,npp,ncpp)

function [npp,ncpp] = o2n2_rosen98_xxx(pdrykpa,vx,ekpa,frq);

% CYH ***************************************************************

% C      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
      F = [118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,...
            59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,...
            56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,...
            55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,...
            53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,...
            52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7632,...
           487.2494, 715.3931, 773.8397, 834.1458];
        
      S300 = [.2936E-14,.8079E-15, .2480E-14,.2228E-14,...
              .3351E-14,.3292E-14, .3721E-14,.3891E-14,...
              .3640E-14,.4005E-14, .3227E-14,.3715E-14,...
              .2627E-14,.3156E-14, .1982E-14,.2477E-14,...
              .1391E-14,.1808E-14, .9124E-15,.1230E-14,...
              .5603E-15,.7842E-15, .3228E-15,.4689E-15,...
              .1748E-15,.2632E-15, .8898E-16,.1389E-15,...
              .4264E-16,.6899E-16, .1924E-16,.3229E-16,...
              .8191E-17,.1423E-16, .6494E-15, .7083E-14, .3025E-14,...
              .1835E-14, .1158E-13, .3993E-14];
          
      BE = [.009, .015, .083, .084, .212, .212, .391, .391, .626, .626,...
            .915, .915, 1.260, 1.260, 1.660, 1.665, 2.119, 2.115, 2.624, 2.625,...
           3.194, 3.194, 3.814, 3.814, 4.484, 4.484, 5.224, 5.224, 6.004, 6.004, 6.844, 6.844,...
           7.744, 7.744, .048, .044, .049, .145, .141, .145];

       % C      WIDTHS IN MHZ/MB
      WB300 = .56; 
      X = .8;

      W300 = [1.63, 1.646, 1.468, 1.449, 1.382, 1.360,...
              1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171,...
              1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 1.05, 1.05,...
              1.02, 1.02, 1.00, 1.00, .97, .97, .94, .94, .92, .92,...
               .89, .89, 1.92, 1.92, 1.92, 1.81, 1.81, 1.81];
              
      Y300 = [-0.0233,  0.2408, -0.3486,  0.5227,...
              -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311,...
               0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599,...
              0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246,...
              0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546,...
              0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529,... 
              0., 0., 0., 0., 0., 0.,];
          
      V = [  0.0079, -0.0978,  0.0844, -0.1273,...
             0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,...
             0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,...
             0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,...
             0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,...
             0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545,...
             0., 0., 0., 0., 0., 0.];


      db2np = log(10.) * 0.1;
      rvap = 0.01 * 8.314510 / 18.01528;
      factor = .182 * frq;
      TEMP = 300./vx;
      PRES = (pdrykpa+ekpa)*10.;
      VAPDEN = ekpa*10./(rvap*TEMP);
      FREQ = frq;

      TH = 300./TEMP;
      TH1 = TH-1.;
      B = TH^X;
      PRESWV = VAPDEN*TEMP/217.;
      PRESDA = PRES -PRESWV;
      DEN = .001*(PRESDA*B + 1.1*PRESWV*TH);
      DENS = .001*(PRESDA + 1.1*PRESWV)*TH;
      DFNR = WB300*DEN;

%      SUM = 1.6E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR));
      SUM = 0.0;

      NLINES = length(F);
      for K = 1:NLINES
         if K == 1 %!exception for 1- line
            DF = W300(1)*DENS;
         else
            DF = W300(K)*DEN;
         end
         Y = .001*PRES*B*(Y300(K)+V(K)*TH1);
         STR = S300(K)*exp(-BE(K)*TH1);
         SF1 = (DF + (FREQ-F(K))*Y)/((FREQ-F(K))^2 + DF*DF);
         SF2 = (DF - (FREQ+F(K))*Y)/((FREQ+F(K))^2 + DF*DF);
         SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))^2;
      end
      
%      O2ABS = .5034E12*SUM*PRESDA*TH^3/3.14159;
      npp = .5034E12*SUM*PRESDA*TH^3/3.14159;
      ncpp = 1.6E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR));
      ncpp = .5034E12*ncpp*PRESDA*TH^3/3.14159;
% C    add N2 term
      ncpp = ncpp + ABSN2(TEMP,PRES,FREQ);
% C     change the units from np/km to ppm
      npp = (npp /db2np)/factor;
      ncpp = (ncpp / db2np)/factor;

return
      
end
      