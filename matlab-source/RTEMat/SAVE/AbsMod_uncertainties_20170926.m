% This function provides the uncertainties affecting absorption model
% coefficients I was able to find in litterature.
% The baseline are the routines of Rosenkranz 2016 + modification to water-2-air by Koshelev et al 2015.
%
% References:
%   Koshelev et al., JQSRT, 154, 24-27, 2015
%   Koshelev et al., in preparation, 2017 (22 GHz)
%   Koshelev et al., JQSRT, 196, 78?86, 2017 (118 GHz)
%   Turner et al., TGRSS, 47, 10, 3326-37, 2009
%   Tretyakov, JMS, 2016
%
% Hystory:
% 2016/12/05 - Nico - First created

% Quotes:

% Phil:
% The parameters that contribute most uncertainty are related to pressure-broadening: widths, mixing, & their T dependence.
%
% The intensities and their temperature dependence are obtained (at least in my model) from HITRAN2012. B2=(Ef+Ei)/2kTo,
% where Ef, Ei=upper and lower energy levels of the line; k=Boltzmann const; To=296K.
% Although HITRAN has uncertainties for six parameters, the energy level is not among them. 
% I expect that it does not contribute significant uncertainty to the model.

%From Verdes et al 2005;
% For the pressure shifts, an estimated uncertainty for all lines belonging 
% to the same molecular species is assumed. This is a poor approximation as 
% it is well known that the pressure shift is usually very different (with 
% possible change of sign) from one given transition to the next. 
% Thus, the investigated pressure shift uncertainties were 300 kHz/Torr for 
% H2O and 50 kHz/Torr for O2 lines
        
% From Makarov et al., 2011
% "The fidelity of the new model to the spectrometer data is generally better than 2% between 54 and 65 GHz."
% "It is therefore important to distinguish between the uncertainty of the calculated absorption, and uncertainties in the values of the model?s coefficients. The coefficients are adjusted to fit measurements, and underestimation of some coefficients can be counterbalanced by excess values in other coefficients when absorption is calculated."

function AMU = AbsMod_uncertainties(dum);

unknown = 0.0; % unknown uncertainty
mb2torr = 0.750062; % mb2torr conversion factor; MHz/Torr -> MHz/mb (the conversion factor is 1/torr2mb, i.e. mb2torr)
useKoshelev2017 = 1; % if 1 uses the values of Koshelev et al. 2017 (in preparation) for gamma_a, gamma_w, delta_a, delta_w.  - only for 22 GHz!

% Water-to-air broadening ratio 
% Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AMU.w2a.value = 1.20; 
AMU.w2a.uncer = 0.05;
AMU.w2a.units = 'unitless';
AMU.w2a.refer = 'Koshelev et al., JQSRT, 2015';
% Stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% WV Continuum 
% Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AMU.con_Cf.value = 5.43e-10;
AMU.con_Cf.uncer = unknown;
AMU.con_Cf.units = '1/(km*(mb^2*GHz^2))';
AMU.con_Cf.refer = 'RTE Rosen98 routine';

AMU.con_Cf_factr.value = 1.105;
AMU.con_Cf_factr.uncer = sqrt(0.098^2 + 0.030^2); % I added the two contributions in Table IV (Turner et al. 2009)
AMU.con_Cf_factr.units = 'unitless';
AMU.con_Cf_factr.refer = 'Turner et al., TGRSS, 2009';
% The value used by P. Rosenkranz is slightly different 
% The multipliers in Turner et al are referenced to Ros 1998 model. 
% Most of the line parameters have been updated since then, causing small changes to absorption in the 
% spectrum windows. Therefore, to keep the new model consistent with Turner's results, it's necessary to 
% modify the multipliers slightly, but mainly in the foreign-continuum; 
AMU.con_Cf_factr.value = 1.0976; % (P. Rosenkranz, personal communication)

AMU.con_Cs.value = 1.80e-8;
AMU.con_Cs.uncer = unknown;
AMU.con_Cs.units = '1/(km*(mb^2*GHz^2))';
AMU.con_Cs.refer = 'RTE Rosen98 routine';

AMU.con_Cs_factr.value = 0.79;
AMU.con_Cs_factr.uncer = sqrt(0.17^2 + 0.06^2); % I added the two contributions in Table IV (Turner et al. 2009)
AMU.con_Cs_factr.units = 'unitless';
AMU.con_Cs_factr.refer = 'Turner et al., TGRSS, 2009';

AMU.con_Xf.value = 3.0;
AMU.con_Xf.uncer = 0.6; % I take the maximum value in Section 3.2
AMU.con_Xf.units = 'unitless';
AMU.con_Xf.refer = 'Tretyakov, JMS, 2016';

AMU.con_Xs.value = 7.5;
AMU.con_Xs.uncer = 0.6; % I take the maximum value in Section 3.2
AMU.con_Xs.units = 'unitless';
AMU.con_Xs.refer = 'Tretyakov, JMS, 2016';
% Stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% WV Line parameters 
% Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Line frequency FL (freq,GHz)
AMU.FL.value = [22.23507985 183.310087 repmat(unknown,1,13)];
AMU.FL.uncer = [ 0.00000005   0.000001 repmat(unknown,1,13)];
AMU.FL.units = 'GHz'; 
AMU.FL.refer = 'Tretyakov, JMS, 2016';

% Line intensity (or strength) S
AMU.S.value = [4.390 774.6 repmat(unknown,1,13)] * 1e-25;
AMU.S.uncer = [0.043   7.7 repmat(unknown,1,13)] * 1e-25; % This is about 1% - thus it also covers for the Speed-dependent (SD) effect (see To_Nico_Aug15.doc)  
AMU.S.units = 'cm/mol'; % Change units, see below
AMU.S.refer = 'Tretyakov, JMS, 2016';
% units in Ros. model are [Hz*cm^2]; the conversion factor is just speed of light in cm (P. Rosenkranz, personal communication)
c_cm = 2.99792458e10;
AMU.S.value = AMU.S.value * c_cm;
AMU.S.uncer = AMU.S.uncer * c_cm;
AMU.S.units = 'Hz*cm^2'; 

% T COEFF. OF INTENSITIES (B2 in Ros model)
% uncertainty?
% From M. Tretyakov: The intensity T-coefficients should be even more accurate (than 1%). 
AMU.B2.value = [2.1440 0.6680 6.1790 1.5410 1.0480 3.5950 5.0480 1.4050 3.5970 2.3790 2.8520 0.1590 2.3910 0.3960 1.4410]; % as in Ros2016 model
AMU.B2.uncer = AMU.B2.value / 100 * 1; % 1% - following Tretyakov's suggestion
AMU.B2.units = 'Unitless';
AMU.B2.refer = 'Rosenkranz, 2016 model + Tretyakov pers. comm.';

% Line shape (i.e. pressure broadening parameters)
AMU.gamma_a.value = [3.63 3.926 repmat(unknown,1,13)];
AMU.gamma_a.uncer = [0.14 0.020 repmat(unknown,1,13)];
AMU.gamma_a.units = 'MHz/Torr';
AMU.gamma_a.refer = 'Tretyakov, JMS, 2016';
if useKoshelev2017 % NB use the values of Koshelev et al. 2017 (in preparation) - only for 22 GHz!
    AMU.gamma_a.value = [3.584 3.926 repmat(unknown,1,13)];
    AMU.gamma_a.uncer = [0.052 0.020 repmat(unknown,1,13)];
    AMU.gamma_a.units = 'MHz/Torr';
    AMU.gamma_a.refer = 'Tretyakov, JMS, 2016 + Koshelev et al., in prep., 2017 (22 GHz)';
end
% units in Ros. model are [MHz/mb]; the conversion factor is 1/torr2mb, i.e. mb2torr
AMU.gamma_a.value = AMU.gamma_a.value*mb2torr;
AMU.gamma_a.uncer = AMU.gamma_a.uncer*mb2torr;
AMU.gamma_a.units = 'MHz/mb';

AMU.gamma_w.value = [17.6 19.7 repmat(unknown,1,13)];
AMU.gamma_w.uncer = [ 0.5  0.5 repmat(unknown,1,13)];
AMU.gamma_w.units = 'MHz/Torr';
AMU.gamma_w.refer = 'Tretyakov, JMS, 2016';
if useKoshelev2017 % NB use the values of Koshelev et al. 2017 (in preparation) - only for 22 GHz!
    AMU.gamma_w.value = [17.707 19.7 repmat(unknown,1,13)];
    AMU.gamma_w.uncer = [ 0.045  0.5 repmat(unknown,1,13)];
    AMU.gamma_w.units = 'MHz/Torr';
    AMU.gamma_w.refer = 'Tretyakov, JMS, 2016 + Koshelev et al., in prep., 2017 (22 GHz)';
end
% units in Ros. model are [MHz/mb]; the conversion factor is 1/torr2mb, i.e. mb2torr
AMU.gamma_w.value = AMU.gamma_w.value*mb2torr;
AMU.gamma_w.uncer = AMU.gamma_w.uncer*mb2torr;
AMU.gamma_w.units = 'MHz/mb';

AMU.n_a.value = [0.70 0.74 repmat(unknown,1,13)];
AMU.n_a.uncer = [0.05 0.03 repmat(unknown,1,13)];
AMU.n_a.units = 'unitless';
AMU.n_a.refer = 'Tretyakov, JMS, 2016';
% Note: these numbers are different wrt Rosenkranz 2016 model ([0.76 0.77])

AMU.n_w.value = [1.2 0.78 repmat(unknown,1,13)];
AMU.n_w.uncer = [0.5 0.08 repmat(unknown,1,13)];
AMU.n_w.units = 'unitless';
AMU.n_w.refer = 'Tretyakov, JMS, 2016';
% Note: these numbers are different wrt Rosenkranz 2016 model ([1.0 0.85])

% Line position 
% Line center linear with pressure frequency shifting coefficients (air)
AMU.delta_a.value = [0.0 -0.096 repmat(unknown,1,13)]; % Note there's a typo on Tretyakov 2016 (-0.96 in Table 4 should be -0.096, see Table 3)
AMU.delta_a.uncer = [0.1  0.010 repmat(unknown,1,13)]; % Note there's a typo on Tretyakov 2016 (0.1 in Table 4 should be 0.01, see above)
AMU.delta_a.units = 'MHz/Torr';
AMU.delta_a.refer = 'Tretyakov, JMS, 2016';
% I replace first value with a number within uncertainty, otherwise 0.0
% would cause the uncertainty of SR (see below) go to NaN
AMU.delta_a.value(1) = 0.001; % totally arbitrary
if useKoshelev2017 % NB use the values of Koshelev et al. 2017 (in preparation) - only for 22 GHz!
    AMU.delta_a.value = [-0.032 -0.096 repmat(unknown,1,13)]; 
    AMU.delta_a.uncer = [ 0.038  0.010 repmat(unknown,1,13)]; 
    AMU.delta_a.units = 'MHz/Torr';
    AMU.delta_a.refer = 'Tretyakov, JMS, 2016 + Koshelev et al., in prep., 2017 (22 GHz)';
end
% this is not used in Ros. model, but units should be consistent with gamma_a, i.e. from [MHz/Torr] to [MHz/mb]
AMU.delta_a.value = AMU.delta_a.value*mb2torr;
AMU.delta_a.uncer = AMU.delta_a.uncer*mb2torr;
AMU.delta_a.units = 'MHz/mb';

% Line center linear with pressure frequency shifting coefficients (water)
AMU.delta_w.value = [-0.4 0.23 repmat(unknown,1,13)];
AMU.delta_w.uncer = [ 1.2 0.03 repmat(unknown,1,13)];
AMU.delta_w.units = 'MHz/Torr';
AMU.delta_w.refer = 'Tretyakov, JMS, 2016';
if useKoshelev2017 % NB use the values of Koshelev et al. 2017 (in preparation) - only for 22 GHz!
    AMU.delta_w.value = [1.099 0.23 repmat(unknown,1,13)];   
    AMU.delta_w.uncer = [0.079 0.03 repmat(unknown,1,13)];
    AMU.delta_w.units = 'MHz/Torr';
    AMU.delta_w.refer = 'Tretyakov, JMS, 2016 + Koshelev et al., in prep., 2017 (22 GHz)';
end
% this is not used in Ros. model, but units should be consistent with gamma_w, i.e. from [MHz/Torr] to [MHz/mb]
AMU.delta_w.value = AMU.delta_w.value*mb2torr;
AMU.delta_w.uncer = AMU.delta_w.uncer*mb2torr;
AMU.delta_w.units = 'MHz/mb';

% Shift to width ratio (this is used in Ros 2016 model instead of delta_a/w)
% In Ros 2016: SR = delta_a/gamma_a
for il = 1:2
    [AMU.SR.value(il),AMU.SR.uncer(il)] = uncertainty_propagation('A/B',AMU.delta_a.value(il),AMU.gamma_a.value(il),AMU.delta_a.uncer(il),AMU.gamma_a.uncer(il));
end
AMU.SR.value = [AMU.SR.value repmat(unknown,1,13)];
AMU.SR.uncer = [AMU.SR.uncer repmat(unknown,1,13)];
AMU.SR.units = 'unitless';
AMU.SR.refer = 'Tretyakov, JMS, 2016 - Eq. 3 and 5';

%

% Stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% O2 Line parameters 
% Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NB: I should try to read HITRAN uncertainty! At least for FL and S !

% Line frequency FL (freq,GHz)
% I take values from Ros 2016 and uncertainty from Tretyakov et al., 2005 (See Table 1) 
% Note Table 1 shows the first 26 lines in Rosenk 2015, while 27+ is the the 28th of Rosenk 2015
% For all higher freq lines I assume 17*1e-6 uncertainty (i.e. the largest value)
AMU.O2FL.value = [118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,...
                   59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,...
                   56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,...
                   55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,...
                   53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,...
                   52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.4310,...
                   50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.7630,...
                   487.2493, 566.8956, 715.3929, 731.1866,...
                   773.8395, 834.1455, 895.0710];
AMU.O2FL.uncer = 1e-6*[ 7 12 6 5 5 5 4 4 4 4 4 6 4 4 5 4 6 4 5 4 7 10 15 5 17 12 12 12 ...
                         repmat(17,1,21)];  % 17*1-6 uncertainty for lines 29-49
AMU.O2FL.units = 'GHz'; 
AMU.O2FL.refer = 'Tretyakov et al., 2005';


% Line intensity (or strength) S
% From M. Tretyakov:
%    Line strengths were measured by Liebe, et al. 1977. 
%    We believe that HITRAN data are accurate to about 1%. 
%    We indirectly confirmed this by measuring N=1- intensity (paper in preparation). 
%    The coincidence with HITRAN is less than 0.5%. 
% Nico: I take values from Ros 2016 with 1% uncertainty
AMU.O2S.value = [0.2906E-14,0.7957E-15,0.2444E-14,0.2194E-14,...
                 0.3301E-14,0.3243E-14,0.3664E-14,0.3834E-14,...
                 0.3588E-14,0.3947E-14,0.3179E-14,0.3661E-14,...
                 0.2590E-14,0.3111E-14,0.1954E-14,0.2443E-14,...
                 0.1373E-14,0.1784E-14,0.9013E-15,0.1217E-14,...
                 0.5545E-15,0.7766E-15,0.3201E-15,0.4651E-15,...
                 0.1738E-15,0.2619E-15,0.8880E-16,0.1387E-15,...
                 0.4272E-16,0.6923E-16,0.1939E-16,0.3255E-16,...
                 0.8301E-17,0.1445E-16,0.3356E-17,0.6049E-17,...
                 0.1280E-17,0.2394E-17,...
                 0.3287E-16,0.6463E-15,0.1334E-16,0.7049E-14,...
                 0.3011E-14,0.1797E-16,0.1826E-14,0.2193E-16,...
                 0.1153E-13,0.3974E-14,0.2512E-16];
AMU.O2S.uncer = AMU.O2S.value / 100;  % 1%
AMU.O2S.units = 'cm2*Hz'; 
AMU.O2S.refer = 'Tretyakov, Personal communication, 2016';


% Temperature exponent for intensity
% From M. Tretyakov:
%    The intensity T-coefficients should be even more accurate (than 1%, see Line intensity above) 
% Nico: I take values from Ros 2016 with 1% uncertainty
AMU.O2BE.value = [ .010,  .014,  .083,  .083,  .207,  .207,  .387,  .387,  .621,  .621,...
                   .910,  .910, 1.255, 1.255, 1.654, 1.654, 2.109, 2.109, 2.618, 2.618,...
                  3.182, 3.182, 3.800, 3.800, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 6.819, 6.819,... 
                  7.709, 7.709, 8.653, 8.653, 9.651, 9.651,...
                   .019,  .048,  .045,  .044,  .049,  .084,  .145,  .136,  .141,  .145,  .201];
AMU.O2BE.uncer = AMU.O2BE.value / 100;  % 1%
AMU.O2BE.units = 'unitless'; 
AMU.O2BE.refer = 'Tretyakov, Personal communication, 2016';


% Line shape (i.e. pressure broadening parameters)
% See Table 5 on Tretyakov et al., 2005 - a3 is (nearly) the same as W300 in Rosenk 2016 [GHz/(1e5 Pa)]
% Only the first line differs, plus Ros has more lines at higher freq 
AMU.O2gamma.value = [1.688, 1.703, 1.513, 1.491, 1.415, 1.408, 1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217, 1.189, 1.174, 1.134, 1.134, 1.089, 1.088, 1.037, 1.038, 0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858, 0.811, 0.811, 0.764, 0.764, 0.717, 0.717, 0.669, 0.669, 2.78, 1.64, 1.64, 1.64, 1.60, 1.60, 1.60, 1.60, 1.62, 1.47, 1.47]; % From Ros 2016
% Table 5 does not report uncertainties; thus I use values on Table 1 and converts units
% Note that in Table 1 are the first 26 lines in Rosenk 2015, while 27+ is the the 28th of Rosenk 2015 - below I assumed 27- is the same as 27+
% For all higher freq lines I assume 30*1e-3 uncertainty
% units in Ros. model are [MHz/mb]; the conversion factor is 1/torr2mb, i.e. mb2torr
%AMU.O2gamma.uncer = 1e-3*[10 17 10 10 10 10 8 12 15 10 10 10 10 10 10 15 10 10 15 15 20 12 40 20 25 20 30 30 ...
%                      repmat(30,1,21)];  % 0.03 uncertainty for lines 29-49
%AMU.O2gamma.uncer = AMU.O2gamma.uncer*mb2torr;
% broadnening parameter uncertainty (see Table 1 of Tretyakov)
% I changed this to cope with Phil's suggestion:
% the uncertainties from table 1 of Tretyakov et al (2005) have to be combined; I did that for the 20 strong lines, and the result is listed in table 1 in my memo of July 11 (O2_off-diagonal.docx).
c4unc = 1e-3*[10 17 10 10 10 10 8 12 15 10 10 10 10 10 10 15 10 10 15 15]; % O2 - unc of col 4 in Table 1
c5unc = 1e-3*[20 30 25 15 15 12 10 10 13 30 10 10 17 15 10 10 18 10 40 15]; % N2 - unc of col 5 in Table 1
unc_comb = (0.21*c4unc + 0.79*c5unc); % combined (air) uncertainty [MHz/Torr]
AMU.O2gamma.uncer =  = unc_comb * mb2torr; % MHz/Torr -> MHz/mb == GHz/bar (this is the same as in table 1 of O2_off-diagonal.docx)
AMU.O2gamma.units = 'MHz/mb';
AMU.O2gamma.refer = 'Tretyakov, JMS, 2005';
% From M. Tretyakov: Additionally you may use Koshelev et al., JQSRT 2016 (Table 1, gamma_self)
% Nico: measurement by Koshelev have smaller uncertainty, so are covered

% DA RIVEDERE !
% After suggestion of Phil Rosenkranz (2017/07/17-18) for Weak Lines (WL)
% to test 0.05 GHz/bar uncertainty on weaker lines
% It's added to O2gamma, so its center value is zero
%AMU.O2gamma_WL.value = zeros(size(AMU.O2gamma.value));
%AMU.O2gamma_WL.uncer = zeros(size(AMU.O2gamma.uncer)); 
%AMU.O2gamma_WL.uncer(21:38) = 0.05; % 0.05 GHz/bar = 0.05 MHz/mb uncertainty for weaker lines (21-38)
AMU.O2gamma_WL.value = AMU.O2gamma.value;
AMU.O2gamma_WL.uncer = AMU.O2gamma.uncer; 
%AMU.O2gamma_WL.uncer([21:38]) = 0.0; % suppress uncertainty of weaker lines
AMU.O2gamma_WL.uncer([1:20 39:end]) = 0.0; % suppress uncertainty of stronger lines
AMU.O2gamma_WL.units = 'MHz/mb';
AMU.O2gamma_WL.refer = 'Rosenkranz, pers. comm., 2017';
% DA RIVEDERE !

% O2 line mixing temperature dependance
% Table 5, last column of Tretyakov 2005. This is called V in MPM
AMU.O2_V.value = [  0.0079, -0.0978,  0.0844, -0.1273,...
                    0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,...
                    0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,...
                    0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,...
                    0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,...
                    0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545,...
                    0.680,  -0.660,   0.685,  -0.665,...
                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.];
AMU.O2_V.uncer = AMU.O2_V.value * 0.2; % Phil suggested "perturb all of them by 20%, simultaneously", see email on 2017/07/18
AMU.O2_V.units = 'unitless';
AMU.O2_V.refer = 'Tretyakov et al., 2005. Uncertainty from Rosenkranz, pers. comm., 2017';


% Temperature dependence of broadening coefficient (as in Makarov et al. 2011; 2008)
AMU.X11.value = 0.785;
AMU.X11.uncer = 0.035;
AMU.X11.units = 'unitless';
AMU.X11.refer = 'Makarov et al. JQSRT 2011 -> Makarov et al. JQSRT 2008';
% same parameter but for values in Koshelev et al., 2016 (Table 1, n_gamma_self)
AMU.X16.value = 0.765;
AMU.X16.uncer = 0.011;
AMU.X16.units = 'unitless';
AMU.X16.refer = 'Koshelev et al., 2016';
% same parameter but for values in Tretyakov et al., 2005 
% we assume 0.05 uncertainty to fit with latest results of Koshelev et al. 2016 
% (see email from Phil of 2017/07/17)
AMU.X05.value = 0.80;
AMU.X05.uncer = 0.05;
AMU.X05.units = 'unitless';
AMU.X05.refer = 'Tretyakov et al., 2005 + Koshelev et al., 2016';

% Line mixing (2% lump uncertainty?)
% Lump absorption percentage uncertainty (APU) due to line mixing parameters
AMU.APU.value = 1.00; % 
%AMU.APU.uncer = 0.02; % if not 1.00, indicates percentage value (e.g. 0.02 = 2%)
AMU.APU.uncer = 1.00; % if 1.00, uncertainty is frq/T/P dependent, see APU_line_mixing.m
AMU.APU.units = 'unitless';
AMU.APU.refer = 'Makarov et al. JQSRT 2011';
% From Makarov et al. 2011
% Systematic discrepancies between the extended model and experiment are still noticeable, 
% especially at lower temperatures, where the mixing effect should be stronger. 
% These discrepancies are similar for records repeated at the same experimental conditions, 
% so they are not caused by random experimental errors.

% Stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
