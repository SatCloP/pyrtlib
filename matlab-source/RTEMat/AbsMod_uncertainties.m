% This function provides the uncertainties affecting absorption model
% coefficients I was able to find in litterature.
% The baseline are the routines of Rosenkranz 2016 + modification to water-2-air by Koshelev et al 2015.
%
% References:
%   Koshelev et al., JQSRT, 112, 2704?2712, 2011
%   Koshelev et al., JQSRT, 154, 24-27, 2015
%   Koshelev et al., in preparation, 2018 (22 GHz)
%   Koshelev et al., JQSRT, 196, 78?86, 2017 (118 GHz)
%   Turner et al., TGRSS, 47, 10, 3326-37, 2009
%   Tretyakov, JMS, 2016
%
% Hystory:
% 2016/12/05 - Nico - First created
% 2018/12/19 - Nico - Modified for adding 658 GHz line as in Rosenkranz 2018 (see if strcmp(mdl,'r18') at the end of WV parameters)
% 2020/06/25 - Nico - Modified to account SD 183 parameters only (see if strcmp(mdl,'r20sd'))

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

function AMU = AbsMod_uncertainties(mdl);

unknown = 0.0; % unknown uncertainty
mb2torr = 0.750062; % mb2torr conversion factor; MHz/Torr -> MHz/mb (the conversion factor is 1/torr2mb, i.e. mb2torr)
useKoshelev2017 = 1; % if 1 uses the values of Koshelev et al. 2018 (see also .doc draft from 2017) for gamma_a, gamma_w, delta_a, delta_w.  - only for 22 GHz!
useKoshelev2017_what = 'RAD'; % either 'RAD', 'video', or 'combo'

% Water-to-air broadening ratio 
% Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AMU.w2a.value = 1.20; 
AMU.w2a.uncer = 0.05;
AMU.w2a.units = 'unitless';
AMU.w2a.refer = 'Koshelev et al., JQSRT, 2015';
% Stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% WV Continuum 
% Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AMU.con_Cf.value = 5.43e-10; % [1/(km*(mb^2*GHz^2))] - divide by db2np*1e-2 to compare with Tretyakov 2016, which uses [dB/(km*(kPa^2*GHz^2))]
AMU.con_Cf.uncer = unknown;
AMU.con_Cf.units = '1/(km*(mb^2*GHz^2))';
AMU.con_Cf.refer = 'RTE Rosen98 routine';

AMU.con_Cf_factr.value = 1.105;
AMU.con_Cf_factr.uncer = sqrt(0.098^2 + 0.030^2); % I added the two contributions in Table IV (Turner et al. 2009)
%AMU.con_Cf_factr.uncer = AMU.con_Cf_factr.uncer * 1.214; % see last paragraph of H2O_cross-covariance_rev.docx (5.56*1.214=6.75)
AMU.con_Cf_factr.uncer = AMU.con_Cf_factr.uncer; % I apply the correction outside, only when it's needed
AMU.con_Cf_factr.units = 'unitless';
AMU.con_Cf_factr.refer = 'Turner et al., TGRSS, 2009';
% The value used by P. Rosenkranz is slightly different 
% The multipliers in Turner et al are referenced to Ros 1998 model. 
% Most of the line parameters have been updated since then, causing small changes to absorption in the 
% spectrum windows. Therefore, to keep the new model consistent with Turner's results, it's necessary to 
% modify the multipliers slightly, but mainly in the foreign-continuum; 
AMU.con_Cf_factr.value = 1.0976; % (P. Rosenkranz, personal communication, to account for line updates)

AMU.con_Cs.value = 1.80e-8; % [1/(km*(mb^2*GHz^2))] - divide by db2np*1e-2 to compare with Tretyakov 2016, which uses [dB/(km*(kPa^2*GHz^2))]
AMU.con_Cs.uncer = unknown;
AMU.con_Cs.units = '1/(km*(mb^2*GHz^2))';
AMU.con_Cs.refer = 'RTE Rosen98 routine';

AMU.con_Cs_factr.value = 0.79;
AMU.con_Cs_factr.uncer = sqrt(0.17^2 + 0.06^2); % I added the two contributions in Table IV (Turner et al. 2009)
AMU.con_Cs_factr.units = 'unitless';
AMU.con_Cs_factr.refer = 'Turner et al., TGRSS, 2009';

AMU.con_Xf.value = 3.0;
%AMU.con_Xf.uncer = 0.6; % I take the maximum value in Section 3.2
AMU.con_Xf.uncer = 0.8; % to have some overlap with value and uncertainty given by Koshelev et al 2011, i.e. 0.91(17)+3.0 = 3.91(17)
AMU.con_Xf.units = 'unitless';
AMU.con_Xf.refer = 'Tretyakov, JMS, 2016; Koshelev et al. 2011';

AMU.con_Xs.value = 7.5;
AMU.con_Xs.uncer = 0.6; % I take the maximum value in Section 3.2, which overlaps with the value given by Koshelev et al 2011, i.e. 5.24(21)+3 = 8.24(21)
%AMU.con_Xs.uncer = 3.1; % this value reaches 10.6, which works better at low temperatures (see Section 3.2, page 17, 1st column)
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
% Actually, Misha (see AMU_paper_Tables_Figures_References_V03m.docx) say the uncertainty is <1e-6 % for 22 and 183 GHz (refer to Tennyson et al., 2013)
AMU.B2.value = [2.1440 0.6680 6.1790 1.5410 1.0480 3.5950 5.0480 1.4050 3.5970 2.3790 2.8520 0.1590 2.3910 0.3960 1.4410]; % as in Ros2016 model
AMU.B2.uncer = AMU.B2.value / 100 * 1; % 1% - this should be a large overestimation (according to Tretyakov)
AMU.B2.units = 'Unitless';
AMU.B2.refer = 'Rosenkranz, 2016 model + Tretyakov pers. comm.';

% Line shape (i.e. pressure broadening parameters)
AMU.gamma_a.value = [3.63 3.926 repmat(unknown,1,13)];
AMU.gamma_a.uncer = [0.14 0.020 repmat(unknown,1,13)];
AMU.gamma_a.units = 'MHz/Torr';
AMU.gamma_a.refer = 'Tretyakov, JMS, 2016';
if useKoshelev2017 % NB use the values of Koshelev et al. 2017 (in preparation) - only for 22 GHz!
    switch useKoshelev2017_what
        case 'RAD'
            AMU.gamma_a.value(1) = 3.598;
            AMU.gamma_a.uncer(1) = 0.022;
            AMU.gamma_a.units = 'MHz/Torr';
            AMU.gamma_a.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)';            
        case 'video'
            AMU.gamma_a.value(1) = 3.397;
            AMU.gamma_a.uncer(1) = 0.079;
            AMU.gamma_a.units = 'MHz/Torr';
            AMU.gamma_a.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)';                        
        case 'combo'
            AMU.gamma_a.value(1) = 3.584;
            AMU.gamma_a.uncer(1) = 0.052;
            AMU.gamma_a.units = 'MHz/Torr';
            AMU.gamma_a.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)';
    end
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
    switch useKoshelev2017_what
        case 'RAD'
            AMU.gamma_w.value(1) = 17.713;
            AMU.gamma_w.uncer(1) =  0.015;
            AMU.gamma_w.units = 'MHz/Torr';
            AMU.gamma_w.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)';            
        case 'video'
            AMU.gamma_w.value(1) = 17.35;
            AMU.gamma_w.uncer(1) =  0.12;
            AMU.gamma_w.units = 'MHz/Torr';
            AMU.gamma_w.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)';                        
        case 'combo'
            AMU.gamma_w.value(1) = 17.707;
            AMU.gamma_w.uncer(1) =  0.045;
            AMU.gamma_w.units = 'MHz/Torr';
            AMU.gamma_w.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)';
    end
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
if useKoshelev2017 % NB use the values of Koshelev et al. 2018 (see also .doc draft from 2017) - only for 22 GHz!
    switch useKoshelev2017_what
        case 'RAD'
            AMU.delta_a.value(1) = -0.044;
            AMU.delta_a.uncer(1) =  0.005;
            AMU.delta_a.units = 'MHz/Torr';
            AMU.delta_a.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)';            
        case 'video'
            AMU.delta_a.value(1) =  0.081;
            AMU.delta_a.uncer(1) =  0.015;
            AMU.delta_a.units = 'MHz/Torr';
            AMU.delta_a.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)';                        
        case 'combo'
            AMU.delta_a.value(1) = -0.032;
            AMU.delta_a.uncer(1) =  0.038;
            AMU.delta_a.units = 'MHz/Torr';
            AMU.delta_a.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)';
    end
end
% this is not used in Ros models (before Ros18), but units should be consistent with gamma_a, i.e. from [MHz/Torr] to [MHz/mb]
AMU.delta_a.value = AMU.delta_a.value*mb2torr;
AMU.delta_a.uncer = AMU.delta_a.uncer*mb2torr;
AMU.delta_a.units = 'MHz/mb';

% Line center linear with pressure frequency shifting coefficients (water)
AMU.delta_w.value = [-0.4 0.23 repmat(unknown,1,13)];
AMU.delta_w.uncer = [ 1.2 0.03 repmat(unknown,1,13)];
AMU.delta_w.units = 'MHz/Torr';
AMU.delta_w.refer = 'Tretyakov, JMS, 2016';
if useKoshelev2017 % NB use the values of Koshelev et al. 2017 (in preparation) - only for 22 GHz!
    switch useKoshelev2017_what
        case 'RAD'
            AMU.delta_w.value(1) = 1.085;
            AMU.delta_w.uncer(1) = 0.011;
            AMU.delta_w.units = 'MHz/Torr';
            AMU.delta_w.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)';            
        case 'video'
            AMU.delta_w.value(1) = 1.530;
            AMU.delta_w.uncer(1) = 0.060;
            AMU.delta_w.units = 'MHz/Torr';
            AMU.delta_w.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)';                        
        case 'combo'
            AMU.delta_w.value(1) = 1.099;
            AMU.delta_w.uncer(1) = 0.079;
            AMU.delta_w.units = 'MHz/Torr';
            AMU.delta_w.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)';
    end
end
% this is not used in Ros models (before Ros18), but units should be consistent with gamma_a, i.e. from [MHz/Torr] to [MHz/mb]
AMU.delta_w.value = AMU.delta_w.value*mb2torr;
AMU.delta_w.uncer = AMU.delta_w.uncer*mb2torr;
AMU.delta_w.units = 'MHz/mb';

% Shift to width ratio (this is used in Ros 2016 & 2017 models instead of delta_a/w)
% In Ros 2016: SR = delta_a/gamma_a
for il = 1:2
    [AMU.SR.value(il),AMU.SR.uncer(il)] = uncertainty_propagation('A/B',AMU.delta_a.value(il),AMU.gamma_a.value(il),AMU.delta_a.uncer(il),AMU.gamma_a.uncer(il));
end
AMU.SR.value = [AMU.SR.value repmat(unknown,1,13)];
AMU.SR.uncer = [AMU.SR.uncer repmat(unknown,1,13)];
AMU.SR.units = 'unitless';
AMU.SR.refer = 'Tretyakov, JMS, 2016 - Eq. 3 and 5';

% Intensity temperature-dependence exponent
% See Partition_sum.doc
AMU.wv_nS.value = 2.5;
%AMU.wv_nS.uncer = AMU.wv_nS.value * 0.06; % value computed by Tretyakov (Partition_sum.doc)
AMU.wv_nS.uncer = AMU.wv_nS.value * 0.005; % value computed by Rosenkranz (email of 2018/03/08 and h2o_partition_sum.txt)
AMU.wv_nS.units = 'unitless';
AMU.wv_nS.refer = 'Rosenkranz, email of 2018/03/08, based on Gamache et al., 2017';
%

if strcmp(mdl,'r18')
   
   % r18 has one more line at 658 GHz 
   
   % Read the list of parameters
   h2o_list_r18
   
   % Continuum parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   AMU.con_Cf.value = 5.96e-10; % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Cf.value = CF;       % this is from Ros 2018 (see h2o_list_r18)
   AMU.con_Cf.uncer = 5.50e-11; % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Cf.units = '1/(km*(mb^2*GHz^2))';
   AMU.con_Cf.refer = 'Rosenkranz, 2018; Cimini et al ACP 2018';

   AMU.con_Cs.value = 1.42e-8;  % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Cs.value = CS;       % this is from Ros 2018 (see h2o_list_r18)
   AMU.con_Cs.uncer = 3.20e-9;  % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Cs.units = '1/(km*(mb^2*GHz^2))';
   AMU.con_Cs.refer = 'Rosenkranz, 2018; Cimini et al ACP 2018';
 
   AMU.con_Xf.value = XCF;      % this is from Ros 2018 (see h2o_list_r18)
   AMU.con_Xf.uncer = 0.8;      % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Xf.units = 'unitless';
   AMU.con_Xf.refer = 'Tretyakov, JMS, 2016; Koshelev et al. 2011'; % I keep the original references
   
   %AMU.con_Xs.value = XCF;      % this is from Ros 2018 (see h2o_list_r18)
   AMU.con_Xs.value = XCS;      % NB: Bug fixed on 2020/07/02 !!!!
   AMU.con_Xs.uncer = 0.6;      % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Xs.units = 'unitless';
   AMU.con_Xs.refer = 'Tretyakov, JMS, 2016'; % I keep the original references

   % With r18 I perturb the continuum parameters directly, not the factor
   AMU.con_Cf_factr.value = 1; % 
   AMU.con_Cf_factr.uncer = 0; % 
   AMU.con_Cf_factr.units = 'unitless';
   AMU.con_Cf_factr.refer = 'With R18 I perturb the continuum parameters, not the factor';

   AMU.con_Cs_factr.value = 1; %
   AMU.con_Cs_factr.uncer = 0; % 
   AMU.con_Cs_factr.units = 'unitless';
   AMU.con_Cs_factr.refer = 'With R18 I perturb the continuum parameters, not the factor';
   
   % Lines' parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   NL = length(FL);
   
   % Line frequency FL (freq,GHz)
   AMU.FL.value = FL;                           % this is from Ros 2018 (see h2o_list_r18)
   AMU.FL.uncer = repmat(unknown,1,NL);
   AMU.FL.uncer(1:2) = [0.00000005   0.000001]; % this is from Cimini et al ACP 2018, Table 1
   AMU.FL.units = 'GHz';
   AMU.FL.refer = 'Tretyakov, JMS, 2016';       % I keep the original references
   
   % Line intensity (or strength) S
   AMU.S.value = S1;                     % this is from Ros 2018 (see h2o_list_r18)
   AMU.S.uncer = repmat(unknown,1,NL);
   AMU.S.uncer(1:2) = S1(1:2) / 100;     % 1% is from Cimini et al ACP 2018, Table 1
   AMU.S.units = 'Hz*cm^2';
   AMU.S.refer = 'Tretyakov, JMS, 2016'; % I keep the original references
   
   % T COEFF. OF INTENSITIES (B2 in Ros model)
   AMU.B2.value = B2;                     % this is from Ros 2018 (see h2o_list_r18)
   AMU.B2.uncer = B2 / 100 * 0.5;         % 0.5% is from Cimini et al ACP 2018, Table 1
   AMU.B2.units = 'Unitless';
   AMU.B2.refer = 'Rosenkranz, 2016 + Tretyakov pers. comm.'; % I keep the original references

   % Line shape (i.e. pressure broadening parameters)
   % AIR-BROADENED WIDTH PARAMETERS AT REFTLINE [MHz/mb -> GHz/mb]
   AMU.gamma_a.value = W3 * 1e3;                 % this is from Ros 2018 (see h2o_list_r18) [GHz/mb -> MHz/mb]
   AMU.gamma_a.uncer = repmat(unknown,1,NL); 
   %AMU.gamma_a.uncer(1:2) = [0.039 0.015] * 1e3; % this is from Cimini et al ACP 2018, Table 1 [GHz/mb -> MHz/mb]
   AMU.gamma_a.uncer(1:2) = [0.039 0.015]; % this is from Cimini et al ACP 2018, Table 1 [GHz/bar == MHz/mb]
   AMU.gamma_a.units = 'MHz/mb';
   AMU.gamma_a.refer = 'Rosenkranz, 2018; Cimini et al ACP 2018';
   
   % T-EXPONENT OF AIR-BROADENING [unitless]
   % Note: X at 22 and 183 is different wrt Tretyakov 2016 ([0.70 0.74]), which have been used in Cimini et al. ACP 2018
   AMU.n_a.value = X;                       % this is from Ros 2018 (see h2o_list_r18)
   AMU.n_a.uncer = repmat(unknown,1,NL);
   AMU.n_a.uncer(1:2) = [0.05 0.03];        % this is from Cimini et al ACP 2018, Table 1
   AMU.n_a.units = 'unitless';
   AMU.n_a.refer = 'Rosenkranz, 2018';
   
   % Line shape (i.e. pressure broadening parameters)
   % SELF-BROADENED WIDTH PARAMETERS AT REFTLINE [MHz/mb -> GHz/mb]
   AMU.gamma_w.value = WS * 1e3;                 % this is from Ros 2018 (see h2o_list_r18) [GHz/mb -> MHz/mb]
   AMU.gamma_w.uncer = repmat(unknown,1,NL); 
   %AMU.gamma_w.uncer(1:2) = [0.034 0.015] * 1e3; % this is from Cimini et al ACP 2018, Table 1 [GHz/mb -> MHz/mb]
   AMU.gamma_w.uncer(1:2) = [0.034 0.015]; % this is from Cimini et al ACP 2018, Table 1 [GHz/bar == MHz/mb]
   % I think there is a typo in ACP paper (0.039 should be 0.034), i.e. from 0.045 (combo uncertainty of Koshelev et al. 2018) * mb2torr = 0.0338
   AMU.gamma_w.units = 'MHz/mb';
   AMU.gamma_w.refer = 'Rosenkranz, 2018; Cimini et al ACP 2018';
   
   % T-EXPONENT OF SELF-BROADENING [unitless]
   AMU.n_w.value = XS;                     % this is from Ros 2018 (see h2o_list_r18)
   AMU.n_w.uncer = repmat(unknown,1,NL);
   AMU.n_w.uncer(1:2) = [0.5 0.08];        % this is from Cimini et al ACP 2018, Table 1
   AMU.n_w.units = 'unitless';
   AMU.n_w.refer = 'Tretyakov, JMS, 2016'; % I keep the original references
  
   % Line position 
   % Line center linear with pressure frequency shifting coefficients (air)
   % AIR-BROADENED SHIFT PARAMETERS [MHz/mb -> GHz/mb]
   AMU.delta_a.value = SH * 1e3;                 % this is from Ros 2018 (see h2o_list_r18) [GHz/mb -> MHz/mb]
   AMU.delta_a.uncer = repmat(unknown,1,NL); 
   AMU.delta_a.uncer(1:2) = [0.0038 0.0075]; % values from Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD), [MHz/Torr -> MHz/mb]
   AMU.delta_a.units = 'MHz/mb';
   AMU.delta_a.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)';

   % T-EXPONENT OF AIR-SHIFTING [unitless]
   AMU.n_da.value = XH;                      % this is from Ros 2018 (see h2o_list_r18)
   AMU.n_da.uncer = repmat(unknown,1,NL);
   AMU.n_da.units = 'unitless';
   AMU.n_da.refer = 'Rosenkranz 2018 and references therein';
   
   % Line center linear with pressure frequency shifting coefficients (water)
   % SELF-BROADENED SHIFT PARAMETERS [MHz/mb -> GHz/mb]
   AMU.delta_w.value = SHS * 1e3;                 % this is from Ros 2018 (see h2o_list_r18) [GHz/mb -> MHz/mb]
   AMU.delta_w.uncer = repmat(unknown,1,NL); 
   AMU.delta_w.uncer(1:2) = [0.0083 0.0225]; % values from Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD), [MHz/Torr -> MHz/mb]
   AMU.delta_w.units = 'MHz/mb';
   AMU.delta_w.refer = 'Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)';

   % T-EXPONENT OF AIR-SHIFTING [unitless]
   AMU.n_dw.value = XHS;                      % this is from Ros 2018 (see h2o_list_r18)
   AMU.n_dw.uncer = repmat(unknown,1,NL);
   AMU.n_dw.units = 'unitless';
   AMU.n_dw.refer = 'Rosenkranz 2018 and references therein';
    
end % end of if mdl == 'r18'

% Stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% O2 Line parameters 
% Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load uncertainty estimaes from Phil
%UCM = load('/Users/Nico/PROGETTI/IMAA/GAIACLIM/MFILES/FM_Sensitivity/UCM_O2.mat');               
%load('UCM_O2.mat','u');
load('/Users/Nico/PROGETTI/IMAA/GAIACLIM/MFILES/FM_Sensitivity/UCM_O2.mat','u');
u.uSpc   = u.all(1);      % percentage
u.uX     = u.all(2);      % unitless
u.uWB300 = u.all(3);      % GHz/bar
u.uW300  = u.all(4:37);   % GHz/bar
u.uY300  = u.all(38:71);  % 1/bar
u.uV     = u.all(72:105); % 1/bar

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
% Nico: I take values from Ros 2016 with 0.25% uncertainty
AMU.O2BE.value = [ .010,  .014,  .083,  .083,  .207,  .207,  .387,  .387,  .621,  .621,...
                   .910,  .910, 1.255, 1.255, 1.654, 1.654, 2.109, 2.109, 2.618, 2.618,...
                  3.182, 3.182, 3.800, 3.800, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 6.819, 6.819,... 
                  7.709, 7.709, 8.653, 8.653, 9.651, 9.651,...
                   .019,  .048,  .045,  .044,  .049,  .084,  .145,  .136,  .141,  .145,  .201];
AMU.O2BE.uncer = 0.25 * AMU.O2BE.value / 100;  % 0.25%
AMU.O2BE.units = 'unitless'; 
AMU.O2BE.refer = 'Tretyakov, Personal communication, 2016';

% Line shape (i.e. pressure broadening parameters)
% See Table 5 on Tretyakov et al., 2005 - a3 is (nearly) the same as W300 in Rosenk 2016 [GHz/(1e5 Pa)]
% Only the first line differs, plus Ros has more lines at higher freq 
% Table 5 does not report uncertainties; thus I use values on Table 1 and converts units
% The uncertainties from table 1 have to be combined.
% Note that Table 1 allows computation of combined uncertanity only for the first 20 (stronger) lines
% In fact, only the first 20 lines have been investigated independently.
% Since each line was measured independently, their uncertainties can be assumed as uncorrelated
% Following Phil's suggestion, the combined uncertainty is obtained as sqrt(0.21*uO2^2 + 0.79*uN2^2).
% These do not correspond to Phil's results listed in table 1 in his memo of July 11 (O2_off-diagonal.docx), because there they were computed as (0.21*uO2 + 0.79*uN2), which is not correct.
% Units in Ros. model are [MHz/mb]; the conversion factor is 1/torr2mb, i.e. mb2torr
% For other lines I here assume 0 uncertainty, as this is correlated and it's treated as O2gammaWL (weaker lines) below
AMU.O2gamma.value = [1.688, 1.703, 1.513, 1.491, 1.415, 1.408, 1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217, 1.189, 1.174, 1.134, 1.134, 1.089, 1.088, 1.037, 1.038, 0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858, 0.811, 0.811, 0.764, 0.764, 0.717, 0.717, 0.669, 0.669, 2.78, 1.64, 1.64, 1.64, 1.60, 1.60, 1.60, 1.60, 1.62, 1.47, 1.47]; % From Ros 2016
AMU.O2gamma.uncer = zeros(size(AMU.O2gamma.value));                         % start with zero uncer
%c4unc = 1e-3*[10 17 10 10 10 10 8 12 15 10 10 10 10 10 10 15 10 10 15 15];  % O2 - unc of col 4 in Table 1 (MHz/Torr)
%c5unc = 1e-3*[20 30 25 15 15 12 10 10 13 30 10 10 17 15 10 10 18 10 40 15]; % N2 - unc of col 5 in Table 1 (MHz/Torr)
%unc_comb = (0.21*c4unc + 0.79*c5unc);                                       % combined (air) uncertainty [MHz/Torr]
%AMU.O2gamma.uncer(1:20) = unc_comb;                                         % fill first 20 with combined uncert. the rest remains to 0
%AMU.O2gamma.uncer = AMU.O2gamma.uncer * mb2torr;                            % MHz/Torr -> MHz/mb == GHz/bar (this is the same as in table 1 of O2_off-diagonal.docx)
%MTX = load('/Users/Nico/PROGETTI/IMAA/GAIACLIM/MFILES/FM_Sensitivity/DATA/FromPhil/sigma_widths_revised.txt');
%AMU.O2gamma.uncer(1:20) = MTX(1:20,3); % GHz/bar == MHz/mb (no need to change units!) only first 20 are uncorrelated; the other are treated below (O2gamma_WL and O2gamma_mmW)
%AMU.O2gamma.uncer(1:38) = MTX(1:38,3); % GHz/bar == MHz/mb (no need to change units!) per Phil suggestion (2017/10/03), all lines can be perturbed independently (the mmW lines are treated below O2gamma_mmW)
AMU.O2gamma.uncer(1:34) = u.uW300; % lines with N>33 are neglected (GHz/bar == MHz/mb)
AMU.O2gamma.units = 'MHz/mb';
AMU.O2gamma.refer = 'Tretyakov, JMS, 2005 + Rosenkranz Pers. Comm. 2017';
% From M. Tretyakov: Additionally you may use Koshelev et al., JQSRT 2016 (Table 1, gamma_self)
% Nico: measurement by Koshelev have smaller uncertainty, so are covered

% Line shape (i.e. pressure broadening parameters) for Weak Lines (WL)
% From Phil (2017/09/22):
% The width values for the weaker lines (21-38) in Tretyakov et al. (a3 in table 5) are extrapolated linearly
% with quantum number N from the combined O2- and N2-broadening measured values of the stronger lines. 
% That extrapolation introduces correlated uncertainties.
% Per Phil suggestion (2017/07/17-18) I test 0.05 GHz/bar uncertainty for WL (21-38) and leave 0 uncertainty to the other
% This is an arbitrary - likely overestimated - value: 0.05 GHz/bar = 0.05 MHz/mb
% Phil recomputed the uncertainty considering linear regression (2017/10/02)
% Per Phil suggestion (2017/10/03), after his new uncertainty calculations, these lines can be perturbed independently, so are treated above (O2gamma)
% O2gamma_WL should only be used if one wants to compute the impact of uncertainties from all weaker lines together
MTX = load('/Users/Nico/PROGETTI/IMAA/GAIACLIM/MFILES/FM_Sensitivity/DATA/FromPhil/sigma_widths_revised.txt');
AMU.O2gamma_WL.value = AMU.O2gamma.value;
AMU.O2gamma_WL.uncer = zeros(size(AMU.O2gamma.value)); 
%AMU.O2gamma_WL.uncer([21:38]) = 0.05; %  arbitrary value 0.05 GHz/bar = 0.05 MHz/mb (no need to change units!)
AMU.O2gamma_WL.uncer([21:38]) = MTX(21:38,3); % GHz/bar == MHz/mb (no need to change units!) - these are totally correlated, uncertainty of the linear regression
AMU.O2gamma_WL.units = 'MHz/mb';
AMU.O2gamma_WL.refer = 'Rosenkranz, pers. comm., 2017';

% Line shape (i.e. pressure broadening parameters) for millimeter-wave (mmW) lines
% This is to check Phil's assumption (2017/09/27):
% This assumes that the millimeter-wave lines don't have an impact at 20-60 GHz.
% I test 0.05 GHz/bar uncertainty for mmW (39-49) and leave 0 uncertainty to the other
% This is an arbitrary value: 0.05 GHz/bar = 0.05 MHz/mb
AMU.O2gamma_mmW.value = AMU.O2gamma.value;
AMU.O2gamma_mmW.uncer = zeros(size(AMU.O2gamma.value)); 
%AMU.O2gamma_mmW.uncer([39:49]) = 0.05; %  arbitrary value 0.05 GHz/bar = 0.05 MHz/mb (no need to change units!)
AMU.O2gamma_mmW.uncer([39:49]) = 0.1 * AMU.O2gamma_mmW.value([39:49]); %  arbitrary value 10% uncertainty (this should really be overestimated, 3-5 times larger than 0.05 MHz/mb)
AMU.O2gamma_mmW.units = 'MHz/mb';
AMU.O2gamma_mmW.refer = 'Rosenkranz, pers. comm., 2017';

% Line shape (i.e. pressure broadening parameters) for neglected lines (NL)
% This is to show that lines with quantum number higher than 33+ (i.e. i>34) make negligible impact at 20-60 GHz.
% I use 10% uncertainty (which should really be overestimated)
AMU.O2gamma_NL.value = AMU.O2gamma.value;
AMU.O2gamma_NL.uncer = zeros(size(AMU.O2gamma.value)); 
AMU.O2gamma_NL.uncer([35:49]) = 0.1 * AMU.O2gamma_NL.value([35:49]); %  arbitrary value 10% uncertainty (this should really be overestimated, 3-5 times larger than 0.05 MHz/mb)
AMU.O2gamma_NL.units = 'MHz/mb';
AMU.O2gamma_NL.refer = 'Rosenkranz, pers. comm., 2017';

% Line intensity of non-resonant pseudo-line (see Janssen TABLE 2A.1 caption, page 84)
% This is hard-coded in MPM (so I had to change it - I called it Snr)
% I assume 5% following suggestions from Misha (5-10%) and HITRAN
% uncertainty code (5-> 5-10%). More conservative would have been 10%, but
% since we assume 15% for the broadening, Phil says it's OK, as these are highly negatively correlated 
AMU.Snr.value = 1.584e-17;  % 1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
AMU.Snr.uncer = AMU.Snr.value * 5/100; % 5%
AMU.Snr.units = 'Hz*cm2/GHz2'; % it has not the same units as other lines because of different definition (to avoid 1/frq -> inf for frq->0)
AMU.Snr.refer = 'Pickett et al., 1998, i.e. JPL line compilation';

% Pressure broadening of non-resonant pseudo-line (see Janssen TABLE 2A.1 caption, page 84)
% This is called WB300 in MPM
% Units in Ros. model are [MHz/mb]; so no need of the conversion factor mb2torr here
AMU.WB300.value = 0.56;  % MHz/mb
%AMU.WB300.uncer = 0.084; % GHz/bar (last entry of table 1 in his memo O2_off-diagonal.docx) == MHz/mb
AMU.WB300.uncer = 0.05; % GHz/bar (last entry of table 1 in his memo O2_off-diagonal.docx) == MHz/mb
AMU.WB300.units = 'MHz/mb';
AMU.WB300.refer = 'Rosenkranz, pers. comm., 2017';

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

% Line mixing coefficients (single lines)
% Table 5, last 2 columns of Tretyakov 2005 (i.e. a5+a6). This is called Y300 in MPM
% Uncertainty was computed by Phil Rosenkranz, through full covariance (see O2_off-diagonal_NC.docx)
% Here I take the values in u.all(38:71) after >load UCM_O2.mat
AMU.Y300.value = [-0.0360, 0.2547, -0.3655,  0.5495,...
                  -0.5696,  0.6181, -0.4252,  0.3517, -0.1496,  0.0430,...
                   0.0640, -0.1605,  0.2906, -0.3730,  0.4169, -0.4819,...
                   0.4963, -0.5481,  0.5512, -0.5931,  0.6212, -0.6558,...
                   0.6920, -0.7208,  0.7312, -0.7550,  0.7555, -0.7751,...
                   0.7914, -0.8073,  0.8307, -0.8431,  0.8676, -0.8761,...
                   0.9046, -0.9092,  0.9416, -0.9423,...  
                   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.];
AMU.Y300.uncer = zeros(size(AMU.Y300.value));           
AMU.Y300.uncer(1:34) = u.uY300;  % lines with N>33 are neglected 
AMU.Y300.units = '1/bar == 1/1e5Pa =~ 1/atm'; % 1/1e5Pa == 1/1e3hPa == 1/1e3mb == 1/bar =~ 1/atm
AMU.Y300.refer = 'Tretyakov et al., 2005; Uncertainty from Rosenkranz, pers. comm., 2017';


% Line mixing coefficients for neglected lines (NL)
% This is to show that lines with quantum number higher than 33+ (i.e. i>34) make negligible impact at 20-60 GHz.
% I use 20% uncertainty (which should really be overestimated)
AMU.Y300_NL.value = AMU.Y300.value;
AMU.Y300_NL.uncer = zeros(size(AMU.Y300.value)); 
AMU.Y300_NL.uncer([35:49]) = 0.2 * AMU.Y300.value([35:49]); %  arbitrary value 20% uncertainty (this should be overestimated)
AMU.Y300_NL.units = '1/mb';
AMU.Y300_NL.refer = 'Rosenkranz, pers. comm., 2017';


% O2 line mixing temperature dependence
% Table 5, last column of Tretyakov 2005. This is called V in MPM
AMU.O2_V.value = [  0.0079, -0.0978,  0.0844, -0.1273,...
                    0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,...
                    0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,...
                    0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,...
                    0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,...
                    0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545,...
                    0.680,  -0.660,   0.685,  -0.665,...
                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.];
AMU.O2_V.uncer = zeros(size(AMU.O2_V.value));           
%AMU.O2_V.uncer = AMU.O2_V.value * 0.2; % Phil suggested "perturb all line v-param 1-38 by 20%, simultaneously", see email on 2017/07/18
%AMU.O2_V.uncer(39:end) = 0.0;          % Phil suggested "perturb all line v-param 1-38 by 20%, simultaneously", see email on 2017/07/18
AMU.O2_V.uncer(1:34) = u.uV;  % lines with N>33 are neglected;
AMU.O2_V.units = '1/bar == 1/1e5Pa =~ 1/atm'; % 1/1e5Pa == 1/1e3hPa == 1/1e3mb == 1/bar =~ 1/atm
AMU.O2_V.refer = 'Tretyakov et al., 2005. Uncertainty from Rosenkranz, pers. comm., 2017';

% O2 line mixing temperature dependence for neglected lines (NL)
% This is to show that lines with quantum number higher than 33+ (i.e. i>34) make negligible impact at 20-60 GHz.
% I use 20% uncertainty (which should really be overestimated)
AMU.O2_V_NL.value = AMU.O2_V.value;
AMU.O2_V_NL.uncer = zeros(size(AMU.O2gamma.value)); 
AMU.O2_V_NL.uncer([35:49]) = 0.2 * AMU.O2gamma.value([35:49]); %  arbitrary value 20% uncertainty (this should be overestimated)
AMU.O2_V_NL.units = '1/mb';
AMU.O2_V_NL.refer = 'Rosenkranz, pers. comm., 2017';

% Intensity temperature-dependence exponent
% See Partition_sum.doc and related discussions
% From Phil: I checked the power-law approximation for Q(T) against the exact computation for oxygen and got 0.1% agreement (recall my e-mail of Oct. 6).
AMU.O2_nS.value = 2.0;
AMU.O2_nS.uncer = AMU.O2_nS.value * 0.001; % value computed by Rosenkranz (email of 2017/10/06 and Partition_sum.doc)
AMU.O2_nS.units = 'unitless';
AMU.O2_nS.refer = 'Rosenkranz, email of 2017/10/06, based on Gamache et al., 2017';

% Stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only 183 GHz SD parameters (Koshelev et al., 2020)
if strcmp(mdl,'r20sd') 
      
   % start with a clean AMU 
   % it may have to be removed if parameters other than SD are to be considered
   clear AMU
   
   % Read the list of parameters
   h2o_sdlist_r20a % this is the same as h2o_sdlist_r19 but the two coefficients W2air W2self at 22.2 GHz (which were missing in h2o_sdlist_r19)

   % SD line parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % given in MHz/Torr in Koshelev et al, 2018 & 2020 while units in Ros. model are [MHz/mb]; 
   % the conversion factor is 1/torr2mb, i.e. mb2torr
   
   % BROADENING PARAMETERS 
   
   % gamma0_a: air-broadening width
   AMU.gamma0_a.value = [3.653 4.044 repmat(unknown,1,14)];
   AMU.gamma0_a.uncer = [0.022 0.007 repmat(unknown,1,14)];
   AMU.gamma0_a.units = 'MHz/Torr';
   AMU.gamma0_a.refer = 'Koshelev et al, 2018; 2020';
   AMU.gamma0_a.value = AMU.gamma0_a.value*mb2torr;  % [MHz/Torr] -> [MHz/mb]
   AMU.gamma0_a.uncer = AMU.gamma0_a.uncer*mb2torr; % [MHz/Torr] -> [MHz/mb]
   AMU.gamma0_a.units = 'MHz/mb';

   % ng0_a: T-exponent of air-broadening width [unitless]
   AMU.ng0_a.value = [unknown 0.617 repmat(unknown,1,14)];
   AMU.ng0_a.uncer = [unknown 0.011 repmat(unknown,1,14)];
   AMU.ng0_a.units = 'unitless';
   AMU.ng0_a.refer = 'Koshelev et al, 2020';

   % gamma0_s: self-broadening width
   AMU.gamma0_w.value = [18.17 20.013 repmat(unknown,1,14)];
   AMU.gamma0_w.uncer = [0.03 0.019 repmat(unknown,1,14)];
   AMU.gamma0_w.units = 'MHz/Torr';
   AMU.gamma0_w.refer = 'Koshelev et al, 2018; 2020';
   AMU.gamma0_w.value = AMU.gamma0_w.value*mb2torr;  % [MHz/Torr] -> [MHz/mb]
   AMU.gamma0_w.uncer = AMU.gamma0_w.uncer*mb2torr; % [MHz/Torr] -> [MHz/mb]
   AMU.gamma0_w.units = 'MHz/mb';

   % ng0_s: T-exponent of self-broadening width [unitless]
   AMU.ng0_w.value = [unknown 0.821 repmat(unknown,1,14)];
   AMU.ng0_w.uncer = [unknown 0.008 repmat(unknown,1,14)];
   AMU.ng0_w.units = 'unitless';
   AMU.ng0_w.refer = 'Koshelev et al, 2020';
   
   % gamma2_a: SD air-broadening width
   AMU.gamma2_a.value = [0.580 0.543 repmat(unknown,1,14)];
   AMU.gamma2_a.uncer = [0.019 0.004 repmat(unknown,1,14)];
   AMU.gamma2_a.units = 'MHz/Torr';
   AMU.gamma2_a.refer = 'Koshelev et al, 2018; 2020';
   AMU.gamma2_a.value = AMU.gamma2_a.value*mb2torr;  % [MHz/Torr] -> [MHz/mb]
   AMU.gamma2_a.uncer = AMU.gamma2_a.uncer*mb2torr; % [MHz/Torr] -> [MHz/mb]
   AMU.gamma2_a.units = 'MHz/mb';

   % ng2_a: T-exponent of SD air-broadening width [unitless]
   AMU.ng2_a.value = [unknown 0.412 repmat(unknown,1,14)];
   AMU.ng2_a.uncer = [unknown 0.055 repmat(unknown,1,14)];
   AMU.ng2_a.units = 'unitless';
   AMU.ng2_a.refer = 'Koshelev et al, 2020';

   % gamma2_s: SD self-broadening width
   AMU.gamma2_w.value = [2.55 1.946 repmat(unknown,1,14)];
   AMU.gamma2_w.uncer = [0.06 0.023 repmat(unknown,1,14)];
   AMU.gamma2_w.units = 'MHz/Torr';
   AMU.gamma2_w.refer = 'Koshelev et al, 2018; 2020';
   AMU.gamma2_w.value = AMU.gamma2_w.value*mb2torr;  % [MHz/Torr] -> [MHz/mb]
   AMU.gamma2_w.uncer = AMU.gamma2_w.uncer*mb2torr; % [MHz/Torr] -> [MHz/mb]
   AMU.gamma2_w.units = 'MHz/mb';

   % ng2_s: T-exponent of SD self-broadening width [unitless]
   AMU.ng2_w.value = [unknown 0.571 repmat(unknown,1,14)];
   AMU.ng2_w.uncer = [unknown 0.095 repmat(unknown,1,14)];
   AMU.ng2_w.units = 'unitless';
   AMU.ng2_w.refer = 'Koshelev et al, 2020';

   % SHIFTING PARAMETERS - given in MHz/Torr in Koshelev et al, 2018 & 2020
   
   % delta0_a: air-broadening shift
   AMU.delta0_a.value = [-0.045 -0.0990 repmat(unknown,1,14)];
   AMU.delta0_a.uncer = [0.006 0.0020 repmat(unknown,1,14)];
   AMU.delta0_a.units = 'MHz/Torr';
   AMU.delta0_a.refer = 'Koshelev et al, 2018; 2020';
   AMU.delta0_a.value = AMU.delta0_a.value*mb2torr;  % [MHz/Torr] -> [MHz/mb]
   AMU.delta0_a.uncer = AMU.delta0_a.uncer*mb2torr; % [MHz/Torr] -> [MHz/mb]
   AMU.delta0_a.units = 'MHz/mb';

   % nd0_a: T-exponent of air-broadening shift [unitless]
   AMU.nd0_a.value = [unknown 1.82 repmat(unknown,1,14)];
   AMU.nd0_a.uncer = [unknown 0.12 repmat(unknown,1,14)];
   AMU.nd0_a.units = 'unitless';
   AMU.nd0_a.refer = 'Koshelev et al, 2020';
   
   % delta0_s: self-broadening shift
   AMU.delta0_w.value = [0.968 0.181 repmat(unknown,1,14)];
   AMU.delta0_w.uncer = [0.010 0.016 repmat(unknown,1,14)];
   AMU.delta0_w.units = 'MHz/Torr';
   AMU.delta0_w.refer = 'Koshelev et al, 2018; 2020';
   AMU.delta0_w.value = AMU.delta0_w.value*mb2torr;  % [MHz/Torr] -> [MHz/mb]
   AMU.delta0_w.uncer = AMU.delta0_w.uncer*mb2torr; % [MHz/Torr] -> [MHz/mb]
   AMU.delta0_w.units = 'MHz/mb';

   % nd0_s: T-exponent of self-broadening shift [unitless]
   AMU.nd0_w.value = [unknown 0.98 repmat(unknown,1,14)];
   AMU.nd0_w.uncer = [unknown 0.44 repmat(unknown,1,14)];
   AMU.nd0_w.units = 'unitless';
   AMU.nd0_w.refer = 'Koshelev et al, 2020';

   % A_a: Frost parameter air-broadening shift [unitless]
   AMU.A_a.value = [unknown 0.0 repmat(unknown,1,14)];
   AMU.A_a.uncer = [unknown 0.1 repmat(unknown,1,14)]; % assuming 0.1 uncertainty (it should be 0, as it's actually kept fixed, but I need it for computing Kb as Aair is correlated to delta0_a and nd0_a)
   AMU.A_a.units = 'unitless';
   AMU.A_a.refer = 'Koshelev et al, 2020';
   
   % A_s: Frost parameter self-broadening shift [unitless]
   AMU.A_w.value = [unknown 12.57 repmat(unknown,1,14)];
   AMU.A_w.uncer = [unknown 1.21 repmat(unknown,1,14)];
   AMU.A_w.units = 'unitless';
   AMU.A_w.refer = 'Koshelev et al, 2020';
   
    % delta2_a: SD air-broadening shift
   AMU.delta2_a.value = [-0.002 -0.021 repmat(unknown,1,14)];
   AMU.delta2_a.uncer = [0.008 0.038 repmat(unknown,1,14)];
   AMU.delta2_a.units = 'MHz/Torr';
   AMU.delta2_a.refer = 'Koshelev et al, 2018; 2020';
   AMU.delta2_a.value = AMU.delta2_a.value*mb2torr;  % [MHz/Torr] -> [MHz/mb]
   AMU.delta2_a.uncer = AMU.delta2_a.uncer*mb2torr; % [MHz/Torr] -> [MHz/mb]
   AMU.delta2_a.units = 'MHz/mb';

   % delta2_s: SD self-broadening shift
   AMU.delta2_w.value = [0.060 0.21 repmat(unknown,1,14)];
   AMU.delta2_w.uncer = [0.060 0.11 repmat(unknown,1,14)];
   AMU.delta2_w.units = 'MHz/Torr';
   AMU.delta2_w.refer = 'Koshelev et al, 2018; 2020';
   AMU.delta2_w.value = AMU.delta2_w.value*mb2torr;  % [MHz/Torr] -> [MHz/mb]
   AMU.delta2_w.uncer = AMU.delta2_w.uncer*mb2torr; % [MHz/Torr] -> [MHz/mb]
   AMU.delta2_w.units = 'MHz/mb';  

   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Other parameters are as 'r18' above
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   % Continuum parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   AMU.con_Cf.value = 5.96e-10; % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Cf.value = CF;       % this is from latest h2o_sdlist_r20a
   AMU.con_Cf.uncer = 5.50e-11; % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Cf.units = '1/(km*(mb^2*GHz^2))';
   AMU.con_Cf.refer = 'Rosenkranz, 2018; Cimini et al ACP 2018';

   AMU.con_Cs.value = 1.42e-8;  % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Cs.value = CS;       % this is from latest h2o_sdlist_r20a
   AMU.con_Cs.uncer = 3.20e-9;  % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Cs.units = '1/(km*(mb^2*GHz^2))';
   AMU.con_Cs.refer = 'Rosenkranz, 2018; Cimini et al ACP 2018';
 
   AMU.con_Xf.value = XCF;      % this is from latest h2o_sdlist_r20a
   AMU.con_Xf.uncer = 0.8;      % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Xf.units = 'unitless';
   AMU.con_Xf.refer = 'Tretyakov, JMS, 2016; Koshelev et al. 2011'; % I keep the original references
   
   AMU.con_Xs.value = XCS;      % this is from latest h2o_sdlist_r20a
   AMU.con_Xs.uncer = 0.6;      % this is from Cimini et al ACP 2018, Table 1
   AMU.con_Xs.units = 'unitless';
   AMU.con_Xs.refer = 'Tretyakov, JMS, 2016'; % I keep the original references

   % Since r18 I perturb the continuum parameters directly, not the factor
   AMU.con_Cf_factr.value = 1; % 
   AMU.con_Cf_factr.uncer = 0; % 
   AMU.con_Cf_factr.units = 'unitless';
   AMU.con_Cf_factr.refer = 'With R18 I perturb the continuum parameters, not the factor';

   AMU.con_Cs_factr.value = 1; %
   AMU.con_Cs_factr.uncer = 0; % 
   AMU.con_Cs_factr.units = 'unitless';
   AMU.con_Cs_factr.refer = 'With R18 I perturb the continuum parameters, not the factor';
   
   % Lines' parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   NL = length(FL);
   
   % Line frequency FL (freq,GHz)
   AMU.FL.value = FL;                           % this is from latest h2o_sdlist_r20a
   AMU.FL.uncer = repmat(unknown,1,NL);
   %AMU.FL.uncer(1:2) = [0.00000005   0.000001]; % this is from Cimini et al ACP 2018, Table 1
   AMU.FL.uncer(1:2) = [0.00000005   0.000002]; % this is from Koshelev et al, 2020 (covariance table)
   AMU.FL.units = 'GHz';
   AMU.FL.refer = 'Tretyakov, JMS, 2016; uncert from Koshelev et al, 2020 (covariance table)';
   
   % Line intensity (or strength) S
   AMU.S.value = S1;                     % this is from latest h2o_sdlist_r20a
   AMU.S.uncer = repmat(unknown,1,NL);
   AMU.S.uncer(1:2) = S1(1:2) / 100;     % 1% is from Cimini et al ACP 2018, Table 1
   AMU.S.units = 'Hz*cm^2';
   AMU.S.refer = 'Tretyakov, JMS, 2016'; % I keep the original references
   
   % T COEFF. OF INTENSITIES (B2 in Ros model)
   AMU.B2.value = B2;                     % this is from latest h2o_sdlist_r20a
   AMU.B2.uncer = B2 / 100 * 0.5;         % 0.5% is from Cimini et al ACP 2018, Table 1
   AMU.B2.units = 'Unitless';
   AMU.B2.refer = 'Rosenkranz, 2016 + Tretyakov pers. comm.'; % I keep the original references
    
   % Intensity temperature-dependence exponent (see Partition_sum.doc)
   AMU.wv_nS.value = 2.5;
   %AMU.wv_nS.uncer = AMU.wv_nS.value * 0.06; % value computed by Tretyakov (Partition_sum.doc)
   AMU.wv_nS.uncer = AMU.wv_nS.value * 0.005; % value computed by Rosenkranz (email of 2018/03/08 and h2o_partition_sum.txt)
   AMU.wv_nS.units = 'unitless';
   AMU.wv_nS.refer = 'Rosenkranz, email of 2018/03/08, based on Gamache et al., 2017';
    
end % end if strcmp(mdl,'r20sd')

return
