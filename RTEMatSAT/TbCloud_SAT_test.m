addpath(['/Users/Nico/Matools/FromDCT/']);
addpath(['/Users/Nico/Matools/RTEMat/']);

% Forward model absorption settings
% Model, frequency and angle settings
AM = 'ros19sd';
ang = 90.;
frq = [20:1:200]; % low res
nf = length(frq);
Es = 1.0;     % from sat surface emissivity is needed (0.6 was used by Payne et al, 2011 (Figure 11))
RT = 1; % use ray tracing

% Model atmosphere
ATM(1).name = 'Tropical';           ATM(1).clr = 'r';
ATM(2).name = 'Midlatitude Summer'; ATM(2).clr = 'm';
ATM(3).name = 'Midlatitude Winter'; ATM(3).clr = 'g';
ATM(4).name = 'Subarctic Summer';   ATM(4).clr = 'b';
ATM(5).name = 'Subarctic Winter';   ATM(5).clr = 'c';
ATM(6).name = 'U.S. Standard';      ATM(6).clr = 'k';

figure

% Processing & printout
for ipro = 1:6
%for ipro = 6 % US standard
%for ipro = 1 % Tropical

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loading Atmosphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Load standard atmosphere (low res at lower levels, only 1 level within 1 km)
    [z,p,t,d,md,gasids] = glatm(ipro);
    % Cutting at 40 km
    %indx = 1:32; z = z(indx); p = p(indx); t = t(indx); md = md(indx,:);
    gkg = ppmv2gkg(md(:,1),gasids(1));
    rh = mr2rh(p,t,gkg)/100;

    % Making higher resolution
    %zlr = z;
    %z = ([0:10:100 150:50:500 600:100:1000 1500:500:5000 6000:1000:10000 12000:2000:40000]/1000)';
    %p = interp1(zlr,p,z);
    %t = interp1(zlr,t,z);
    %rh = interp1(zlr,rh,z);
    
    % default no cloud
    denliq = zeros(size(z));
    denice = zeros(size(z));
    cldh = zeros(2,0);

    O1 = TbCloud_SAT(z,p,t,rh,denliq,denice,cldh,frq,ang,AM,RT,Es);
    
    plot(frq,O1.tbtotal,ATM(ipro).clr); hold on;
    addlegend(ATM(ipro).name);
    
end


xlabel('Frequency (GHz)'); ylabel('T_B(K)'); grid on;
%legend(ATM(1).name,ATM(2).name,ATM(3).name,ATM(4).name,ATM(5).name,ATM(6).name,'Location','Best');
format4paper(gcf);
%saveas(gcf,['/Users/Nico/PROGETTI/IMAA/GAIACLIM/MFILES/FM_Sensitivity/SAT/UpwellingTb_all_e' num2str(Es) '_20_200.png'],'png');

