function [rho,a,b,r0,r] = polyTEOS10_bsq(SA,CT,p)
% polyTEOS10_bsq              in-situ density (55-term polynomial equation)
%==========================================================================
%
% USAGE:
%     [rho,a,b,r0,r] = polyTEOS10_bsq(SA,CT,p)
%
% DESCRIPTION:
%  Calculates in-situ density from Absolute Salinity, Conservative
%  Temperature and pressure, using the computationally-efficient 55-term
%  polynomial expression for density (Roquet et al., 2014).
%
%  Note that the 55-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is available to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT & p need to have the same dimensions.
%
% OUTPUT:
%  rho   =  in situ density                                  [ kg/m^3 ]
%  a     =  Boussinesq thermal expansion (=-dr/dCT)        [ kg/m^3/K ]
%  b     =  Boussinesq haline contraction (=dr/dSA)   [ kg/m^3/(g/kg) ]
%  r0    =  vertical reference density                       [ kg/m^3 ]
%  r     =  density anomaly (rho = r0 + r)                   [ kg/m^3 ]
%
% CHECK VALUES (SA=30g/kg, CT=10degC, p=1e3dbar):
%  rho = 1027.45140
%  a   = 0.179646281
%  b   = 0.765555368
%  r0  = 4.59763035
%  r   = 1022.85377
%
% REFERENCES:
%  Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate 
%   polynomial expressions for the density and specific volume of 
%   seawater using the TEOS-10 standard. Ocean Modelling.
%  McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  Journal of Atmospheric and 
%   Oceanic Technology, 20, 730-741. 
%
% AUTHOR:
%  Fabien Roquet
%  jan 2015

% reduced variables
SAu = 40*35.16504/35; CTu = 40; Zu=1e4; deltaS = 32;
ss   = sqrt ( (SA+deltaS)/SAu );
tt   = CT / CTu;
pp   =   p / Zu;

% vertical reference profile of density
R00 = 4.6494977072e+01; R01 = -5.2099962525e+00; R02 = 2.2601900708e-01;
R03 = 6.4326772569e-02; R04 = 1.5616995503e-02; R05 = -1.7243708991e-03;
r0 = ( ( ( ( ( R05*pp + R04 ).*pp + R03 ).*pp + R02 ).*pp + R01 ).*pp + R00 ).*pp ;

% density anomaly
R000 = 8.0189615746e+02; R100 = 8.6672408165e+02; R200 = -1.7864682637e+03;
R300 = 2.0375295546e+03; R400 = -1.2849161071e+03; R500 = 4.3227585684e+02;
R600 = -6.0579916612e+01; R010 = 2.6010145068e+01; R110 = -6.5281885265e+01;
R210 = 8.1770425108e+01; R310 = -5.6888046321e+01; R410 = 1.7681814114e+01;
R510 = -1.9193502195e+00; R020 = -3.7074170417e+01; R120 = 6.1548258127e+01;
R220 = -6.0362551501e+01; R320 = 2.9130021253e+01; R420 = -5.4723692739e+00;
R030 = 2.1661789529e+01; R130 = -3.3449108469e+01; R230 = 1.9717078466e+01;
R330 = -3.1742946532e+00; R040 = -8.3627885467e+00; R140 = 1.1311538584e+01;
R240 = -5.3563304045e+00; R050 = 5.4048723791e-01; R150 = 4.8169980163e-01;
R060 = -1.9083568888e-01; R001 = 1.9681925209e+01; R101 = -4.2549998214e+01;
R201 = 5.0774768218e+01; R301 = -3.0938076334e+01; R401 = 6.6051753097e+00;
R011 = -1.3336301113e+01; R111 = -4.4870114575e+00; R211 = 5.0042598061e+00;
R311 = -6.5399043664e-01; R021 = 6.7080479603e+00; R121 = 3.5063081279e+00;
R221 = -1.8795372996e+00; R031 = -2.4649669534e+00; R131 = -5.5077101279e-01;
R041 = 5.5927935970e-01; R002 = 2.0660924175e+00; R102 = -4.9527603989e+00;
R202 = 2.5019633244e+00; R012 = 2.0564311499e+00; R112 = -2.1311365518e-01;
R022 = -1.2419983026e+00; R003 = -2.3342758797e-02; R103 = -1.8507636718e-02;
R013 = 3.7969820455e-01;
rz3 =       R013.*tt + R103.*ss + R003;
rz2 =      (R022.*tt+R112.*ss+R012).*tt+(R202.*ss+R102).*ss+R002;
rz1 =    (((R041.*tt+R131.*ss+R031).*tt ...
    +   (R221.*ss+R121).*ss+R021).*tt ...
    +  ((R311.*ss+R211).*ss+R111).*ss+R011).*tt ...
    + (((R401.*ss+R301).*ss+R201).*ss+R101).*ss+R001;
rz0 =  (((((R060.*tt+R150.*ss+R050).*tt ...
    +    (R240.*ss+R140).*ss+R040).*tt ...
    +   ((R330.*ss+R230).*ss+R130).*ss+R030).*tt ...
    +  (((R420.*ss+R320).*ss+R220).*ss+R120).*ss+R020).*tt ...
    + ((((R510.*ss+R410).*ss+R310).*ss+R210).*ss+R110).*ss+R010).*tt ...
    +(((((R600.*ss+R500).*ss+R400).*ss+R300).*ss+R200).*ss+R100).*ss+R000;
r = ( ( rz3.*pp + rz2 ).*pp + rz1 ).*pp + rz0;

% in-situ density
rho = r + r0;

% thermal expansion a
ALP000 = -6.5025362670e-01; ALP100 = 1.6320471316e+00; ALP200 = -2.0442606277e+00;
ALP300 = 1.4222011580e+00; ALP400 = -4.4204535284e-01; ALP500 = 4.7983755487e-02;
ALP010 = 1.8537085209e+00; ALP110 = -3.0774129064e+00; ALP210 = 3.0181275751e+00;
ALP310 = -1.4565010626e+00; ALP410 = 2.7361846370e-01; ALP020 = -1.6246342147e+00;
ALP120 = 2.5086831352e+00; ALP220 = -1.4787808849e+00; ALP320 = 2.3807209899e-01;
ALP030 = 8.3627885467e-01; ALP130 = -1.1311538584e+00; ALP230 = 5.3563304045e-01;
ALP040 = -6.7560904739e-02; ALP140 = -6.0212475204e-02; ALP050 = 2.8625353333e-02;
ALP001 = 3.3340752782e-01; ALP101 = 1.1217528644e-01; ALP201 = -1.2510649515e-01;
ALP301 = 1.6349760916e-02; ALP011 = -3.3540239802e-01; ALP111 = -1.7531540640e-01;
ALP211 = 9.3976864981e-02; ALP021 = 1.8487252150e-01; ALP121 = 4.1307825959e-02;
ALP031 = -5.5927935970e-02; ALP002 = -5.1410778748e-02; ALP102 = 5.3278413794e-03;
ALP012 = 6.2099915132e-02; ALP003 = -9.4924551138e-03;
a = ((ALP003.*pp	...
    + ALP012.*tt + ALP102.*ss+ALP002).*pp	...
    + ((ALP031.*tt	...
    + ALP121.*ss+ALP021).*tt	...
    + (ALP211.*ss+ALP111).*ss+ALP011).*tt	...
    + ((ALP301.*ss+ALP201).*ss+ALP101).*ss+ALP001).*pp	...
    + ((((ALP050.*tt	...
    + ALP140.*ss+ALP040).*tt	...
    + (ALP230.*ss+ALP130).*ss+ALP030).*tt	...
    + ((ALP320.*ss+ALP220).*ss+ALP120).*ss+ALP020).*tt	...
    + (((ALP410.*ss+ALP310).*ss+ALP210).*ss+ALP110).*ss+ALP010).*tt	...
    + ((((ALP500.*ss+ALP400).*ss+ALP300).*ss+ALP200).*ss+ALP100).*ss+ALP000;

% haline contraction b
BET000 = 1.0783203594e+01; BET100 = -4.4452095908e+01; BET200 = 7.6048755820e+01;
BET300 = -6.3944280668e+01; BET400 = 2.6890441098e+01; BET500 = -4.5221697773e+00;
BET010 = -8.1219372432e-01; BET110 = 2.0346663041e+00; BET210 = -2.1232895170e+00;
BET310 = 8.7994140485e-01; BET410 = -1.1939638360e-01; BET020 = 7.6574242289e-01;
BET120 = -1.5019813020e+00; BET220 = 1.0872489522e+00; BET320 = -2.7233429080e-01;
BET030 = -4.1615152308e-01; BET130 = 4.9061350869e-01; BET230 = -1.1847737788e-01;
BET040 = 1.4073062708e-01; BET140 = -1.3327978879e-01; BET050 = 5.9929880134e-03;
BET001 = -5.2937873009e-01; BET101 = 1.2634116779e+00; BET201 = -1.1547328025e+00;
BET301 = 3.2870876279e-01; BET011 = -5.5824407214e-02; BET111 = 1.2451933313e-01;
BET211 = -2.4409539932e-02; BET021 = 4.3623149752e-02; BET121 = -4.6767901790e-02;
BET031 = -6.8523260060e-03; BET002 = -6.1618945251e-02; BET102 = 6.2255521644e-02;
BET012 = -2.6514181169e-03; BET003 = -2.3025968587e-04;
b = ((BET003.*pp	...
    + BET012.*tt + BET102.*ss+BET002).*pp	...
    + ((BET031.*tt	...
    + BET121.*ss+BET021).*tt	...
    + (BET211.*ss+BET111).*ss+BET011).*tt	...
    + ((BET301.*ss+BET201).*ss+BET101).*ss+BET001).*pp	...
    + ((((BET050.*tt	...
    + BET140.*ss+BET040).*tt	...
    + (BET230.*ss+BET130).*ss+BET030).*tt	...
    + ((BET320.*ss+BET220).*ss+BET120).*ss+BET020).*tt	...
    + (((BET410.*ss+BET310).*ss+BET210).*ss+BET110).*ss+BET010).*tt	...
    + ((((BET500.*ss+BET400).*ss+BET300).*ss+BET200).*ss+BET100).*ss+BET000;
b = b ./ ss ;



