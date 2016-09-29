function [specvol,alpha,beta,v0,delta] = polyTEOS10_55t(SA,CT,p)
% polyTEOS10_55t              specific volume (55-term polynomial equation)
%==========================================================================
%
% USAGE:
%     [specvol,alpha,beta,v0,delta] = polyTEOS10_55t(SA,CT,p)
%
% DESCRIPTION:
%  Calculates specific volume from Absolute Salinity, Conservative
%  Temperature and pressure, using the computationally-efficient 55-term
%  polynomial expression for specific volume (Roquet et al., 2014).
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
%  specvol   =  specific volume                                  [ m^3/kg ]
%  alpha     =  Thermal expansion  (=v.dv/dCT)                      [ 1/K ]
%  beta      =  Haline contraction (=-v.dv/dSA)                [ 1/(g/kg) ]
%  v0        =  vertical reference specific volume               [ m^3/kg ]
%  delta     =  specific volume anomaly (specvol = v0 + delta)   [ m^3/kg ]
%
% CHECK VALUES (SA=30g/kg, CT=10degC, p=1e3dbar):
%  specvol = 9.732820466e-04
%  alpha   = 1.748553121e-04
%  beta    = 7.450974025e-04
%  v0      = -4.333016903e-06
%  delta   = 9.776150635e-04
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
SAu = 40*35.16504/35; CTu = 40; Pu=1e4; deltaS = 24;
ss  = sqrt ( (SA+deltaS)/SAu );
tt  = CT / CTu ;
pp  = p / Pu ;

% vertical reference profile of specific volume
V00 = -4.4015007269e-05; V01 = 6.9232335784e-06; V02 = -7.5004675975e-07; 
V03 = 1.7009109288e-08; V04 = -1.6884162004e-08; V05 = 1.9613503930e-09; 
v0  = (((((V05.*pp+V04).*pp+V03 ).*pp+V02 ).*pp+V01).*pp+V00).*pp ;

% specific volume anomaly
V000 = 1.0772899069e-03; V100 = -3.1263658781e-04; V200 = 6.7615860683e-04; 
V300 = -8.6127884515e-04; V400 = 5.9010812596e-04; V500 = -2.1503943538e-04; 
V600 = 3.2678954455e-05; V010 = -1.4949652640e-05; V110 = 3.1866349188e-05; 
V210 = -3.8070687610e-05; V310 = 2.9818473563e-05; V410 = -1.0011321965e-05; 
V510 = 1.0751931163e-06; V020 = 2.7546851539e-05; V120 = -3.6597334199e-05; 
V220 = 3.4489154625e-05; V320 = -1.7663254122e-05; V420 = 3.5965131935e-06; 
V030 = -1.6506828994e-05; V130 = 2.4412359055e-05; V230 = -1.4606740723e-05; 
V330 = 2.3293406656e-06; V040 = 6.7896174634e-06; V140 = -8.7951832993e-06; 
V240 = 4.4249040774e-06; V050 = -7.2535743349e-07; V150 = -3.4680559205e-07; 
V060 = 1.9041365570e-07; V001 = -1.6889436589e-05; V101 = 2.1106556158e-05; 
V201 = -2.1322804368e-05; V301 = 1.7347655458e-05; V401 = -4.3209400767e-06; 
V011 = 1.5355844621e-05; V111 = 2.0914122241e-06; V211 = -5.7751479725e-06; 
V311 = 1.0767234341e-06; V021 = -9.6659393016e-06; V121 = -7.0686982208e-07; 
V221 = 1.4488066593e-06; V031 = 3.1134283336e-06; V131 = 7.9562529879e-08; 
V041 = -5.6590253863e-07; V002 = 1.0500241168e-06; V102 = 1.9600661704e-06; 
V202 = -2.1666693382e-06; V012 = -3.8541359685e-06; V112 = 1.0157632247e-06; 
V022 = 1.7178343158e-06; V003 = -4.1503454190e-07; V103 = 3.5627020989e-07; 
V013 = -1.1293871415e-07; 

vp3 =       V013.*tt+V103.*ss+V003 ;
vp2 =      (V022.*tt+V112.*ss+V012).*tt+(V202*ss+V102).*ss+V002 ;
vp1 =    (((V041.*tt+V131.*ss+V031).*tt ...
     +   (V221.*ss+V121).*ss+V021).*tt ...
     +  ((V311.*ss+V211).*ss+V111).*ss+V011).*tt ...
     + (((V401.*ss+V301).*ss+V201).*ss+V101).*ss+V001 ;
vp0 =  (((((V060.*tt+V150.*ss+V050).*tt ...
     +    (V240.*ss+V140).*ss+V040).*tt ...
     +   ((V330.*ss+V230).*ss+V130).*ss+V030).*tt ...
     +  (((V420.*ss+V320).*ss+V220).*ss+V120).*ss+V020).*tt ...
     + ((((V510.*ss+V410).*ss+V310).*ss+V210).*ss+V110).*ss+V010).*tt ...
     +(((((V600.*ss+V500).*ss+V400).*ss+V300).*ss+V200).*ss+V100).*ss+V000 ;
delta = ( ( vp3.*pp + vp2 ).*pp + vp1 ).*pp + vp0 ;

% specific volume
specvol = v0 + delta;

% alpha
A000 = -3.7374131601e-07; A100 = 7.9665872970e-07; A200 = -9.5176719025e-07; 
A300 = 7.4546183908e-07; A400 = -2.5028304913e-07; A500 = 2.6879827908e-08; 
A010 = 1.3773425769e-06; A110 = -1.8298667100e-06; A210 = 1.7244577313e-06; 
A310 = -8.8316270612e-07; A410 = 1.7982565968e-07; A020 = -1.2380121746e-06; 
A120 = 1.8309269291e-06; A220 = -1.0955055542e-06; A320 = 1.7470054992e-07; 
A030 = 6.7896174634e-07; A130 = -8.7951832993e-07; A230 = 4.4249040774e-07; 
A040 = -9.0669679187e-08; A140 = -4.3350699006e-08; A050 = 2.8562048354e-08; 
A001 = 3.8389611552e-07; A101 = 5.2285305603e-08; A201 = -1.4437869931e-07; 
A301 = 2.6918085852e-08; A011 = -4.8329696508e-07; A111 = -3.5343491104e-08; 
A211 = 7.2440332965e-08; A021 = 2.3350712502e-07; A121 = 5.9671897409e-09; 
A031 = -5.6590253863e-08; A002 = -9.6353399212e-08; A102 = 2.5394080617e-08; 
A012 = 8.5891715792e-08; A003 = -2.8234678537e-09; 

ap3 = A003;
ap2 = A012.*tt + A102.*ss+A002;
ap1 = ((A031.*tt+A121.*ss+A021).*tt	...
	+ (A211.*ss+A111).*ss+A011).*tt	...
	+ ((A301.*ss+A201).*ss+A101).*ss+A001;
ap0 = ((((A050.*tt+A140.*ss+A040).*tt	...
	+ (A230.*ss+A130).*ss+A030).*tt	...
	+ ((A320.*ss+A220).*ss+A120).*ss+A020).*tt	...
	+ (((A410.*ss+A310).*ss+A210).*ss+A110).*ss+A010).*tt	...
	+ ((((A500.*ss+A400).*ss+A300).*ss+A200).*ss+A100).*ss+A000;
a = ( ( ap3.*pp + ap2 ).*pp + ap1 ).*pp + ap0 ;
alpha = a ./ specvol;

% beta
B000 = 3.8896161405e-06; B100 = -1.6824629831e-05; B200 = 3.2146372768e-05; 
B300 = -2.9366928644e-05; B400 = 1.3376886957e-05; B500 = -2.4394186796e-06; 
B010 = -3.9645988657e-07; B110 = 9.4730026353e-07; B210 = -1.1129447472e-06; 
B310 = 4.9821679257e-07; B410 = -6.6884182186e-08; B020 = 4.5531965020e-07; 
B120 = -8.5818216892e-07; B220 = 6.5926332049e-07; B320 = -1.7898168433e-07; 
B030 = -3.0372230734e-07; B130 = 3.6345467351e-07; B230 = -8.6940314119e-08; 
B040 = 1.0942381108e-07; B140 = -1.1010341714e-07; B050 = 4.3147241272e-09; 
B001 = -2.6259371009e-07; B101 = 5.3056825249e-07; B201 = -6.4748391553e-07; 
B301 = 2.1503303093e-07; B011 = -2.6019957550e-08; B111 = 1.4370108710e-07; 
B211 = -4.0187626893e-08; B021 = 8.7944033950e-09; B121 = -3.6050174460e-08; 
B031 = -9.8986399054e-10; B002 = -2.4385837456e-08; B102 = 5.3912512852e-08; 
B012 = -1.2637449319e-08; B003 = -4.4324765969e-09; 

bp3 = B003;
bp2 = B012.*tt + B102.*ss+B002;
bp1 = ((B031.*tt+B121.*ss+B021).*tt	...
	+ (B211.*ss+B111).*ss+B011).*tt	...
	+ ((B301.*ss+B201).*ss+B101).*ss+B001;
bp0 = ((((B050.*tt+B140.*ss+B040).*tt	...
	+ (B230.*ss+B130).*ss+B030).*tt	...
	+ ((B320.*ss+B220).*ss+B120).*ss+B020).*tt	...
	+ (((B410.*ss+B310).*ss+B210).*ss+B110).*ss+B010).*tt	...
	+ ((((B500.*ss+B400).*ss+B300).*ss+B200).*ss+B100).*ss+B000;
b = ( ( bp3.*pp + bp2 ).*pp + bp1 ).*pp + bp0 ;
beta = b ./ ss ./ specvol;


