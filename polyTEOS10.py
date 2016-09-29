import math
import numpy as npy


# polyTEOS10_bsq              in-situ density (55-term polynomial equation)
#==========================================================================
#
# USAGE:
#     [rho,a,b,r0,r] = polyTEOS10_bsq(SA,CT,p)
#
# DESCRIPTION:
#  Calculates in-situ density from Absolute Salinity, Conservative
#  Temperature and pressure, using the computationally-efficient 55-term
#  polynomial expression for density (Roquet et al., 2014).
#
#  Note that the 55-term equation has been fitted in a restricted range of
#  parameter space, and is most accurate inside the "oceanographic funnel"
#  described in McDougall et al. (2011).  The GSW library function
#  "gsw_infunnel(SA,CT,p)" is available to be used if one wants to test if
#  some of one's data lies outside this "funnel".
#
# INPUT:
#  SA  =  Absolute Salinity                                        [ g/kg ]
#  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
#  p   =  sea pressure                                             [ dbar ]
#         ( i.e. absolute pressure - 10.1325 dbar )
#
#  SA & CT & p need to have the same dimensions.
#
# OUTPUT:
#  rho   =  in situ density                                  [ kg/m^3 ]
#  a     =  Boussinesq thermal expansion (=-dr/dCT)        [ kg/m^3/K ]
#  b     =  Boussinesq haline contraction (=dr/dSA)   [ kg/m^3/(g/kg) ]
#  r0    =  vertical reference density                       [ kg/m^3 ]
#  r     =  density anomaly (rho = r0 + r)                   [ kg/m^3 ]
#
# CHECK VALUES (SA=30g/kg, CT=10degC, p=1e3dbar):
#  rho = 1027.45140
#  a   = 0.179646281
#  b   = 0.765555368
#  r0  = 4.59763035
#  r   = 1022.85377
#
# REFERENCES:
#  Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate 
#   polynomial expressions for the density and specific volume of 
#   seawater using the TEOS-10 standard. Ocean Modelling.
#  McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003: 
#   Accurate and computationally efficient algorithms for potential 
#   temperature and density of seawater.  Journal of Atmospheric and 
#   Oceanic Technology, 20, 730-741. 
#
# AUTHOR:
#  Fabien Roquet
#  dec 2015

def polyTEOS10_bsq(SA,CT,p):
    
    # reduced variables
    SAu = 40.*35.16504/35.; CTu = 40.; Zu=1e4; deltaS = 32.
    ss   = npy.sqrt ( (SA+deltaS)/SAu )
    tt   = CT / CTu
    pp   =   p / Zu
    
    # vertical reference profile of density
    R00 = 4.6494977072e+01; R01 = -5.2099962525e+00; R02 = 2.2601900708e-01
    R03 = 6.4326772569e-02; R04 = 1.5616995503e-02; R05 = -1.7243708991e-03
    r0 = ( ( ( ( ( R05*pp + R04 )*pp + R03 )*pp + R02 )*pp + R01 )*pp + R00 )*pp
    
    # density anomaly
    R000 = 8.0189615746e+02; R100 = 8.6672408165e+02; R200 = -1.7864682637e+03
    R300 = 2.0375295546e+03; R400 = -1.2849161071e+03; R500 = 4.3227585684e+02
    R600 = -6.0579916612e+01; R010 = 2.6010145068e+01; R110 = -6.5281885265e+01
    R210 = 8.1770425108e+01; R310 = -5.6888046321e+01; R410 = 1.7681814114e+01
    R510 = -1.9193502195e+00; R020 = -3.7074170417e+01; R120 = 6.1548258127e+01
    R220 = -6.0362551501e+01; R320 = 2.9130021253e+01; R420 = -5.4723692739e+00
    R030 = 2.1661789529e+01; R130 = -3.3449108469e+01; R230 = 1.9717078466e+01
    R330 = -3.1742946532e+00; R040 = -8.3627885467e+00; R140 = 1.1311538584e+01
    R240 = -5.3563304045e+00; R050 = 5.4048723791e-01; R150 = 4.8169980163e-01
    R060 = -1.9083568888e-01; R001 = 1.9681925209e+01; R101 = -4.2549998214e+01
    R201 = 5.0774768218e+01; R301 = -3.0938076334e+01; R401 = 6.6051753097e+00
    R011 = -1.3336301113e+01; R111 = -4.4870114575e+00; R211 = 5.0042598061e+00
    R311 = -6.5399043664e-01; R021 = 6.7080479603e+00; R121 = 3.5063081279e+00
    R221 = -1.8795372996e+00; R031 = -2.4649669534e+00; R131 = -5.5077101279e-01
    R041 = 5.5927935970e-01; R002 = 2.0660924175e+00; R102 = -4.9527603989e+00
    R202 = 2.5019633244e+00; R012 = 2.0564311499e+00; R112 = -2.1311365518e-01
    R022 = -1.2419983026e+00; R003 = -2.3342758797e-02; R103 = -1.8507636718e-02
    R013 = 3.7969820455e-01
    rz3 =       R013*tt + R103*ss + R003
    rz2 =       (R022*tt+R112*ss+R012)*tt+(R202*ss+R102)*ss+R002
    rz1 =       (((R041*tt+R131*ss+R031)*tt + \
                (R221*ss+R121)*ss+R021)*tt + \
                ((R311*ss+R211)*ss+R111)*ss+R011)*tt + \
                (((R401*ss+R301)*ss+R201)*ss+R101)*ss+R001
    rz0 =       (((((R060*tt+R150*ss+R050)*tt + \
                (R240*ss+R140)*ss+R040)*tt + \
                ((R330*ss+R230)*ss+R130)*ss+R030)*tt + \
                (((R420*ss+R320)*ss+R220)*ss+R120)*ss+R020)*tt + \
                ((((R510*ss+R410)*ss+R310)*ss+R210)*ss+R110)*ss+R010)*tt + \
                (((((R600*ss+R500)*ss+R400)*ss+R300)*ss+R200)*ss+R100)*ss+R000
    r = ( ( rz3*pp + rz2 )*pp + rz1 )*pp + rz0;
    
    # in-situ density
    rho = r + r0
    
    # thermal expansion a
    ALP000 = -6.5025362670e-01; ALP100 = 1.6320471316e+00; ALP200 = -2.0442606277e+00
    ALP300 = 1.4222011580e+00; ALP400 = -4.4204535284e-01; ALP500 = 4.7983755487e-02
    ALP010 = 1.8537085209e+00; ALP110 = -3.0774129064e+00; ALP210 = 3.0181275751e+00
    ALP310 = -1.4565010626e+00; ALP410 = 2.7361846370e-01; ALP020 = -1.6246342147e+00
    ALP120 = 2.5086831352e+00; ALP220 = -1.4787808849e+00; ALP320 = 2.3807209899e-01
    ALP030 = 8.3627885467e-01; ALP130 = -1.1311538584e+00; ALP230 = 5.3563304045e-01
    ALP040 = -6.7560904739e-02; ALP140 = -6.0212475204e-02; ALP050 = 2.8625353333e-02
    ALP001 = 3.3340752782e-01; ALP101 = 1.1217528644e-01; ALP201 = -1.2510649515e-01
    ALP301 = 1.6349760916e-02; ALP011 = -3.3540239802e-01; ALP111 = -1.7531540640e-01
    ALP211 = 9.3976864981e-02; ALP021 = 1.8487252150e-01; ALP121 = 4.1307825959e-02
    ALP031 = -5.5927935970e-02; ALP002 = -5.1410778748e-02; ALP102 = 5.3278413794e-03
    ALP012 = 6.2099915132e-02; ALP003 = -9.4924551138e-03
    
    a = ((ALP003*pp + \
          ALP012*tt + ALP102*ss+ALP002)*pp	+ \
         ((ALP031*tt + \
           ALP121*ss+ALP021)*tt + \
          (ALP211*ss+ALP111)*ss+ALP011)*tt	+ \
         ((ALP301*ss+ALP201)*ss+ALP101)*ss+ALP001)*pp + \
        ((((ALP050*tt + \
            ALP140*ss+ALP040)*tt + \
           (ALP230*ss+ALP130)*ss+ALP030)*tt	+ \
          ((ALP320*ss+ALP220)*ss+ALP120)*ss+ALP020)*tt + \
         (((ALP410*ss+ALP310)*ss+ALP210)*ss+ALP110)*ss+ALP010)*tt + \
        ((((ALP500*ss+ALP400)*ss+ALP300)*ss+ALP200)*ss+ALP100)*ss+ALP000
    
    
    # haline contraction b
    BET000 = 1.0783203594e+01; BET100 = -4.4452095908e+01; BET200 = 7.6048755820e+01
    BET300 = -6.3944280668e+01; BET400 = 2.6890441098e+01; BET500 = -4.5221697773e+00
    BET010 = -8.1219372432e-01; BET110 = 2.0346663041e+00; BET210 = -2.1232895170e+00
    BET310 = 8.7994140485e-01; BET410 = -1.1939638360e-01; BET020 = 7.6574242289e-01
    BET120 = -1.5019813020e+00; BET220 = 1.0872489522e+00; BET320 = -2.7233429080e-01
    BET030 = -4.1615152308e-01; BET130 = 4.9061350869e-01; BET230 = -1.1847737788e-01
    BET040 = 1.4073062708e-01; BET140 = -1.3327978879e-01; BET050 = 5.9929880134e-03
    BET001 = -5.2937873009e-01; BET101 = 1.2634116779e+00; BET201 = -1.1547328025e+00
    BET301 = 3.2870876279e-01; BET011 = -5.5824407214e-02; BET111 = 1.2451933313e-01
    BET211 = -2.4409539932e-02; BET021 = 4.3623149752e-02; BET121 = -4.6767901790e-02
    BET031 = -6.8523260060e-03; BET002 = -6.1618945251e-02; BET102 = 6.2255521644e-02
    BET012 = -2.6514181169e-03; BET003 = -2.3025968587e-04
    
    b = ((BET003*pp + \
          BET012*tt + BET102*ss+BET002)*pp + \
         ((BET031*tt + \
           BET121*ss+BET021)*tt + \
          (BET211*ss+BET111)*ss+BET011)*tt + \
         ((BET301*ss+BET201)*ss+BET101)*ss+BET001)*pp + \
        ((((BET050*tt + \
            BET140*ss+BET040)*tt + \
           (BET230*ss+BET130)*ss+BET030)*tt + \
          ((BET320*ss+BET220)*ss+BET120)*ss+BET020)*tt + \
         (((BET410*ss+BET310)*ss+BET210)*ss+BET110)*ss+BET010)*tt + \
        ((((BET500*ss+BET400)*ss+BET300)*ss+BET200)*ss+BET100)*ss+BET000
    
    b = b / ss ;
    
    return rho,a,b,r0,r
    
    
# polyTEOS10_stif             in-situ density (55-term polynomial equation)
#==========================================================================
#
# USAGE:
#     [rho,a,b,r1,rdot] = polyTEOS10_stif(SA,CT,p)
#
# DESCRIPTION:
#  Calculates in-situ density from Absolute Salinity, Conservative
#  Temperature and pressure, using the computationally-efficient 55-term
#  polynomial expression for density (Roquet et al., 2014).
#
#  Note that the 55-term equation has been fitted in a restricted range of
#  parameter space, and is most accurate inside the "oceanographic funnel"
#  described in McDougall et al. (2011).  The GSW library function
#  "gsw_infunnel(SA,CT,p)" is available to be used if one wants to test if
#  some of one's data lies outside this "funnel".
#
# INPUT:
#  SA  =  Absolute Salinity                                        [ g/kg ]
#  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
#  p   =  sea pressure                                             [ dbar ]
#         ( i.e. absolute pressure - 10.1325 dbar )
#
#  SA & CT & p need to have the same dimensions.
#
# OUTPUT:
#  rho   =  in situ density                                  [ kg/m^3 ]
#  a     =  Boussinesq thermal expansion (=-drho/dCT)      [ kg/m^3/K ]
#  b     =  Boussinesq haline contraction (=drho/dSA) [ kg/m^3/(g/kg) ]
#  r1    =  vertical reference ratio                             [ SI ]
#  rdot  =  stiffened density (rho = r1 x rdot)              [ kg/m^3 ]
#
# CHECK VALUES (SA=30g/kg, CT=10degC, p=1e3dbar):
#  rho = 1027.45140
#  a   = 0.179649406
#  b   = 0.765554495
#  r1  = 1.00447333
#  rdot= 1022.87574
#
# REFERENCES:
#  Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate 
#   polynomial expressions for the density and specific volume of 
#   seawater using the TEOS-10 standard. Ocean Modelling.
#  McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003: 
#   Accurate and computationally efficient algorithms for potential 
#   temperature and density of seawater.  Journal of Atmospheric and 
#   Oceanic Technology, 20, 730-741. 
#
# AUTHOR:
#  Fabien Roquet
#  jan 2015

def polyTEOS10_stif(SA,CT,p):
    
    # reduced variables
    SAu = 40.*35.16504/35.; CTu = 40.; Zu=1e4; deltaS = 32.
    ss   = npy.sqrt ( (SA+deltaS)/SAu )
    tt   = CT / CTu
    pp   =   p / Zu
    
    # vertical reference profile of density
    R10 = 4.5238001132e-02; R11 = -5.0691457704e-03; R12 = 2.1990865986e-04; 
    R13 = 6.2587720090e-05; R14 = 1.5194795322e-05; R15 = -1.6777531159e-06;
    r1 = ( ( ( ( ( R15*pp + R14 )*pp + R13 )*pp + R12 )*pp + R11 )*pp + R10 )*pp + 1. ;

    # stiffened density
    R000 = 8.0185969881e+02; R100 = 8.6694399997e+02; R200 = -1.7869886805e+03; 
    R300 = 2.0381548497e+03; R400 = -1.2853207957e+03; R500 = 4.3240996619e+02; 
    R600 = -6.0597695001e+01; R010 = 2.6018938392e+01; R110 = -6.5349779146e+01; 
    R210 = 8.1938301569e+01; R310 = -5.7075042739e+01; R410 = 1.7778970855e+01; 
    R510 = -1.9385269480e+00; R020 = -3.7047586837e+01; R120 = 6.1469677558e+01; 
    R220 = -6.0273564480e+01; R320 = 2.9086147388e+01; R420 = -5.4641145446e+00; 
    R030 = 2.1645370860e+01; R130 = -3.3415215649e+01; R230 = 1.9694119706e+01; 
    R330 = -3.1710494147e+00; R040 = -8.3587258634e+00; R140 = 1.1301873278e+01; 
    R240 = -5.3494903247e+00; R050 = 5.4258499460e-01; R150 = 4.7964098705e-01; 
    R060 = -1.9098981559e-01; R001 = 2.1989266031e+01; R101 = -4.2043785414e+01; 
    R201 = 4.8565183521e+01; R301 = -3.0473875108e+01; R401 = 6.5025796369e+00; 
    R011 = -1.3731593003e+01; R111 = -4.3667263842e+00; R211 = 5.2899298884e+00; 
    R311 = -7.1323826203e-01; R021 = 7.4843325711e+00; R121 = 3.1442996192e+00; 
    R221 = -1.8141771987e+00; R031 = -2.6010182316e+00; R131 = -4.9866739215e-01; 
    R041 = 5.5882364387e-01; R002 = 1.1144125393e+00; R102 = -4.5413502768e+00; 
    R202 = 2.7242121539e+00; R012 = 2.8508446713e+00; R112 = -4.4471361300e-01; 
    R022 = -1.5059302816e+00; R003 = 1.9817079368e-01; R103 = -1.7905369937e-01; 
    R013 = 2.5254165600e-01; 
    rz3 =       R013*tt + R103*ss + R003;
    rz2 =      (R022*tt+R112*ss+R012)*tt+(R202*ss+R102)*ss+R002;
    rz1 =    (((R041*tt+R131*ss+R031)*tt \
        +   (R221*ss+R121)*ss+R021)*tt \
        +  ((R311*ss+R211)*ss+R111)*ss+R011)*tt \
        + (((R401*ss+R301)*ss+R201)*ss+R101)*ss+R001;
    rz0 =  (((((R060*tt+R150*ss+R050)*tt \
        +    (R240*ss+R140)*ss+R040)*tt \
        +   ((R330*ss+R230)*ss+R130)*ss+R030)*tt \
        +  (((R420*ss+R320)*ss+R220)*ss+R120)*ss+R020)*tt \
        + ((((R510*ss+R410)*ss+R310)*ss+R210)*ss+R110)*ss+R010)*tt \
        +(((((R600*ss+R500)*ss+R400)*ss+R300)*ss+R200)*ss+R100)*ss+R000;
    rdot = ( ( rz3*pp + rz2 )*pp + rz1 )*pp + rz0;

    # in-situ density
    rho = r1 * rdot ;

    # thermal expansion a
    ALP000 = -6.5047345980e-01; ALP100 = 1.6337444787e+00; ALP200 = -2.0484575392e+00; 
    ALP300 = 1.4268760685e+00; ALP400 = -4.4447427136e-01; ALP500 = 4.8463173700e-02; 
    ALP010 = 1.8523793418e+00; ALP110 = -3.0734838779e+00; ALP210 = 3.0136782240e+00; 
    ALP310 = -1.4543073694e+00; ALP410 = 2.7320572723e-01; ALP020 = -1.6234028145e+00; 
    ALP120 = 2.5061411737e+00; ALP220 = -1.4770589780e+00; ALP320 = 2.3782870611e-01; 
    ALP030 = 8.3587258634e-01; ALP130 = -1.1301873278e+00; ALP230 = 5.3494903247e-01; 
    ALP040 = -6.7823124325e-02; ALP140 = -5.9955123381e-02; ALP050 = 2.8648472338e-02; 
    ALP001 = 3.4328982507e-01; ALP101 = 1.0916815960e-01; ALP201 = -1.3224824721e-01; 
    ALP301 = 1.7830956551e-02; ALP011 = -3.7421662855e-01; ALP111 = -1.5721498096e-01; 
    ALP211 = 9.0708859933e-02; ALP021 = 1.9507636737e-01; ALP121 = 3.7400054411e-02; 
    ALP031 = -5.5882364387e-02; ALP002 = -7.1271116782e-02; ALP102 = 1.1117840325e-02; 
    ALP012 = 7.5296514078e-02; ALP003 = -6.3135413999e-03; 
    a = r1 * ( \
        ((ALP003*pp	\
        + ALP012*tt + ALP102*ss+ALP002)*pp	\
        + ((ALP031*tt	\
        + ALP121*ss+ALP021)*tt	\
        + (ALP211*ss+ALP111)*ss+ALP011)*tt	\
        + ((ALP301*ss+ALP201)*ss+ALP101)*ss+ALP001)*pp	\
        + ((((ALP050*tt	\
        + ALP140*ss+ALP040)*tt	\
        + (ALP230*ss+ALP130)*ss+ALP030)*tt	\
        + ((ALP320*ss+ALP220)*ss+ALP120)*ss+ALP020)*tt	\
        + (((ALP410*ss+ALP310)*ss+ALP210)*ss+ALP110)*ss+ALP010)*tt	\
        + ((((ALP500*ss+ALP400)*ss+ALP300)*ss+ALP200)*ss+ALP100)*ss+ALP000 );

    # haline contraction b
    BET000 = 1.0785939671e+01; BET100 = -4.4465045269e+01; BET200 = 7.6072094337e+01; 
    BET300 = -6.3964420131e+01; BET400 = 2.6898783594e+01; BET500 = -4.5234968986e+00; 
    BET010 = -8.1303841476e-01; BET110 = 2.0388435182e+00; BET210 = -2.1302689715e+00; 
    BET310 = 8.8477644261e-01; BET410 = -1.2058930400e-01; BET020 = 7.6476477580e-01; 
    BET120 = -1.4997670675e+00; BET220 = 1.0856114040e+00; BET320 = -2.7192349143e-01; 
    BET030 = -4.1572985119e-01; BET130 = 4.9004223351e-01; BET230 = -1.1835625260e-01; 
    BET040 = 1.4061037779e-01; BET140 = -1.3310958936e-01; BET050 = 5.9673736141e-03; 
    BET001 = -5.2308076768e-01; BET101 = 1.2084313165e+00; BET201 = -1.1374069553e+00; 
    BET301 = 3.2360305476e-01; BET011 = -5.4327900468e-02; BET111 = 1.3162756682e-01; 
    BET211 = -2.6620905846e-02; BET021 = 3.9119281064e-02; BET121 = -4.5141568126e-02; 
    BET031 = -6.2040874705e-03; BET002 = -5.6500454603e-02; BET102 = 6.7785665385e-02; 
    BET012 = -5.5328304955e-03; BET003 = -2.2276668383e-03; 
    b = ((BET003*pp	\
        + BET012*tt + BET102*ss+BET002)*pp	\
        + ((BET031*tt	\
        + BET121*ss+BET021)*tt	\
        + (BET211*ss+BET111)*ss+BET011)*tt	\
        + ((BET301*ss+BET201)*ss+BET101)*ss+BET001)*pp	\
        + ((((BET050*tt	\
        + BET140*ss+BET040)*tt	\
        + (BET230*ss+BET130)*ss+BET030)*tt	\
        + ((BET320*ss+BET220)*ss+BET120)*ss+BET020)*tt	\
        + (((BET410*ss+BET310)*ss+BET210)*ss+BET110)*ss+BET010)*tt	\
        + ((((BET500*ss+BET400)*ss+BET300)*ss+BET200)*ss+BET100)*ss+BET000;
    b = b * r1 / ss ;

    return rho,a,b,r1,rdot
    
    
# polyTEOS10_55t              specific volume (55-term polynomial equation)
#==========================================================================
#
# USAGE:
#     [specvol,alpha,beta,v0,delta] = polyTEOS10_55t(SA,CT,p)
#
# DESCRIPTION:
#  Calculates specific volume from Absolute Salinity, Conservative
#  Temperature and pressure, using the computationally-efficient 55-term
#  polynomial expression for specific volume (Roquet et al., 2014).
#
#  Note that the 55-term equation has been fitted in a restricted range of
#  parameter space, and is most accurate inside the "oceanographic funnel"
#  described in McDougall et al. (2011).  The GSW library function
#  "gsw_infunnel(SA,CT,p)" is available to be used if one wants to test if
#  some of one's data lies outside this "funnel".
#
# INPUT:
#  SA  =  Absolute Salinity                                        [ g/kg ]
#  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
#  p   =  sea pressure                                             [ dbar ]
#         ( i.e. absolute pressure - 10.1325 dbar )
#
#  SA & CT & p need to have the same dimensions.
#
# OUTPUT:
#  specvol   =  specific volume                                  [ m^3/kg ]
#  alpha     =  Thermal expansion  (=v.dv/dCT)                      [ 1/K ]
#  beta      =  Haline contraction (=-v.dv/dSA)                [ 1/(g/kg) ]
#  v0        =  vertical reference specific volume               [ m^3/kg ]
#  delta     =  specific volume anomaly (specvol = v0 + delta)   [ m^3/kg ]
#
# CHECK VALUES (SA=30g/kg, CT=10degC, p=1e3dbar):
#  specvol = 9.732820466e-04
#  alpha   = 1.748553121e-04
#  beta    = 7.450974025e-04
#  v0      = -4.333016903e-06
#  delta   = 9.776150635e-04
#
# REFERENCES:
#  Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate 
#   polynomial expressions for the density and specific volume of 
#   seawater using the TEOS-10 standard. Ocean Modelling.
#  McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003: 
#   Accurate and computationally efficient algorithms for potential 
#   temperature and density of seawater.  Journal of Atmospheric and 
#   Oceanic Technology, 20, 730-741. 
#
# AUTHOR:
#  Fabien Roquet
#  jan 2015

def polyTEOS10_55t(SA,CT,p):
    
    # reduced variables
    SAu = 40.*35.16504/35.; CTu = 40.; Zu=1e4; deltaS = 32.
    ss   = npy.sqrt ( (SA+deltaS)/SAu )
    tt   = CT / CTu
    pp   =   p / Zu
    
    # vertical reference profile of specific volume
    V00 = -4.4015007269e-05; V01 = 6.9232335784e-06; V02 = -7.5004675975e-07; 
    V03 = 1.7009109288e-08; V04 = -1.6884162004e-08; V05 = 1.9613503930e-09; 
    v0  = (((((V05*pp+V04)*pp+V03 )*pp+V02 )*pp+V01)*pp+V00)*pp ;

    # specific volume anomaly
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

    vp3 =       V013*tt+V103*ss+V003 ;
    vp2 =      (V022*tt+V112*ss+V012)*tt+(V202*ss+V102)*ss+V002 ;
    vp1 =    (((V041*tt+V131*ss+V031)*tt \
         +   (V221*ss+V121)*ss+V021)*tt \
         +  ((V311*ss+V211)*ss+V111)*ss+V011)*tt \
         + (((V401*ss+V301)*ss+V201)*ss+V101)*ss+V001 ;
    vp0 =  (((((V060*tt+V150*ss+V050)*tt \
         +    (V240*ss+V140)*ss+V040)*tt \
         +   ((V330*ss+V230)*ss+V130)*ss+V030)*tt \
         +  (((V420*ss+V320)*ss+V220)*ss+V120)*ss+V020)*tt \
         + ((((V510*ss+V410)*ss+V310)*ss+V210)*ss+V110)*ss+V010)*tt \
         +(((((V600*ss+V500)*ss+V400)*ss+V300)*ss+V200)*ss+V100)*ss+V000 ;
    delta = ( ( vp3*pp + vp2 )*pp + vp1 )*pp + vp0 ;

    # specific volume
    specvol = v0 + delta;

    # alpha
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
    ap2 = A012*tt + A102*ss+A002;
    ap1 = ((A031*tt+A121*ss+A021)*tt	\
        + (A211*ss+A111)*ss+A011)*tt	\
        + ((A301*ss+A201)*ss+A101)*ss+A001;
    ap0 = ((((A050*tt+A140*ss+A040)*tt	\
        + (A230*ss+A130)*ss+A030)*tt	\
        + ((A320*ss+A220)*ss+A120)*ss+A020)*tt	\
        + (((A410*ss+A310)*ss+A210)*ss+A110)*ss+A010)*tt	\
        + ((((A500*ss+A400)*ss+A300)*ss+A200)*ss+A100)*ss+A000;
    a = ( ( ap3*pp + ap2 )*pp + ap1 )*pp + ap0 ;
    alpha = a / specvol;

    # beta
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
    bp2 = B012*tt + B102*ss+B002;
    bp1 = ((B031*tt+B121*ss+B021)*tt	\
        + (B211*ss+B111)*ss+B011)*tt	\
        + ((B301*ss+B201)*ss+B101)*ss+B001;
    bp0 = ((((B050*tt+B140*ss+B040)*tt	\
        + (B230*ss+B130)*ss+B030)*tt	\
        + ((B320*ss+B220)*ss+B120)*ss+B020)*tt	\
        + (((B410*ss+B310)*ss+B210)*ss+B110)*ss+B010)*tt	\
        + ((((B500*ss+B400)*ss+B300)*ss+B200)*ss+B100)*ss+B000;
    b = ( ( bp3*pp + bp2 )*pp + bp1 )*pp + bp0 ;
    beta = b / ss / specvol;

    return specvol,alpha,beta,v0,delta
    
    
# polyTEOS10_75t              specific volume (75-term polynomial equation)
#==========================================================================
#
# USAGE:
#     [specvol,alpha,beta,v0,delta] = polyTEOS10_75t(SA,CT,p)
#
# DESCRIPTION:
#  Calculates specific volume from Absolute Salinity, Conservative
#  Temperature and pressure, using the computationally-efficient 75-term
#  polynomial expression for density (Roquet et al., 2014).
#
#  Note that the 75-term equation has been fitted in a restricted range of
#  parameter space, and is most accurate inside the "oceanographic funnel"
#  described in McDougall et al. (2011).  The GSW library function
#  "gsw_infunnel(SA,CT,p)" is available to be used if one wants to test if
#  some of one's data lies outside this "funnel".
#
# INPUT:
#  SA  =  Absolute Salinity                                        [ g/kg ]
#  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
#  p   =  sea pressure                                             [ dbar ]
#         ( i.e. absolute pressure - 10.1325 dbar )
#
#  SA & CT & p need to have the same dimensions.
#
# OUTPUT:
#  specvol   =  specific volume                                  [ m^3/kg ]
#  alpha     =  Thermal expansion  (=v.dv/dCT)                      [ 1/K ]
#  beta      =  Haline contraction (=-v.dv/dSA)                [ 1/(g/kg) ]
#  v0        =  vertical reference specific volume               [ m^3/kg ]
#  delta     =  specific volume anomaly (specvol = v0 + delta)   [ m^3/kg ]
#
# CHECK VALUES (SA=30g/kg, CT=10degC, p=1e3dbar):
#  specvol = 9.732819628e-04
#  alpha   = 1.748439401e-04
#  beta    = 7.451213159e-04
#  v0      = -4.333016903e-06
#  delta   = 9.776149797e-04
#
# REFERENCES:
#  Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate 
#   polynomial expressions for the density and specific volume of 
#   seawater using the TEOS-10 standard. Ocean Modelling.
#  McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003: 
#   Accurate and computationally efficient algorithms for potential 
#   temperature and density of seawater.  Journal of Atmospheric and 
#   Oceanic Technology, 20, 730-741. 
#
# AUTHOR:
#  Fabien Roquet
#  jan 2015

def polyTEOS10_75t(SA,CT,p):
    
    # reduced variables
    SAu = 40.*35.16504/35.; CTu = 40.; Zu=1e4; deltaS = 32.
    ss   = npy.sqrt ( (SA+deltaS)/SAu )
    tt   = CT / CTu
    pp   =   p / Zu
    
    # vertical reference profile of specific volume
    V00 = -4.4015007269e-05; V01 = 6.9232335784e-06; V02 = -7.5004675975e-07;
    V03 = 1.7009109288e-08; V04 = -1.6884162004e-08; V05 = 1.9613503930e-09;
    v0  = (((((V05*pp+V04)*pp+V03 )*pp+V02 )*pp+V01)*pp+V00)*pp ;

    # specific volume anomaly
    V000 = 1.0769995862e-03; V100 = -3.1038981976e-04; V200 = 6.6928067038e-04; 
    V300 = -8.5047933937e-04; V400 = 5.8086069943e-04; V500 = -2.1092370507e-04; 
    V600 = 3.1932457305e-05; V010 = -1.5649734675e-05; V110 = 3.5009599764e-05; 
    V210 = -4.3592678561e-05; V310 = 3.4532461828e-05; V410 = -1.1959409788e-05; 
    V510 = 1.3864594581e-06; V020 = 2.7762106484e-05; V120 = -3.7435842344e-05; 
    V220 = 3.5907822760e-05; V320 = -1.8698584187e-05; V420 = 3.8595339244e-06; 
    V030 = -1.6521159259e-05; V130 = 2.4141479483e-05; V230 = -1.4353633048e-05; 
    V330 = 2.2863324556e-06; V040 = 6.9111322702e-06; V140 = -8.7595873154e-06; 
    V240 = 4.3703680598e-06; V050 = -8.0539615540e-07; V150 = -3.3052758900e-07; 
    V060 = 2.0543094268e-07; V001 = -1.6784136540e-05; V101 = 2.4262468747e-05; 
    V201 = -3.4792460974e-05; V301 = 3.7470777305e-05; V401 = -1.7322218612e-05; 
    V501 = 3.0927427253e-06; V011 = 1.8505765429e-05; V111 = -9.5677088156e-06; 
    V211 = 1.1100834765e-05; V311 = -9.8447117844e-06; V411 = 2.5909225260e-06; 
    V021 = -1.1716606853e-05; V121 = -2.3678308361e-07; V221 = 2.9283346295e-06; 
    V321 = -4.8826139200e-07; V031 = 7.9279656173e-06; V131 = -3.4558773655e-06; 
    V231 = 3.1655306078e-07; V041 = -3.4102187482e-06; V141 = 1.2956717783e-06; 
    V051 = 5.0736766814e-07; V002 = 3.0623833435e-06; V102 = -5.8484432984e-07; 
    V202 = -4.8122251597e-06; V302 = 4.9263106998e-06; V402 = -1.7811974727e-06; 
    V012 = -1.1736386731e-06; V112 = -5.5699154557e-06; V212 = 5.4620748834e-06; 
    V312 = -1.3544185627e-06; V022 = 2.1305028740e-06; V122 = 3.9137387080e-07; 
    V222 = -6.5731104067e-07; V032 = -4.6132540037e-07; V132 = 7.7618888092e-09; 
    V042 = -6.3352916514e-08; V003 = -3.8088938393e-07; V103 = 3.6310188515e-07; 
    V203 = 1.6746303780e-08; V013 = -3.6527006553e-07; V113 = -2.7295696237e-07; 
    V023 = 2.8695905159e-07; V004 = 8.8302421514e-08; V104 = -1.1147125423e-07; 
    V014 = 3.1454099902e-07; V005 = 4.2369007180e-09; 

    vp5 = V005 ;
    vp4 = V014*tt + V104*ss + V004 ;
    vp3 = ( V023*tt + V113*ss + V013 )*tt \
        + ( V203*ss + V103 )*ss + V003 ;
    vp2 = ( ( ( V042*tt + V132*ss + V032 )*tt \
        + ( V222*ss + V122 )*ss + V022 )*tt \
        + ( ( V312*ss + V212 )*ss + V112 )*ss + V012 )*tt \
        + ( ( ( V402*ss + V302 )*ss + V202 )*ss + V102 )*ss + V002 ;
    vp1 = ( ( ( ( V051*tt + V141*ss + V041 )*tt \
        + ( V231*ss + V131 )*ss + V031 )*tt \
        + ( ( V321*ss + V221 )*ss + V121 )*ss + V021 )*tt \
        + ( ( ( V411*ss + V311 )*ss + V211 )*ss + V111 )*ss + V011 )*tt \
        + ( ( ( ( V501*ss + V401 )*ss + V301 )*ss + V201 )*ss + V101 )*ss + V001 ;
    vp0 = ( ( ( ( ( V060*tt + V150*ss + V050 )*tt \
        + ( V240*ss + V140 )*ss + V040 )*tt \
        + ( ( V330*ss + V230 )*ss + V130 )*ss + V030 )*tt \
        + ( ( ( V420*ss + V320 )*ss + V220 )*ss + V120 )*ss + V020 )*tt \
        + ( ( ( ( V510*ss + V410 )*ss + V310 )*ss + V210 )*ss + V110 )*ss + V010 )*tt \
        + ((((( V600*ss + V500 )*ss + V400 )*ss + V300 )*ss + V200 )*ss + V100 )*ss + V000 ;
    delta = ( ( ( ( vp5*pp + vp4 )*pp + vp3 )*pp + vp2 )*pp + vp1 )*pp + vp0 ;

    # specific volume
    specvol = v0 + delta ;

    # alpha
    A000 = -3.9124336688e-07; A100 = 8.7523999410e-07; A200 = -1.0898169640e-06; 
    A300 = 8.6331154570e-07; A400 = -2.9898524469e-07; A500 = 3.4661486454e-08; 
    A010 = 1.3881053242e-06; A110 = -1.8717921172e-06; A210 = 1.7953911380e-06; 
    A310 = -9.3492920933e-07; A410 = 1.9297669622e-07; A020 = -1.2390869444e-06; 
    A120 = 1.8106109612e-06; A220 = -1.0765224786e-06; A320 = 1.7147493417e-07; 
    A030 = 6.9111322702e-07; A130 = -8.7595873154e-07; A230 = 4.3703680598e-07; 
    A040 = -1.0067451943e-07; A140 = -4.1315948624e-08; A050 = 3.0814641402e-08; 
    A001 = 4.6264413572e-07; A101 = -2.3919272039e-07; A201 = 2.7752086911e-07; 
    A301 = -2.4611779461e-07; A401 = 6.4773063150e-08; A011 = -5.8583034263e-07; 
    A111 = -1.1839154180e-08; A211 = 1.4641673148e-07; A311 = -2.4413069600e-08; 
    A021 = 5.9459742130e-07; A121 = -2.5919080242e-07; A221 = 2.3741479559e-08; 
    A031 = -3.4102187482e-07; A131 = 1.2956717783e-07; A041 = 6.3420958518e-08; 
    A002 = -2.9340966828e-08; A102 = -1.3924788639e-07; A202 = 1.3655187208e-07; 
    A302 = -3.3860464067e-08; A012 = 1.0652514370e-07; A112 = 1.9568693540e-08; 
    A212 = -3.2865552033e-08; A022 = -3.4599405028e-08; A122 = 5.8214166069e-10; 
    A032 = -6.3352916514e-09; A003 = -9.1317516382e-09; A103 = -6.8239240593e-09; 
    A013 = 1.4347952579e-08; A004 = 7.8635249756e-09; 

    ap4 = A004;
    ap3 = A013*tt+A103*ss+A003;
    ap2 = ((A032*tt+A122*ss+A022)*tt	\
        + (A212*ss+A112)*ss+A012)*tt	\
        + ((A302*ss+A202)*ss+A102)*ss+A002;
    ap1 = (((A041*tt+A131*ss+A031)*tt	\
        + (A221*ss+A121)*ss+A021)*tt	\
        + ((A311*ss+A211)*ss+A111)*ss+A011)*tt	\
        + (((A401*ss+A301)*ss+A201)*ss+A101)*ss+A001;
    ap0 = ((((A050*tt+A140*ss+A040)*tt	\
        + (A230*ss+A130)*ss+A030)*tt	\
        + ((A320*ss+A220)*ss+A120)*ss+A020)*tt	\
        + (((A410*ss+A310)*ss+A210)*ss+A110)*ss+A010)*tt	\
        + ((((A500*ss+A400)*ss+A300)*ss+A200)*ss+A100)*ss+A000;

    a = ( ( ( ap4*pp + ap3 )*pp + ap2 )*pp + ap1 )*pp + ap0 ;
    alpha = a / specvol;

    # beta
    B000 = 3.8616633493e-06; B100 = -1.6653488424e-05; B200 = 3.1743292000e-05; 
    B300 = -2.8906727363e-05; B400 = 1.3120861084e-05; B500 = -2.3836941584e-06; 
    B010 = -4.3556611614e-07; B110 = 1.0847021286e-06; B210 = -1.2888896515e-06; 
    B310 = 5.9516403589e-07; B410 = -8.6247024451e-08; B020 = 4.6575180991e-07; 
    B120 = -8.9348241648e-07; B220 = 6.9790598120e-07; B320 = -1.9207099914e-07; 
    B030 = -3.0035220418e-07; B130 = 3.5715667939e-07; B230 = -8.5335075632e-08; 
    B040 = 1.0898094956e-07; B140 = -1.0874641554e-07; B050 = 4.1122040579e-09; 
    B001 = -3.0185747199e-07; B101 = 8.6572923996e-07; B201 = -1.3985593422e-06; 
    B301 = 8.6204601419e-07; B401 = -1.9238922270e-07; B011 = 1.1903505888e-07; 
    B111 = -2.7621838107e-07; B211 = 3.6744403581e-07; B311 = -1.2893812777e-07; 
    B021 = 2.9458973764e-09; B121 = -7.2864777086e-08; B221 = 1.8223868848e-08; 
    B031 = 4.2995723805e-08; B131 = -7.8766845761e-09; B041 = -1.6119885062e-08; 
    B002 = 7.2762435164e-09; B102 = 1.1974099887e-07; B202 = -1.8386962715e-07; 
    B302 = 8.8641889138e-08; B012 = 6.9297177307e-08; B112 = -1.3591099350e-07; 
    B212 = 5.0552320246e-08; B022 = -4.8692129590e-09; B122 = 1.6355652107e-08; 
    B032 = -9.6568249433e-11; B003 = -4.5174717490e-09; B103 = -4.1669270980e-10; 
    B013 = 3.3959486762e-09; B004 = 1.3868510806e-09; 

    bp4 = B004;
    bp3 = B013*tt+B103*ss+B003;
    bp2 = ((B032*tt+B122*ss+B022)*tt	\
        + (B212*ss+B112)*ss+B012)*tt	\
        + ((B302*ss+B202)*ss+B102)*ss+B002;
    bp1 = (((B041*tt+B131*ss+B031)*tt	\
        + (B221*ss+B121)*ss+B021)*tt	\
        + ((B311*ss+B211)*ss+B111)*ss+B011)*tt	\
        + (((B401*ss+B301)*ss+B201)*ss+B101)*ss+B001;
    bp0 = ((((B050*tt+B140*ss+B040)*tt	\
        + (B230*ss+B130)*ss+B030)*tt	\
        + ((B320*ss+B220)*ss+B120)*ss+B020)*tt	\
        + (((B410*ss+B310)*ss+B210)*ss+B110)*ss+B010)*tt	\
        + ((((B500*ss+B400)*ss+B300)*ss+B200)*ss+B100)*ss+B000;

    b = ( ( ( bp4*pp + bp3 )*pp + bp2 )*pp + bp1 )*pp + bp0 ;
    beta = b / ss / specvol;
    
    return specvol,alpha,beta,v0,delta







    
    










