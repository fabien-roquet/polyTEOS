 MODULE polyTEOS10_bsq
   !!==============================================================================
   !!                       ***  MODULE  polyTEOS10_bsq  ***
   !! Ocean diagnostic variable : equation of state - in situ density
   !!----------------------------------------------------------------------
   !!   polyTEOS10_insitu     : Compute the in situ density
   !!   polyTEOS10_rab        : compute in situ thermal expansion and haline contraction
   !!   polyTEOS10_pt_from_ct : convert CT into PT
   !!   polyTEOS10_init       : set eos parameters
   !!----------------------------------------------------------------------
   !!
   !! DESCRIPTION:
   !!  Calculates in-situ density from Absolute Salinity, Conservative
   !!  Temperature and pressure, using the computationally-efficient 55-term
   !!  polynomial expression for density (Roquet et al., 2014).
   !!
   !!  Note that the 55-term equation has been fitted in a restricted range of
   !!  parameter space, and is most accurate inside the "oceanographic funnel"
   !!  described in McDougall et al. (2003).  The GSW library function
   !!  "gsw_infunnel(SA,CT,p)" is available to be used if one wants to test if
   !!  some of one's data lies outside this "funnel".
   !!
   !! INPUT:
   !!  SA  =  Absolute Salinity                                        [ g/kg ]
   !!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
   !!  p   =  sea pressure                                             [ dbar ]
   !!         ( i.e. absolute pressure - 10.1325 dbar )
   !!
   !!  SA & CT & p need to have the same dimensions.
   !!
   !! OUTPUT:
   !!  rho   =  in situ density                                  [ kg/m^3 ]
   !!  a     =  Boussinesq thermal expansion (=-dr/dCT)        [ kg/m^3/K ]
   !!  b     =  Boussinesq haline contraction (=dr/dSA)   [ kg/m^3/(g/kg) ]
   !!  r0    =  vertical reference density                       [ kg/m^3 ]
   !!  r     =  density anomaly (rho = r0 + r)                   [ kg/m^3 ]
   !!
   !! CHECK VALUES (SA=30g/kg, CT=10degC, p=1e3dbar):
   !!  rho = 1027.45140
   !!  a   = 0.179646281
   !!  b   = 0.765555368
   !!  r0  = 4.59763035
   !!  r   = 1022.85377
   !!
   !! REFERENCES:
   !!  Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate 
   !!   polynomial expressions for the density and specific volume of 
   !!   seawater using the TEOS-10 standard. Ocean Modelling.
   !!
   !! polyTEOS, distributed on GitHub: http://github.com/fabien-roquet/polyTEOS
   !! F. Roquet 2016
   !! GNU General Public License

   IMPLICIT NONE
   PRIVATE

   !                   !! * Interface
   PUBLIC   polyTEOS10_insitu     ! Compute in-situ density
   PUBLIC   polyTEOS10_rab        ! Compute thermal expansion and haline contraction
   PUBLIC   polyTEOS10_pt_from_ct ! Compute potential temperature from conservative temperature
   PUBLIC   polyTEOS10_init       ! Initialize constants

   ! TEOS10 parameters
   REAL(kind=8) ::   r1_S0, r1_T0, r1_P0, rdeltaS
   
   ! EOS parameters
   REAL(kind=8) ::   EOS000 , EOS100 , EOS200 , EOS300 , EOS400 , EOS500 , EOS600
   REAL(kind=8) ::   EOS010 , EOS110 , EOS210 , EOS310 , EOS410 , EOS510
   REAL(kind=8) ::   EOS020 , EOS120 , EOS220 , EOS320 , EOS420
   REAL(kind=8) ::   EOS030 , EOS130 , EOS230 , EOS330
   REAL(kind=8) ::   EOS040 , EOS140 , EOS240
   REAL(kind=8) ::   EOS050 , EOS150
   REAL(kind=8) ::   EOS060
   REAL(kind=8) ::   EOS001 , EOS101 , EOS201 , EOS301 , EOS401
   REAL(kind=8) ::   EOS011 , EOS111 , EOS211 , EOS311
   REAL(kind=8) ::   EOS021 , EOS121 , EOS221
   REAL(kind=8) ::   EOS031 , EOS131
   REAL(kind=8) ::   EOS041
   REAL(kind=8) ::   EOS002 , EOS102 , EOS202
   REAL(kind=8) ::   EOS012 , EOS112
   REAL(kind=8) ::   EOS022
   REAL(kind=8) ::   EOS003 , EOS103
   REAL(kind=8) ::   EOS013 
   
   ! ALPHA parameters
   REAL(kind=8) ::   ALP000 , ALP100 , ALP200 , ALP300 , ALP400 , ALP500
   REAL(kind=8) ::   ALP010 , ALP110 , ALP210 , ALP310 , ALP410
   REAL(kind=8) ::   ALP020 , ALP120 , ALP220 , ALP320
   REAL(kind=8) ::   ALP030 , ALP130 , ALP230
   REAL(kind=8) ::   ALP040 , ALP140
   REAL(kind=8) ::   ALP050
   REAL(kind=8) ::   ALP001 , ALP101 , ALP201 , ALP301
   REAL(kind=8) ::   ALP011 , ALP111 , ALP211
   REAL(kind=8) ::   ALP021 , ALP121
   REAL(kind=8) ::   ALP031
   REAL(kind=8) ::   ALP002 , ALP102
   REAL(kind=8) ::   ALP012
   REAL(kind=8) ::   ALP003
   
   ! BETA parameters
   REAL(kind=8) ::   BET000 , BET100 , BET200 , BET300 , BET400 , BET500
   REAL(kind=8) ::   BET010 , BET110 , BET210 , BET310 , BET410
   REAL(kind=8) ::   BET020 , BET120 , BET220 , BET320
   REAL(kind=8) ::   BET030 , BET130 , BET230
   REAL(kind=8) ::   BET040 , BET140
   REAL(kind=8) ::   BET050
   REAL(kind=8) ::   BET001 , BET101 , BET201 , BET301
   REAL(kind=8) ::   BET011 , BET111 , BET211
   REAL(kind=8) ::   BET021 , BET121
   REAL(kind=8) ::   BET031
   REAL(kind=8) ::   BET002 , BET102
   REAL(kind=8) ::   BET012
   REAL(kind=8) ::   BET003

   ! R parameters
   REAL(kind=8) ::   R00 , R01 , R02 , R03 , R04 , R05

   !!----------------------------------------------------------------------
   !! GNU General Public License
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE polyTEOS10_insitu( sa, ct, pres, rho )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE polyTEOS10_insitu  ***
      !!
      !! ** Purpose :   Compute the in situ density from
      !!       temperature and salinity using the polyTEOS10-55t equation of state
      !!
      !! ** Method  :   eos_insitu = rho(sa, ct, pres)
      !!         with   rho    in situ density              kg/m^3
      !!                ct     Conservative Temperature     Celsius
      !!                sa     Absolute Salinity            g/kg
      !!                pres   pressure                     dbar
      !!
      !!         Check value: rho = 1027.45140 kg/m^3 for p=1000 dbar, ct=10 Celcius, sa=30 g/kg
      !!
      !! ** Action  :   compute rho , the in situ density (kg/m^3)
      !!
      !! References :   Roquet et al, Ocean Modelling (2014)
      !!                TEOS-10 Manual, 2010
      !!----------------------------------------------------------------------
      REAL(kind=8), INTENT(in   ) ::   ct             ! conservative temperature  [Celcius]
      REAL(kind=8), INTENT(in   ) ::   sa             ! absolute salinity         [g/kg]
      REAL(kind=8), INTENT(in   ) ::   pres           ! pressure                  [dbar]
      REAL(kind=8), INTENT(  out) ::   rho            ! in situ density           [kg/m^3]
      !
      REAL(kind=8) ::   zt , zh , zs , zr0            ! local scalars
      REAL(kind=8) ::   zn , zn0, zn1, zn2, zn3       !   -      -
      !!----------------------------------------------------------------------
      !
      zp  = pres * r1_P0                          ! pressure
      zt  = ct   * r1_T0                          ! temperature
      zs  = SQRT( ABS( sa + rdeltaS ) * r1_S0 )   ! square root salinity
      !
      zn3 = EOS013*zt   &
         &   + EOS103*zs+EOS003
      !
      zn2 = (EOS022*zt   &
         &   + EOS112*zs+EOS012)*zt   &
         &   + (EOS202*zs+EOS102)*zs+EOS002
      !
      zn1 = (((EOS041*zt   &
         &   + EOS131*zs+EOS031)*zt   &
         &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
         &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
         &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
      !
      zn0 = (((((EOS060*zt   &
         &   + EOS150*zs+EOS050)*zt   &
         &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
         &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
         &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
         &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
         &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
      !
      zn  = ( ( zn3 * zp + zn2 ) * zp + zn1 ) * zp + zn0
      !
      zr0 = (((((R05 ∗ zp+R04) ∗ zp+R03 ) ∗ zp+R02 ) ∗ zp+R01) ∗ zp+R00) ∗ zp
      !
      rho =  zn + zr0                             ! density
      !
   END SUBROUTINE polyTEOS10_insitu


   SUBROUTINE polyTEOS10_rab( sa, ct, pres, a, b )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE polyTEOS10_rab  ***
      !!
      !! ** Purpose :   Calculates thermal/haline expansion at T-points
      !!
      !! ** Method  :   calculates alpha and beta at T-points
      !!                           a = - d(rho) / d(CT)
      !!                           b =   d(rho) / d(SA)
      !!         Check value: for p=1000 dbar, ct=10 Celcius, sa=30 g/kg
      !!                      a = 0.179646281 kg m−3 K−1
      !!                      b = 0.765555368 kg m−3 (g/kg)−1 .
      !!
      !! ** Action  : - a     : thermal expansion
      !! ** Action  : - b     : haline contraction
      !!----------------------------------------------------------------------
      REAL(kind=8), INTENT(in   ) ::   ct             ! conservative temperature  [Celcius]
      REAL(kind=8), INTENT(in   ) ::   sa             ! absolute salinity         [g/kg]
      REAL(kind=8), INTENT(in   ) ::   pres           ! pressure                  [dbar]
      REAL(kind=8), INTENT(  out) ::   a              ! thermal expansion         [kg/m^3/K]
      REAL(kind=8), INTENT(  out) ::   b              ! haline contraction        [kg/m^3/(g/kg)]
      !
      REAL(kind=8) ::   zt , zh , zs                  ! local scalars
      REAL(kind=8) ::   zn , zn0, zn1, zn2, zn3       !   -      -
      !!----------------------------------------------------------------------
      !
      zp  = pres * r1_P0                          ! pressure
      zt  = ct   * r1_T0                          ! temperature
      zs  = SQRT( ABS( sa + rdeltaS ) * r1_S0 )   ! square root salinity
      !
      ! alpha
      zn3 = ALP003
      !
      zn2 = ALP012*zt + ALP102*zs+ALP002
      !
      zn1 = ((ALP031*zt   &
         &   + ALP121*zs+ALP021)*zt   &
         &   + (ALP211*zs+ALP111)*zs+ALP011)*zt   &
         &   + ((ALP301*zs+ALP201)*zs+ALP101)*zs+ALP001
      !
      zn0 = ((((ALP050*zt   &
         &   + ALP140*zs+ALP040)*zt   &
         &   + (ALP230*zs+ALP130)*zs+ALP030)*zt   &
         &   + ((ALP320*zs+ALP220)*zs+ALP120)*zs+ALP020)*zt   &
         &   + (((ALP410*zs+ALP310)*zs+ALP210)*zs+ALP110)*zs+ALP010)*zt   &
         &   + ((((ALP500*zs+ALP400)*zs+ALP300)*zs+ALP200)*zs+ALP100)*zs+ALP000
      !
      zn  = ( ( zn3 * zp + zn2 ) * zp + zn1 ) * zp + zn0
      !
      a   = zn
      !
      ! beta
      zn3 = BET003
      !
      zn2 = BET012*zt + BET102*zs+BET002
      !
      zn1 = ((BET031*zt   &
         &   + BET121*zs+BET021)*zt   &
         &   + (BET211*zs+BET111)*zs+BET011)*zt   &
         &   + ((BET301*zs+BET201)*zs+BET101)*zs+BET001
         !
      zn0 = ((((BET050*zt   &
         &   + BET140*zs+BET040)*zt   &
         &   + (BET230*zs+BET130)*zs+BET030)*zt   &
         &   + ((BET320*zs+BET220)*zs+BET120)*zs+BET020)*zt   &
         &   + (((BET410*zs+BET310)*zs+BET210)*zs+BET110)*zs+BET010)*zt   &
         &   + ((((BET500*zs+BET400)*zs+BET300)*zs+BET200)*zs+BET100)*zs+BET000
      !
      zn  = ( ( zn3 * zp + zn2 ) * zp + zn1 ) * zp + zn0
      !
      b = zn / zs
      !
   END SUBROUTINE polyTEOS10_rab


   SUBROUTINE polyTEOS10_pt_from_ct( ctmp, psal, ptmp )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE polyTEOS10_pt_from_ct  ***
      !!
      !! ** Purpose :   Compute pot.temp. from cons. temp. [Celcius]
      !!
      !! ** Method  :   rational approximation (5/3th order) of TEOS-10 algorithm
      !!       checkvalue: pt=20.02391895 Celsius for sa=35.7g/kg, ct=20degC
      !!
      !! Reference  :   TEOS-10, UNESCO
      !!                Rational approximation to TEOS10 algorithm (rms error on WOA13 values: 4.0e-5 degC)
      !!----------------------------------------------------------------------
      REAL(kind=8), INTENT(in   ) ::   ctmp   ! Cons. Temp [Celcius]
      REAL(kind=8), INTENT(in   ) ::   psal   ! absolute salinity   [g/kg]
      REAL(kind=8), INTENT(  out) ::   ptmp   ! potential temperature [Celcius]
      !
      REAL(kind=8) ::   zt , zs, zn , zd      ! local scalars
      !!----------------------------------------------------------------------
      !
      zt  = ctmp / 40.
      zs  = SQRT( ABS( psal + 5. ) * 0.875 / 35.16504 )
      !
      zn = ((((-2.1385727895e-01*zt   &
         &   - 2.7674419971e-01*zs+1.0728094330)*zt   &
         &   + (2.6366564313*zs+3.3546960647)*zs-7.8012209473)*zt   &
         &   + ((1.8835586562*zs+7.3949191679)*zs-3.3937395875)*zs-5.6414948432)*zt   &
         &   + (((3.5737370589*zs-1.5512427389e+01)*zs+2.4625741105e+01)*zs   &
         &      +1.9912291000e+01)*zs-3.2191146312e+01)*zt   &
         &   + ((((5.7153204649e-01*zs-3.0943149543)*zs+9.3052495181)*zs   &
         &      -9.4528934807)*zs+3.1066408996)*zs-4.3504021262e-01
      !
      zd = (2.0035003456*zt   &
         &   -3.4570358592e-01*zs+5.6471810638)*zt   &
         &   + (1.5393993508*zs-6.9394762624)*zs+1.2750522650e+01
      !
      ptmp = ctmp + zn / zd
      !
   END FUNCTION polyTEOS10_pt_from_ct


   SUBROUTINE polyTEOS10_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE polyTEOS10_init  ***
      !!
      !! ** Purpose :   initializations for the equation of state
      !!
      !! ** Method  :   Read the namelist nameos and control the parameters
      !!----------------------------------------------------------------------
      !
      rdeltaS = 32.
      r1_S0  = 0.875/35.16504
      r1_T0  = 1./40.
      r1_P0  = 1.e-4
      !
      R00 = 4.6494977072e+01
      R01 = −5.2099962525
      R02 = 2.2601900708e−01
      R03 = 6.4326772569e−02
      R04 = 1.5616995503e−02
      R05 = −1.7243708991e−03
      !
      EOS000 = 8.0189615746e+02
      EOS100 = 8.6672408165e+02
      EOS200 = -1.7864682637e+03
      EOS300 = 2.0375295546e+03
      EOS400 = -1.2849161071e+03
      EOS500 = 4.3227585684e+02
      EOS600 = -6.0579916612e+01
      EOS010 = 2.6010145068e+01
      EOS110 = -6.5281885265e+01
      EOS210 = 8.1770425108e+01
      EOS310 = -5.6888046321e+01
      EOS410 = 1.7681814114e+01
      EOS510 = -1.9193502195
      EOS020 = -3.7074170417e+01
      EOS120 = 6.1548258127e+01
      EOS220 = -6.0362551501e+01
      EOS320 = 2.9130021253e+01
      EOS420 = -5.4723692739
      EOS030 = 2.1661789529e+01
      EOS130 = -3.3449108469e+01
      EOS230 = 1.9717078466e+01
      EOS330 = -3.1742946532
      EOS040 = -8.3627885467
      EOS140 = 1.1311538584e+01
      EOS240 = -5.3563304045
      EOS050 = 5.4048723791e-01
      EOS150 = 4.8169980163e-01
      EOS060 = -1.9083568888e-01
      EOS001 = 1.9681925209e+01
      EOS101 = -4.2549998214e+01
      EOS201 = 5.0774768218e+01
      EOS301 = -3.0938076334e+01
      EOS401 = 6.6051753097
      EOS011 = -1.3336301113e+01
      EOS111 = -4.4870114575
      EOS211 = 5.0042598061
      EOS311 = -6.5399043664e-01
      EOS021 = 6.7080479603
      EOS121 = 3.5063081279
      EOS221 = -1.8795372996
      EOS031 = -2.4649669534
      EOS131 = -5.5077101279e-01
      EOS041 = 5.5927935970e-01
      EOS002 = 2.0660924175
      EOS102 = -4.9527603989
      EOS202 = 2.5019633244
      EOS012 = 2.0564311499
      EOS112 = -2.1311365518e-01
      EOS022 = -1.2419983026
      EOS003 = -2.3342758797e-02
      EOS103 = -1.8507636718e-02
      EOS013 = 3.7969820455e-01
      !
      ALP000 = -6.5025362670e-01
      ALP100 = 1.6320471316
      ALP200 = -2.0442606277
      ALP300 = 1.4222011580
      ALP400 = -4.4204535284e-01
      ALP500 = 4.7983755487e-02
      ALP010 = 1.8537085209
      ALP110 = -3.0774129064
      ALP210 = 3.0181275751
      ALP310 = -1.4565010626
      ALP410 = 2.7361846370e-01
      ALP020 = -1.6246342147
      ALP120 = 2.5086831352
      ALP220 = -1.4787808849
      ALP320 = 2.3807209899e-01
      ALP030 = 8.3627885467e-01
      ALP130 = -1.1311538584
      ALP230 = 5.3563304045e-01
      ALP040 = -6.7560904739e-02
      ALP140 = -6.0212475204e-02
      ALP050 = 2.8625353333e-02
      ALP001 = 3.3340752782e-01
      ALP101 = 1.1217528644e-01
      ALP201 = -1.2510649515e-01
      ALP301 = 1.6349760916e-02
      ALP011 = -3.3540239802e-01
      ALP111 = -1.7531540640e-01
      ALP211 = 9.3976864981e-02
      ALP021 = 1.8487252150e-01
      ALP121 = 4.1307825959e-02
      ALP031 = -5.5927935970e-02
      ALP002 = -5.1410778748e-02
      ALP102 = 5.3278413794e-03
      ALP012 = 6.2099915132e-02
      ALP003 = -9.4924551138e-03
      !
      BET000 = 1.0783203594e+01
      BET100 = -4.4452095908e+01
      BET200 = 7.6048755820e+01
      BET300 = -6.3944280668e+01
      BET400 = 2.6890441098e+01
      BET500 = -4.5221697773
      BET010 = -8.1219372432e-01
      BET110 = 2.0346663041
      BET210 = -2.1232895170
      BET310 = 8.7994140485e-01
      BET410 = -1.1939638360e-01
      BET020 = 7.6574242289e-01
      BET120 = -1.5019813020
      BET220 = 1.0872489522
      BET320 = -2.7233429080e-01
      BET030 = -4.1615152308e-01
      BET130 = 4.9061350869e-01
      BET230 = -1.1847737788e-01
      BET040 = 1.4073062708e-01
      BET140 = -1.3327978879e-01
      BET050 = 5.9929880134e-03
      BET001 = -5.2937873009e-01
      BET101 = 1.2634116779
      BET201 = -1.1547328025
      BET301 = 3.2870876279e-01
      BET011 = -5.5824407214e-02
      BET111 = 1.2451933313e-01
      BET211 = -2.4409539932e-02
      BET021 = 4.3623149752e-02
      BET121 = -4.6767901790e-02
      BET031 = -6.8523260060e-03
      BET002 = -6.1618945251e-02
      BET102 = 6.2255521644e-02
      BET012 = -2.6514181169e-03
      BET003 = -2.3025968587e-04
      !
   END SUBROUTINE polyTEOS10_init

   !!======================================================================
END MODULE polyTEOS10_bsq
