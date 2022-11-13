!*******************************************************************************
! SVN Version ID: $Id: parameters.f90 8130 2019-04-15 08:44:12Z cje012 $
!*******************************************************************************

!Module cointaining all parameter values

module parameter_values

  implicit none

  ! This is the portable declaration of the real type precision kind, in
  ! terms of the number of decimal digits and the maximum range of
  ! the values. This should guarantee replicable results on all systems and
  ! compilers.
  ! For example, the code below makes sure the real values have 15 digits
  ! of precision and can range from 1E-307 to 1E307. This is equivalent to
  ! real kind 8, 64 bit.
  ! and this is for 128 bit reals: selected_real_kind(33, 4931)
  integer, parameter, public :: RP  = selected_real_kind(15,  307)

  !Running only forward (y) or backwards and forward (n)?
  ! NB! If you choose y you need to have the strategy.bin file located in
  ! the folder specified in your folder.f90 file
  character(len=1), parameter, public  :: forwardOnly = "n"

  !Fixed seed for random number generation for forward and backward.
  integer, parameter, public       :: random_seed_nr = 4 


  !Model Resolution
  real(kind=RP), parameter, public :: length_min = 10._RP, length_max = 35._RP, length_step = 0.5_RP    !Max and min body length of fish in model [cm]
  integer,       parameter, public :: L_categories = nint((length_max-length_min)/length_step)+1        !Body length categories
  real(kind=RP), parameter, public :: forward_length_max = 30._RP                                       !Max length for forward used to avoids weird behaviour right before they reach length_max
  integer,       parameter, public :: R_categories = 22                                                 !Energy reserves categories
  real(kind=RP), parameter, public :: THF_min = 0._RP, THF_max = 5._RP, THF_step = 0.25_RP              !T3 concentration in serum [ng ml-1] (Tagawa et al., 1994; Varghese & Oommen, 1999; Varghese et al., 2001; J. Geoffrey Eales, 1988)
  integer,       parameter, public :: Th_categories = nint((THF_max-THF_min)/THF_step)+1                !Thyroid Hormone Function categories
  real(kind=RP), parameter, public :: GHF_min = 0._RP, GHF_max = 200._RP, GHF_step = 1.25_RP            !IGF-1 concentration in serum [ng ml-1] (Cameron et al., 2007; Uchida et al., 2003)
  integer,       parameter, public :: Gh_categories = nint((GHF_max-GHF_min)/GHF_step)+1                !Growth Hormone Function categories
  real(kind=RP), parameter, public :: OXF_min = 0._RP, OXF_max = 2500._RP, OXF_step = 37.5_RP           !Orexin concentration in blood serum [pg ml-1] (from pigs, humans) (Tomasik et al., 2004; Kaminski et al., 2013)
  integer,       parameter, public :: Ox_categories = nint((OXF_max-OXF_min)/OXF_step)+1                !Orexin Hormone Function categories
  integer,       parameter, public :: t_max = 250       !200                                            !Number of timesteps
  real(kind=RP), parameter, public :: t_duration = 7._RP  !7                                            !Days in each timestep
  integer,       parameter, public :: E_categories = 12  !11                                            !Number of food environment categories

  !PARAMETERS
  real(kind=RP), parameter, public :: k_somatic = 0.85e-5_RP, k_max_reserves = 1.2e-5_RP !Minimum (no reserves) and maximum (max reserves) Fulton's condition factor (units kg weight, cm length) (Lambert & Dutil, 1997)
  real(kind=RP), parameter, public :: reserves_energy_density = 5.e6_RP                  !J kg-1 !e is the same as *10^x 
  real(kind=RP), parameter, public :: soma_energy_density     = 4.e6_RP                  !J kg-1 
  real(kind=RP), parameter, public :: SMR_exp = 0.7_RP !0.7                              !SMR scales proportional to W^SMR_Exp (Killen et al., 2007)
  real(kind=RP), parameter, public :: SMR_coeff_weight = 3._RP                           !Weight used for correction of SMR_coeff according to new exponent
  real(kind=RP), parameter, public :: THF_SMR_effect = 0.25_RP                           !Proportion increase in SMR between min and max THF
  real(kind=RP), parameter, public :: THF_O2max_effect = 0.2_RP                          !Proportion increase in SMR at max THF level 
  real(kind=RP), parameter, public :: OXF_effect_on_intake = 8.5_RP                      !Max intake at max OXF in multiples of SMR_std 
  real(kind=RP), parameter, public :: growth_max = 0.04_RP*t_duration                    !Max growth at max GHF, proportion weight per day
!  real(kind=RP), parameter, public :: intake_max_per_kg = 0.08e6_RP*t_duration          !Max intake per day in J kg-1 timestep-1 used Jonsson et al., 2001 fig 2 (p 707) to estimate a value  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(kind=RP), parameter, public :: conversion_efficiency_reserves = 0.85_RP           !Conversion efficiency for intake to reserves and from reserves pay for metabolism and fixed costs (reference???)
  real(kind=RP), parameter, public :: conversion_efficiency_growth = 0.75_RP             !Conversion efficiency for energy taken to reserves to growth 
  real(kind=RP), parameter, public :: SDA_coeff = 0.15_RP                                !Cost of digestion, proportion of intake
  real(kind=RP), parameter, public :: foraging_cost_coeff = 0.2_RP                       !Energetic cost of foraging behaviour, in SMR units per behaviour unit 
  real(kind=RP), parameter, public :: starvation_limit = 0.01_RP                         !When reserves fall below this proportion of reserves*max the survival approaches 0


  !Parameters used for environment
  real(kind=RP), parameter, public :: AutoCorr = 0.80_RP       !The strength of the autocorrelation if AutoCorr = 0 there is no autocorrelation.
                                                               !If 0 < AC < 1 the noise is positively autocorrelated and if -1 < AC < 0 it is negativly correlated. (Ripa & Lundberg, 1996)
  real(kind=RP), parameter, public :: StochScale = 0.35_RP     !Number of standard deviations category units that correspond to E_cat=1 (and max)
  real(kind=RP), parameter, public :: ValueOfE_min = 1._RP - 0.64_RP, ValueOfE_max = 1._RP + 0.4_RP  !The difference between the worst and the best environment
                                                                                                      ! Also defines the shape of the distribution of the environments
  real(kind=RP), parameter, public :: ParasiteInE_min = 0._RP, ParasiteInE_max = ParasiteInE_min      !The difference between the worst and the best environment in terms of parasite costs
                                                                                                      ! Parasite cost is a proportion of intake ParasiteCost = intake*ParasiteInE[E] and it is linear
  integer,       parameter, public :: maxCounter = 1000000      !Counter used for probability calculations of the environment
  
  !AMR - Adapted from Claireaux et al. 2000   ----- AMR or MMR ???
  real(kind=RP), parameter, public :: Temperature = 5._RP                                                                    !degrees C
  !Measured aerobic scope at given temperature (from Claireaux et al. 2000, cod, eq. 9, for 100% O2-saturation)
  real(kind=RP), parameter, public :: O2max_coeff = ((17.29_RP*Temperature**(-0.015_RP*Temperature+1.062_RP)+30.01_RP) *   & !mgO2 h-1 kg-1 in paper converted to J timestep-1 kg-1
                                      (1._RP-exp(-0.035_RP*100+0.34_RP)) * 336.3_RP * t_duration)
  real(kind=RP), parameter, public :: O2max_exp = SMR_exp


  !Mortality
  real(kind=RP), parameter, public :: M_size_coeff = 1.3_RP*(t_duration/365._RP)
  real(kind=RP), parameter, public :: M_size_exp = -0.75_RP                               !M_size=M_size_coeff*L^M_size_exp
  real(kind=RP), parameter, public :: M_sizeindependent = 0.01_RP*(t_duration/365._RP)    !Fixed background mortality [timestep-1]
  real(kind=RP), parameter, public :: M_foraging_fact = 0._RP                             !Always used together with M_foraging_coeff og M_foraging_exp: M_foraging = (fact * forag) + (coeff*(forag^exp))
  real(kind=RP), parameter, public :: M_foraging_coeff = 0.03_RP, M_foraging_exp = 3._RP  !Mortality cost of foraging behaviour !Always used together
  real(kind=RP), parameter, public :: M_O2_exp = 2.7_RP,    &                             !M_O2_exp must be larger than 1 
                                      M_O2_coeff = 1.3_RP                                 !Free-scope mortality !Always used together
  real(kind=RP), parameter, public :: M_interaction  = 1.2_RP                             !Scaling of interaction effects of mortality components)

  !Parameters for FORWARD
  real(kind=RP), parameter, public :: R_categories_init = real(R_categories,RP) * 0.15_RP !Initial reserve category of individual in timestep 1 in forward simulation
  real(kind=RP), parameter, public :: length_init = 10._RP       !Initial length of all individuals [cm]
  integer,       parameter, public :: n_max = 10000              !Max number of individuals in population
  integer,       parameter, public :: r_i_d_c = 1                !Calculation type for integer and decimal parts of R_cat 
                                                                 !... 1. Starts with reserve level 1 as the lowest value, interpolates strategies for all reserve_categories
                                                                 !... 2. For R_cat < 2 individuals only use 1.0 in interpolation of strategy, for R_cat >= 2, it is the same as in option 1
                                                                 !... 3. Individuals only use R_categories as a whole number; 1.0, 2.0, 3.0 etc.


end module parameter_values
