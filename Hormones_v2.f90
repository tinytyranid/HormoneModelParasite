program Hormones_v2

  !Module containing all parameter values
  !!NOTE!! New parameter values have to be added to this module
  !In trunk the module is located in the file named parameters.f90
  use parameter_values

  !Module containing all functions as well as some parameters (like pi)
  ! used for calculations
  use functions

  !Module for printing results to binary files
  !In trunk the module is located in the file named binary_write.f90
  use binary_IO

  implicit none

  !VARIABLES
  integer t, L, R, Th, Gh, Ox                               !Variables used for loops
  real(kind=RP)  RRV                                        !Residual Reproductive Value; future fitness
  real(kind=RP)  reserves, reserves_new                     !Variables for reserves [J]
  real(kind=RP)  reserves_max, reserves_max_new             !Variables for max reserves [J]
  real(kind=RP)  surplus_before_growth                      !Variables for surplus
  real(kind=RP)  weight, weight_somatic, weight_somatic_new !Variables for weight [kg]
  real(kind=RP)  length, length_new                         !Variables for length [cm]
  real(kind=RP)  chance_of_max_length                       !Chance of reaching max length
  real(kind=RP)  growth                                     !Variables for growth [kg]
  real(kind=RP)  growth_cost                                !Variables for costs [J]
  real(kind=RP)  THF, GHF, OXF                              !Hormone function levels
  real(kind=RP)  intake                                     !Intake [J timestep -1]
  real(kind=RP)  target_intake                              !Target intake set by OXF, in multiples of SMR_std
  real(kind=RP)  foraging_required                          !Foraging intensity, in multiples of SMR_std
  real(kind=RP)  foraging_cost                              !Variables for costs [J timestep -1]
  real(kind=RP)  parasite_cost                              !Parasite cost [J timestep-1]
  real(kind=RP)  conversion_cost_via_reserves               !Cost of conversion to reserves, and from reserves to metabolism [J]
  real(kind=RP)  conversion_cost_to_growth                  !Cost of conversion due to growth from reserves [J]
  real(kind=RP)  L_new, R_new, L_cat, R_cat                 !L and R category for new length and reserves, continuous
  real(kind=RP)  dL, dR, dE                                 !Integer part of new length and reserves and environmental category...
  integer        intL, intR, intE                           !...and decimal part

  !Variables for environment
  real(kind=RP), dimension(1:E_categories)               :: ValueOfE     !Multiplicative effect of food environment on intake [-]
  real(kind=RP), dimension(1:E_categories)               :: ParasiteInE  !Parasite cost in the different environments [-]
  real(kind=RP), dimension(1:E_categories,1:E_categories):: ProbOfE      !Probability of autocorrelated environment in next timestep [-]
  real(kind=RP), dimension(1:E_categories)               :: ProbAllE     !Probability of being in a certain environment
  real(kind=RP), dimension(1:size(ProbAllE)+1)           :: CumuProbAllE !Array containg the cumulative values of ProbAllE
  real(kind=RP), dimension(1:3)                          :: phiVE        !Array containing phi, ValueOfE and E_real, and is returned from the autocorrelatedE function

  real(kind=RP)  ActualMeanEnvironment
  real(kind=RP)  E_mean, E_real, nextE_real, sumProbOfE
  integer        E_cat, nextE_cat, sumE_cat
  integer        counter
  real(kind=RP)  ValueOfE_real    !Continuous value of E in forward simulation
  real(kind=RP)  ParasiteInE_real !Continuous parasite proportion of intake [-]
  

  !Variables for metabolic rate
  real(kind=RP)  SMR_std, SMR_somatic_std, SMR_exp_CJ, SMR_coeff_CJ, SMR_coeff, SMR_THF, SDA

  !Variables for Oxygen
  real(kind=RP)  O2max_std, O2max_THF, O2_used

  !Variables for mortality
  real(kind=RP)  M_size, M_foraging, survival, M_O2

  !MATRICES
  integer,        dimension(1:3, 1:(E_categories),  1:(R_categories),  1:(L_categories),  1:t_max)           :: strategy
  real(kind=RP),  dimension(     1:(E_categories),  1:(R_categories),  1:(L_categories),  1:t_max)           :: fitness
  real(kind=RP),  dimension(     1:(Ox_categories), 1:(Gh_categories),  1:(Th_categories),1:(E_categories))  :: fitnessconsequences

  !Variables / matrix: forward simulation
  integer       :: n          !Individual
  real(kind=RP) :: end_status !Status at the end of run: Did the individual die, reach max length or run out of time?
  real(kind=RP), dimension(1:35, 1:t_max+1, 1:n_max) :: ind !Individual matrix
  real(kind=RP), dimension(1:n_max, 1:t_max)         :: environment !Array containing how the individuals move through the environment in the forward part of the model  


  !Array containing the number of individuals that survived, was to slow or died during the forwards part of the model
  !Will not be saved as a binary file but will be printed in the end of the model run
  integer, dimension(1:3) :: survarray !surarray(1) = individuals that died, survarray(2) = individuals that were to slow, survarray(3) = survivers
  
  !Timer and declaration for runtime
  real(kind=RP) :: start, finish

  call cpu_time(start)
  
  survarray(:) = 0 !Filling survarray with 0 



  !CALCULATIONS
  !SMR
  !Clarke & Johnston 1999 - General teleost fish
  SMR_exp_CJ = 0.80_RP
  SMR_coeff_CJ = exp(-5.43_RP)                                         !mmol O2 g-1 h-1 Originally -5.43
  SMR_coeff_CJ = SMR_coeff_CJ * 434._RP * 24._RP                       !J g-1 d-1  Conversion from paper, p 895 %Units conversion
  SMR_coeff_CJ = SMR_coeff_CJ * (0.001_RP**(-SMR_exp_CJ)) * t_duration !J kg-1 timestep-1 %Units conversion
  SMR_coeff = SMR_coeff_CJ * (SMR_coeff_weight**(SMR_exp_CJ-SMR_exp))  !Correction according to new exponent, converted for a 3-kg fish




  !Setting fixed seed for random calculations
  call fixed_seed(random_seed_nr)

  !If forwardOnly is set to "n" we run both environment, backward and forward
  forwardONLY_IF: if (forwardOnly == "n") then

    !ENVIRONMENTAL STOCHASTICITY

    !Initializing values
    E_mean = 1._RP+real((E_categories-1),RP)/2._RP !Mean environment
    E_real = E_mean                     !Environment in real
    E_cat = nint(E_real)                !Environment in categories
    ProbOfE(:,:) = 0._RP                !Probability of autocorrelated environment in next timestep
    ProbAllE(:) = 0._RP
    sumE_cat = 0
    
    do counter = 1, maxCounter
      phiVE = autocorrelatedE(E_real, AutoCorr, StochScale, ValueOfE_min, ValueOfE_max, E_categories)
      nextE_real = phiVE(3)
      nextE_cat  = nint(nextE_real)
      !Count of how many times nextE_cat occurs after one another (next cat compared to current cat)
      ProbOfE(nextE_cat,E_cat) = ProbOfE(nextE_cat,E_cat) + 1._RP
      sumE_cat = sumE_cat + nextE_cat

      !Update current for next timestep: what is current environment
      E_real = nextE_real
      E_cat = nextE_cat

    end do

    ActualMeanEnvironment = real(sumE_cat,RP)/real(maxCounter,RP)

    !Use ProbAllE vector - probability of being in each E category irrespective of preceding category - for
    do E_cat=1, E_categories
        ProbAllE(E_cat) = sum(ProbOfE(E_cat,:)) !For each E_real - sum the number of times that category occurred (irrespective of preceding category)
    end do

    ProbAllE(:) = ProbAllE(:)/sum(ProbAllE(:))  !Probability of being in each E category, E_cat (irrespective of preceding category) !Normalized to sum=1 for probabilities

    !Use pE matrix - probability of next category depending on current category (E_cat*) given by pE(:,E_cat*) - for autocorrelated fitness
    ! calculations dependent on current environment, E_cat*
    do E_cat=1, E_categories
      sumProbOfE = sum(ProbOfE(:,E_cat))
      if (sumProbOfE > 0) then
        ProbOfE(:,E_cat) = ProbOfE(:,E_cat)/sumProbOfE !Probability of the next E_cat being a certain value given current E_cat value
      end if
    end do

    !Used to scale feeding(t) value due to stochasticity of environment
    !!NOTE!! Reusing E_cat from above since we are using environmental categories here as well
    do E_cat = 1, E_categories
      !By using the autocorrelatedE() function with AutoCorr=1. the environment returned is the same as
      ! E_cat (lastE = E_real & lastphi=phi) and ValueOfE returned is therefore also the ValueOfE for E_cat
      phiVE = autocorrelatedE(real(E_cat,RP), 1._RP, StochScale, ValueOfE_min, ValueOfE_max, E_categories)
      ValueOfE(E_cat) = phiVE(2) !Getting ValueOfE for E_cat from the phiVE array 
    end do
    
    !Calculating parasite cost in the different environments
    do E_cat=1, E_categories
      ParasiteInE(E_cat) = (ParasiteInE_min + (((ParasiteInE_max-ParasiteInE_min)/(real(E_categories-1)))*(real(E_cat-1))))
    end do
    
    
    !Printing to binary
    call array_to_binary(ValueOfE,    "ValueOfE.bin"   )
    call array_to_binary(ProbofE,     "ProbOfE.bin"    )
    call array_to_binary(ProbAllE,    "ProbAllE.bin"   )
    call array_to_binary(ParasiteInE, "ParasiteInE.bin")

    print*, "Environmental calculations finished!"




    !!BACKWARDS OPTIMIZATION!!
    !Initialise fitness: 1 at L_max, otherwise 0.
    fitness(:,:,:,:) = 0._RP
    fitness(:,:,L_categories,:) = 1._RP

    TIME_LOOP: do t = t_max-1, 1, -1

      !STATES
      LENGTH_LOOP: do L = 1, L_categories-1                                                               
        length = length_min + (real(L,RP)-1._RP) * length_step                                            !Calculate body length [cm]
        weight_somatic = k_somatic*(length**3._RP)                                                        !Converts length into weight with min reserves [kg]
        !Calculates weight at max reserves and subtracts with weight at min reserves and converts this into energy to get max energy in reserves [J]
        reserves_max = (k_max_reserves*(length**3._RP) - weight_somatic)*reserves_energy_density
        M_size = M_size_coeff*(length**M_size_exp)                                                        !Mortality rate due to size [timestep-1]
        O2max_std = O2max_coeff*(weight_somatic**O2max_exp)                                               !Max aerobic scope
        SMR_somatic_std = SMR_coeff * (weight_somatic**SMR_exp)                                           !Baseline SMR based on weight [J timestep-1]

        RESERVES_LOOP: do R = 1, R_categories                              !Energy reserves (relative to max storing capacity at length) [J]
          reserves = reserves_max * (real(R-1,RP)/real(R_categories-1,RP)) !Calculate energy content in reserves [J]
          weight = weight_somatic + reserves/reserves_energy_density       !Weight based on min reserves + actual reserves converted to [kg]
          SMR_std = SMR_coeff * (weight**SMR_exp)                          !Baseline SMR based on weight [J timestep-1]

          !DECISIONS
          THF_LOOP: do Th = 1, Th_categories
            THF = THF_min + (real(Th,RP)-1._RP)*THF_step                           !Thyroid Hormone concentration in blood [ng ml-1]
            SMR_THF = SMR_std*(1._RP+((THF/THF_max)-0.5_RP)*THF_SMR_effect)        !THF increases SMR up to 50%, baseline SMR at standard THF [J timestep-1]
            O2max_THF = O2max_std*(1._RP+((THF/THF_max)-0.5_RP)*THF_O2max_effect)  !Max aerobic scope based on THF, potential...
                                                                                   !...[the difference between minimum and maximum oxygen consumption rate,...
                                                                                   !...02_(min) - 02_(max)]
            GHF_LOOP: do Gh = 1, Gh_categories
              GHF = GHF_min + (real(Gh,RP)-1._RP)*GHF_step           !IGF-1 concentration in blood [ng ml-1]
              growth = (GHF/GHF_max) * growth_max * weight_somatic   !Growth based on IGF-1 levels [kg]; (hormone concentration in larger blood volume can trigger more growth)
              growth_cost = growth*soma_energy_density               !Growth in J to calculate conversion costs below [J]
              weight_somatic_new = weight_somatic + growth           !New somatic weight based on growth [kg]
              length_new = (weight_somatic_new/k_somatic)**(1._RP/3._RP)   !New length calculated using new somatic weight [cm]
              L_new = ((length_new- length_min)/length_step) + 1._RP       !L category for new length, continuous
              intL = max(1,min(int(L_new),L_categories-1))                 !Integer part of new length category...
              dL = max(0._RP,min(L_new - real(intL,RP),1._RP))             !...and decimal part
              reserves_max_new = (k_max_reserves*(length_new**3._RP) - weight_somatic_new)*reserves_energy_density !Calculating new max reserves [J]

              OXF_LOOP: do Ox = 1, Ox_categories
                OXF = OXF_min + (real(Ox,RP)-1._RP)*OXF_step     !Orexin concentration in serum [pg ml-1]

                !ENVIRONMENT
                ENVIRONMENT_LOOP: do E_cat = 1,E_categories
                  target_intake = (OXF/OXF_max)*OXF_effect_on_intake   !Target intake given OXF, in multiples of SMR_std
                  
                  intake = target_intake * SMR_somatic_std            !Target intake converted to [J timestep-1] Depends on somatic weight because predator-prey interactions are determined by length, not influenced by reserves.
                  foraging_required = target_intake / ValueOfE(E_cat) !Foraging required to get intake in given environment, in multiples of SMR
                  SDA = SDA_coeff*intake                                         !Cost of digestion [J timestep-1]
                  foraging_cost = foraging_cost_coeff*foraging_required*SMR_std  !Cost of foraging for food [J timestep-1]...
                                                                                 !...Foraging is more energetically expensive for fat fish (large reserves) due to hydrodynamics
                  parasite_cost = SMR_somatic_std * ParasiteInE(E_cat)           !Calculating parasite costs based on SMR_somatic_std as this is not dependent on reserves [J timestep-1]
                  surplus_before_growth = intake - SDA - SMR_THF - foraging_cost - parasite_cost  !Surplus energy available for GROWTH, storing, reproduction, etc. [J timestep-1]

                  !From positive surplus, there is conversion loss to intermediate metabolites = reserves. Growth is always taken from intermediate metabolites/reserves.
                  if (surplus_before_growth > 0._RP) then
                     conversion_cost_via_reserves = surplus_before_growth * (1._RP-conversion_efficiency_reserves)   !Conversion cost from intake to reserves [J]
                  !For negative surplus; there is a conversion loss for draining reserves to cover metabolic expenses
                  else
                     conversion_cost_via_reserves =  abs(surplus_before_growth) / conversion_efficiency_reserves * (1._RP-conversion_efficiency_reserves)    !Conversion cost from reserves to metabolism [J]                                             !No positive intake that needs to be converted to reserves
                  end if

                  conversion_cost_to_growth = (growth_cost/conversion_efficiency_growth)*(1._RP-conversion_efficiency_growth)  !Growth cost rescaled to requirement, then only conversion cost is considered [J]
                  reserves_new = min((reserves + surplus_before_growth) - (conversion_cost_via_reserves + growth_cost + conversion_cost_to_growth), reserves_max_new) !New reserves after growth [J]

                  R_new = (reserves_new*real(R_categories-1,RP)/reserves_max_new) + 1._RP !R category for new reserves, continuous
                  intR = max(1,min(int(R_new),R_categories-1))                            !Integer part of new reserves category...
                  dR = max(0._RP,min(R_new - real(intR,RP),1._RP))                        !...and decimal part
                  O2_used = SMR_THF + SDA + foraging_cost + conversion_cost_via_reserves + conversion_cost_to_growth

                  !MORTALITY
                  M_foraging = (M_foraging_fact*foraging_required) + (M_foraging_coeff * (foraging_required**M_foraging_exp))  !Proportional increase in predation due to foraging
                  M_O2 = M_O2_coeff * (O2_used/O2max_THF)**M_O2_exp                        !M_O2_exp must be larger than 1.
                  survival = exp( - M_sizeindependent - M_size - (M_size*M_O2) - (M_size*M_foraging) - M_interaction*(M_size*M_O2*M_foraging))  !Probability of surviving the current timestep

                  !If reserves is less than 0 the survival is set to 0
                  if (reserves_new <= 0._RP) then
                    survival = 0._RP
                  !If reserves fall below starvation_limit*reserves_max the survival approaches 0
                  elseif (reserves_new < starvation_limit*reserves_max_new) then
                    survival = survival * (5._RP*reserves_new/reserves_max_new)
                  end if

                  !Interpolation for fitnessconsequences
                  RRV = 0._RP
                  do nextE_cat = 1, E_categories
                    RRV = RRV + ProbOfE(NextE_cat,E_cat)*( &
                      (1._RP-dR)*(1._RP-dL)*fitness(NextE_cat,intR  ,intL  ,t+1) + &
                      (      dR)*(1._RP-dL)*fitness(NextE_cat,intR+1,intL  ,t+1) + &
                      (1._RP-dR)*(      dL)*fitness(NextE_cat,intR  ,intL+1,t+1) + &
                      (      dR)*(      dL)*fitness(NextE_cat,intR+1,intL+1,t+1))
                  end do !E in next timestep

                  fitnessconsequences(Ox,Gh,Th,E_cat) = survival * RRV

                end do ENVIRONMENT_LOOP

              end do OXF_LOOP

            end do GHF_LOOP

          end do THF_LOOP

          !OPTIMIZE strategy for this state combination
          do E_cat = 1, E_categories
            strategy(1:3,E_cat,R,L,t) = maxloc(fitnessconsequences(:,:,:,E_cat)) !Find hormone strategy that maximizes fitness
            fitness(E_cat,R,L,t) = fitnessconsequences(strategy(1,E_cat,R,L,t), strategy(2,E_cat,R,L,t), strategy(3,E_cat,R,L,t),E_cat) !Store optimal fitness (chance of reaching length of maturity during simulation)

          end do

        end do RESERVES_LOOP

      end do LENGTH_LOOP

      !Printing how many points in the hormone strategies that changes between each loop while running
      print*, "------------ Timestep:", t, "------------"
      do E_cat = 1, E_categories
        if (E_cat==1) print "(4a12)"," Environment", "      Orexin", "      Growth", "     Thyroid" !Header only printed during first timestep
        print "(4i12)", E_cat, sum(abs(strategy(1,E_cat,:,:,t)-strategy(1,E_cat,:,:,t+1))), sum(abs(strategy(2,E_cat,:,:,t)-strategy(2,E_cat,:,:,t+1))), sum(abs(strategy(3,E_cat,:,:,t)-strategy(3,E_cat,:,:,t+1)))
      end do

    end do TIME_LOOP

    !Save fitness, fitnessconsequences, strategy to file
    call array_to_binary(fitness,             "fitness.bin"            )
    call array_to_binary(fitnessconsequences, "fitnessconsequences.bin")
    call array_to_binary(strategy,            "strategy.bin"           )


  !If forwardOnly = y then read the strategy, fitness, ProbAllE and ParasiteInE arrays from binary files
  elseif (forwardOnly == "y") then
    call binary_to_array(strategy,    "strategy.bin"   )
    call binary_to_array(fitness,     "fitness.bin"    )
    call binary_to_array(ProbAllE,    "ProbAllE.bin"   )
    call binary_to_array(ParasiteInE, "ParasiteInE.bin")


  !If forwardOnly is neither y or n then print an error message and terminate program
  else
    print*, "-------------------------- !!!ERROR!!! --------------------------"
    print*, "forwardOnly is neither y or n, please check your parameter file"
    print*, "-------------------------- !!!ERROR!!! --------------------------"
    stop !Terminate code
  end if forwardOnly_IF




  !!CALCULATING THE ENVIRONMENT FOR EVERY INDIVIDUAL FOR EVERY TIME STEP!!
  
  !Calling fixed seed for debugging (makes the environment they start in predictable)
  call fixed_seed(random_seed_nr)

  !Making an array with the cumulative values of ProbAllE
  ! used for initializing the environment
  CumuProbAllE = cumucalc(ProbAllE)
  
  Individ_E_LOOP: do n=1, n_max, 1
    !Initializing the environment of each individual using the cumulative probability values in CumProbAllE
    E_real = real(rcumu(CumuProbAllE),RP)                                                     !Returns environment as a real number without any decimals         
    environment(n,1) = E_real
    
    !Calculate environment for next timestep.
    Time_E_LOOP:do t=2, t_max, 1
      phiVE = autocorrelatedE(E_real, AutoCorr, StochScale, ValueOfE_min, ValueOfE_max, E_categories)   
      E_real = phiVE(3)
      environment(n,t) = E_real
    end do Time_E_LOOP
    
  end do Individ_E_LOOP
  
  
  
  
  !!FORWARD SIMULATION!!
  
  ind(:,:,:) = -1000 !Setting all values in ind array to -1000
  t = 1              !Starting with timestep one to initialize values
  n = 1              !Starting with individual one to initialize values

  !Initial values for each individual
  ind(1,1,:) = length_init          !Initial length [cm]

  !Calculate rest of initial individual characteristics
  ind(4,1,:) = k_somatic * (ind(1,1,:)**3._RP)                              !Weight_somatic [kg]
  ind(25,1,:) = (k_max_reserves * (ind(1,1,:)**3._RP) - ind(4,1,:)) * reserves_energy_density  !Max reserves possible at given size [J]
  ind(3,1,:) = ind(25,1,:) * ((R_categories_init-1._RP) / real(R_categories-1,RP)) !Reserves [J]
  ind(5,1,:) = ind(4,1,:) + ind(3,1,:) / reserves_energy_density            !Weight [kg]
  ind(23,1,:) = 1._RP                                                       !Probability of being alive at the beginning of next timestep

  !Saving initial reserve category in percent to ind array for each individual
  ind(2,1,:) = ((R_categories_init-1._RP) / real(R_categories-1,RP))*100._RP                                                      !Probability of being alive at the beginning of next timestep

  Individ_f_LOOP: do n = 1, n_max, 1
    !Runs over all individuals in population
    length = ind(1,1,n)         !Length [cm]
    weight_somatic = ind(4,1,n) !Somatic weight [kg]
    reserves_max = ind(25,1,n)  !Max reserves possible at given size [J]
    reserves = ind(3,1,n)       !Reserves [J]
    weight = ind(5,1,n)         !Total weight [kg]
    survival = ind(23,1,n)      !The probability of surviving the current timestep

    !Runs over all timesteps
    Time_f_LOOP: do t = 1, t_max, 1
    
      !Finding the environment and the food value of the environment
      E_real = environment(n,t)
      phiVE = autocorrelatedE(E_real, 1._RP, StochScale, ValueOfE_min, ValueOfE_max, E_categories) !Since AutoCorr = 1, lastE=E_real and lastphi=phi
      ValueOfE_real = phiVE(2)  

      !Calculate states / variables for this timestep
      M_size = M_size_coeff * (length**M_size_exp)            !Mortality due to size (per timestep)
      O2max_std = O2max_coeff * (weight_somatic**O2max_exp)   !Max aerobic scope [02_(min) - 02_(max)]
      SMR_somatic_std = SMR_coeff * (weight_somatic**SMR_exp) !Baseline SMR based on somatic weight [J timestep-1]
      SMR_std = SMR_coeff * (weight**SMR_exp)                 !Baseline SMR based on total weight [J timestep-1]

      !Calculate integer and decimal part of length, reserves and environment
      L_cat = ((length - length_min) / length_step) + 1._RP !Find length category from length
      intL = max(1,min(int(L_cat),L_categories-2))          !... integer part
      dL = max(0._RP,min(L_cat - real(intL,RP),1._RP))      !... decimal part

      intE =  max(1,min(int(E_real),E_categories-1)) !Find integer part of environmental category
      dE = max(0.,min(E_real - real(intE,RP),1._RP))  !... and decimal part
      
      R_cat = (reserves * real(R_categories-1,RP) / reserves_max) + 1._RP !Find reserve category from reserves
      intR = max(1,min(int(R_cat),R_categories-1))                        !... integer part
  
      !Options for calculating decimal parts of reserves
      
      !Option 1: Old calculation, no switch and decimal part is continuous
      int_d_opt_IF: if (r_i_d_c==1) then
        dR = max(0._RP,min(R_cat - real(intR,RP),1._RP))
      
      !Option 2: Switch between R_cat 1 and 2, decimal part = 0. if R_cat < 2 if R_cat >= the decimal part is continuous 
      elseif (r_i_d_c==2) then        
        if (R_cat < 2) then
          dR = 0._RP
        else
          dR = max(0._RP,min(R_cat - real(intR,RP),1._RP))
        end if
      
      !Option 3: Decimal part of R_cat = 0.
      elseif (r_i_d_c==3) then
        dR = 0._RP
      
      else
        print*, "-------------------------- !!!ERROR!!! --------------------------"
        print*, "intopt is neither 1, 2 or 3, please check your parameter file"
        print*, "-------------------------- !!!ERROR!!! --------------------------"
        stop !Terminate code
      end if int_d_opt_IF
      
      

      !!Interpolate hormone categories from strategy matrix
      !Individuals with body length corresponding to L_categories-1 interpolate for length and reserves
            OXF = &
                (1._RP - dE) * (1._RP - dR) * (1._RP - dL) * strategy(1, intE    , intR    , intL    , 1) + &
                (1._RP - dE) * (1._RP - dR) * (        dL) * strategy(1, intE    , intR    , intL + 1, 1) + &
                (1._RP - dE) * (        dR) * (1._RP - dL) * strategy(1, intE    , intR + 1, intL    , 1) + &
                (1._RP - dE) * (        dR) * (        dL) * strategy(1, intE    , intR + 1, intL + 1, 1) + &
                (        dE) * (1._RP - dR) * (1._RP - dL) * strategy(1, intE + 1, intR    , intL    , 1) + &
                (        dE) * (1._RP - dR) * (        dL) * strategy(1, intE + 1, intR    , intL + 1, 1) + &
                (        dE) * (        dR) * (1._RP - dL) * strategy(1, intE + 1, intR + 1, intL    , 1) + &
                (        dE) * (        dR) * (        dL) * strategy(1, intE + 1, intR + 1, intL + 1, 1)

            GHF = &
                (1._RP - dE) * (1._RP - dR) * (1._RP - dL) * strategy(2, intE    , intR    , intL    , 1) + &
                (1._RP - dE) * (1._RP - dR) * (        dL) * strategy(2, intE    , intR    , intL + 1, 1) + &
                (1._RP - dE) * (        dR) * (1._RP - dL) * strategy(2, intE    , intR + 1, intL    , 1) + &
                (1._RP - dE) * (        dR) * (        dL) * strategy(2, intE    , intR + 1, intL + 1, 1) + &
                (        dE) * (1._RP - dR) * (1._RP - dL) * strategy(2, intE + 1, intR    , intL    , 1) + &
                (        dE) * (1._RP - dR) * (        dL) * strategy(2, intE + 1, intR    , intL + 1, 1) + &
                (        dE) * (        dR) * (1._RP - dL) * strategy(2, intE + 1, intR + 1, intL    , 1) + &
                (        dE) * (        dR) * (        dL) * strategy(2, intE + 1, intR + 1, intL + 1, 1)

            THF = &
                (1._RP - dE) * (1._RP - dR) * (1._RP - dL) * strategy(3, intE    , intR    , intL    , 1) + &
                (1._RP - dE) * (1._RP - dR) * (        dL) * strategy(3, intE    , intR    , intL + 1, 1) + &
                (1._RP - dE) * (        dR) * (1._RP - dL) * strategy(3, intE    , intR + 1, intL    , 1) + &
                (1._RP - dE) * (        dR) * (        dL) * strategy(3, intE    , intR + 1, intL + 1, 1) + &
                (        dE) * (1._RP - dR) * (1._RP - dL) * strategy(3, intE + 1, intR    , intL    , 1) + &
                (        dE) * (1._RP - dR) * (        dL) * strategy(3, intE + 1, intR    , intL + 1, 1) + &
                (        dE) * (        dR) * (1._RP - dL) * strategy(3, intE + 1, intR + 1, intL    , 1) + &
                (        dE) * (        dR) * (        dL) * strategy(3, intE + 1, intR + 1, intL + 1, 1)

      !Convert from category to hormonvalue concentration
      OXF = OXF_min + (OXF - 1._RP) * OXF_step !Orexin [pg ml-1]
      GHF = GHF_min + (GHF - 1._RP) * GHF_step !IGF concentration [ng ml-1]
      THF = THF_min + (THF - 1._RP) * THF_step !T3 concentration  [ng ml-1]

      !Interpolating parasite cost proportion in environment
      ParasiteInE_real = ((1._RP - dE) * ParasiteInE(intE    )) + ((        dE) * ParasiteInE(intE + 1))  

      !Calculating thyroid dependent variables
      SMR_THF = SMR_std * (1._RP + ((THF / THF_max) - 0.5_RP) * THF_SMR_effect)       !SMR based on thyroid, THF increases SMR up to 50%, baseline SMR at standard THF [J timestep-1]
      O2max_THF = O2max_std * (1._RP + ((THF / THF_max) - 0.5_RP) * THF_O2max_effect) !Max aerobic scope based on thyroid [02_(min) - 02_(max)]

      !Calculate growth hormone - dependent variables
      growth = (GHF / GHF_max) * growth_max * weight_somatic  !Growth [kg]
      growth_cost = growth * soma_energy_density              !Cost of new tissue growth [J]

      !Calculating new somatic weight and length after growth
      weight_somatic = weight_somatic + growth               !New somatic weight [kg]
      length = (weight_somatic / k_somatic) ** (1._RP/3._RP) !New length [cm]
      reserves_max = (k_max_reserves * (length**3._RP) - weight_somatic) * reserves_energy_density !Max reserves possible at given size [J]

      !Calculate orexin - dependent variables
      target_intake = (OXF / OXF_max) * OXF_effect_on_intake   !Foraging intensity given OXF in multiples of SMR_std
                  
      intake = target_intake *  SMR_somatic_std            !Intake in [J timestep-1]
      foraging_required = target_intake / ValueOfE_real    !Foraging required to get intake in given environment in multiples of SMR
      SDA = SDA_coeff * intake                                          !Cost of digestion [J timestep-1]
      foraging_cost = foraging_cost_coeff * foraging_required * SMR_std !Cost of foraging for food [J timestep-1]
      parasite_cost = SMR_somatic_std * ParasiteInE_real                 !Calculating parasite costs based on SMR_somatic_std as this is not dependent on reserves [J timestep-1]
      surplus_before_growth = intake - SDA - SMR_THF - foraging_cost - parasite_cost  !Surplus energy available for GROWTH, storing, reproduction, etc. [J timestep-1]

      !From positive surplus, there is conversion loss to intermediate metabolites = reserves. Growth is always taken from intermediate metabolites/reserves.
      if (surplus_before_growth > 0._RP) then
         conversion_cost_via_reserves = surplus_before_growth * (1._RP-conversion_efficiency_reserves) !Conversion cost from intake to reserves [J]
      !For negative surplus; there is a conversion loss for draining reserves to cover metabolic expenses
      else
         conversion_cost_via_reserves =  abs(surplus_before_growth) / conversion_efficiency_reserves * (1._RP-conversion_efficiency_reserves)    !Conversion cost from reserves to metabolism [J]                                                         !No positive intake that needs to be converted to reserves
      end if

      conversion_cost_to_growth = (growth_cost/conversion_efficiency_growth)*(1._RP-conversion_efficiency_growth)  !Growth cost rescaled to requirement, then only conversion cost is considered  [J timestep-1]
      reserves= min(((reserves + surplus_before_growth) - (conversion_cost_via_reserves + growth_cost + conversion_cost_to_growth)), reserves_max) !New reserves after growth [J]

      weight = weight_somatic + reserves / reserves_energy_density  !Total weight [kg]
      O2_used = SMR_THF + SDA + foraging_cost + conversion_cost_via_reserves + conversion_cost_to_growth  !Remember: growth is added as new tissue and should not be included here

      !Mortality
      M_foraging = (M_foraging_fact*foraging_required) + (M_foraging_coeff * (foraging_required**M_foraging_exp))    !Proportional increase in predation due to foraging activity
      M_O2 = M_O2_coeff * (O2_used / O2max_THF)**M_O2_exp                        !Mortality based on size of aerobic scope (=proxy for escapement probability)
      survival = exp( - M_sizeindependent - M_size - (M_size*M_O2) - (M_size*M_foraging) - M_interaction*(M_size*M_O2*M_foraging))!The probability of surviving the current timestep

      !If reserves are less than -100 the individual dies (survival = 0)
      if (reserves <= -100._RP) then
        survival = 0._RP
      !If reserves fall below starvation_limit*reserves_max the survival approaches 0
      elseif (reserves < starvation_limit*reserves_max) then
        survival = survival * (5._RP*max(0._RP,reserves)/reserves_max)
      end if
      
      !If reserves was less than 0, set them to 0,
      !This should help the individuals that make small mistakes in forward
      reserves = max(0._RP,reserves)

      !Calculating chance of reaching max length
      L_cat = ((length - length_min) / length_step) + 1._RP    !Find length category from length..
      intL = max(1,min(int(L_cat),L_categories-1))             !... integer part
      dL = max(0._RP,min(L_cat - real(intL,RP),1._RP))         !... decimal part
      R_cat = (reserves * real(R_categories-1,RP) / reserves_max) + 1._RP !Find reserve category from reserves...
      intR = max(1,min(int(R_cat),R_categories-1))             !... integer part
      dR = max(0._RP,min(R_cat - real(intR,RP),1._RP))         !... decimal part

      !Chance of reaching max length
      chance_of_max_length = &
          (1._RP - dE) * (1._RP-dR)*(1._RP-dL)*fitness(intE    , intR  , intL  , t) + &
          (1._RP - dE) * (1._RP-dR)*(      dL)*fitness(intE    , intR  , intL+1, t) + &
          (1._RP - dE) * (      dR)*(1._RP-dL)*fitness(intE    , intR+1, intL  , t) + &
          (1._RP - dE) * (      dR)*(      dL)*fitness(intE    , intR+1, intL+1, t) + &
          (        dE) * (1._RP-dR)*(1._RP-dL)*fitness(intE + 1, intR  , intL  , t) + &
          (        dE) * (1._RP-dR)*(      dL)*fitness(intE + 1, intR  , intL+1, t) + &
          (        dE) * (      dR)*(1._RP-dL)*fitness(intE + 1, intR+1, intL  , t) + &
          (        dE) * (      dR)*(      dL)*fitness(intE + 1, intR+1, intL+1, t)

      !Did the individual die, survive until max length or was it unable to reach max_length?
      if (survival == 0) then               !If the individual dies (survival = 0.) = 1
        end_status = 1
        survarray(1) = survarray(1)+1       !For printing summary when the model is finished running
      elseif (length >= forward_length_max) then !If the individual reach max length = 3
        end_status = 3
        survarray(3) = survarray(3)+1       
      elseif (t == t_max) then              !If the individual is unable to reach max = 2
        end_status = 2 
        survarray(2) = survarray(2)+1                           
      else
        end_status = 0
      end if
      
      !Writing to ind array
      !NB please note that the numbers are not in ascending order

      !Save variables in ind-matrix as they would be at the end of the timestep
      ind(1,t+1,n) = length                 !Length [cm]
      ind(2,t+1,n) = ((R_cat-1)/real(R_categories-1,RP)) * 100._RP !R_categories transformed to % of reserve fullness
      ind(3,t+1,n) = reserves               !Reserves [J]
      ind(25,t+1,n) = reserves_max          !Max reserves
      ind(4,t+1,n) = weight_somatic         !Structural weight [kg]
      ind(5,t+1,n) = weight                 !Total weight [kg]
      ind(23,t+1,n) = ind(23,t,n)*survival  !Probability of being alive at the beginning of next timestep
      ind(33,t+1,n) = (ind(3,t+1,n) - ind(3,t,n)) !Change in reserves from the start and the end of the timestep [J]
      
      !Information stored from this timestep
      ind(6,t,n) = OXF                      !Orexin [pg ml-1]
      ind(7,t,n) = GHF                      !IGF [ng ml-1]
      ind(8,t,n) = THF                      !T3 [ng ml-1]
      ind(9,t,n) = SMR_std                  !SMR standard
      ind(10,t,n) = SMR_THF                 !SMR adjusted for thyroid hormones
      ind(11,t,n) = growth                  !Growth [kg]
      ind(12,t,n) = target_intake           !Target intake given OXF in multiples of SMR_std
      ind(31,t,n) = foraging_required       !Foraging required given intake and environment in multiples of SMR_std
      ind(13,t,n) = foraging_cost           !Cost of foraging for food [J timestep-1]
      ind(14,t,n) = intake                  !Intake [J timestep-1]
      ind(32,t,n) = surplus_before_growth   !Surplus energy available for growth, etc [J timestep-1]
      ind(15,t,n) = SDA                     !Cost of digestion [J timestep-1]
      ind(16,t,n) = conversion_cost_via_reserves  !Cost of conversion to reserves (->) and from reserves to metabolism (<-) [J], and...
      ind(34,t,n) = conversion_cost_to_growth     !...growth cost rescaled to requirement  [J]
      ind(17,t,n) = O2_used                       !Oxygen used
      ind(18,t,n) = O2max_THF                     !Maximum O2 potential given thyroid hormones
      ind(19,t,n) = chance_of_max_length          !Chance of reaching max length
      ind(24,t,n) = survival                      !The probability of surviving the current timestep
      ind(28,t,n) = E_real                        !Environment
      ind(35,t,n) = parasite_cost                 !Parasite cost in environment [J timestep-1]
      ind(29,1:t_max+1,n) = end_status            !Status at the end of the run
      ind(30,1:t_max+1,n) = t                     !Saving last timestep when the organism was alive

      
      !Mortality from this timestep
      ind(26,t,n) = M_sizeindependent * (365._RP / t_duration)      !Size independent mortality per year
      ind(20,t,n) = M_size *(365._RP / t_duration)                  !Size dependent mortality per year
      ind(21,t,n) = (M_size * M_O2) * (365._RP / t_duration)        !Size dependent and free-scope mortality per year
      ind(22,t,n) = (M_size * M_foraging) * (365._RP / t_duration)  !Size dependent and foraging mortality per year
      ind(27,t,n) = (M_interaction * (M_size * M_O2 * M_foraging)) * (365._RP / t_duration) !Size dependent, free-scope and foraging mortality per year


      !Print data for individual after it reaches max length, dies or at t_max
      if(n==1 .AND. t==1) print "(5a12)", "Individual", "Timestep", "Environment", " Length", " Alive" !Header
      if(end_status /= 0) print "(5i12)", n, t, nint(E_real), nint(length), nint(end_status)

      if (survival == 0) exit Time_f_LOOP                !Exit Time_LOOP if individual is dead, start with new individual
      if (length >= forward_length_max) exit Time_f_LOOP !Exit Time_LOOP if individual has reached max length, start with a new individual

    end do Time_f_LOOP

  end do Individ_f_LOOP

  !Write ind-matrix to binary file
  call array_to_binary(ind, "ind.bin")
  
  print*, "Died", survarray(1), "To slow", survarray(2), "Survived", survarray(3) 
  
  print*, "Individual 1 started in environment nr: ", ind(28, 1, 1)
  print*, "Second environment:                     ", ind(28, 2, 1)
  
  print*, "Individual 2 started in environment nr: ", ind(28, 1, 2)
  print*, "Second environment:                     ", ind(28, 2, 2)
  
  print*, "Individual 3 started in environment nr: ", ind(28, 1, 3)
  print*, "Second environment:                     ", ind(28, 2, 3)
  
  print*, "Individual 4 started in environment nr: ", ind(28, 1, 4)
  print*, "Second environment:                     ", ind(28, 2, 4)
  
  !Prints run time to screen
  call cpu_time(finish)
  print*, "Runtime in minutes :",nint((finish-start)/60)


end program Hormones_v2
