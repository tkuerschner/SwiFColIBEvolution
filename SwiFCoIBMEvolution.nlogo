;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  SwiFCoIBM_Evolution - Individual-Based Swine Fever Community Model with explicit movement, dynamic landscapes and pathogen evolution (based on model versions published in Kramer-Schadt et al. 2009, Lange et al. 2012, Scherer et al. 2020             ;;
;;  and Kürschner et al. 2021)                                                                                                                                                                                                                               ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  The model is composed of two major components:                                                                                                                                                                                                           ;;
;;    (1) a wild boar demography model considering seasonal reproduction, herd splitting (dispersal) and mortality and                                                                                                                                       ;;
;;    (2) a Classical Swine Fever (CSF) virus model operating on the emerging wild boar population.                                                                                                                                                          ;;
;;                                                                                                                                                                                                                                                           ;;
;;  Wild boar population density and structure are affected by the disease via virus-induced mortality and litter size depression.                                                                                                                           ;;
;;  The crucial model entity is the wil boar individual, characterised by age in weeks (one week represents the approx. CSF incubation time; Artois et al. 2002, Moenning et al. 2003),                                                                      ;;
;;  resulting in age classes piglets (age <= 8 months), yearlings (subadults; 8 months < age <= 2 years) and adults (age > 2 years).                                                                                                                         ;;
;;  Each host has a location, which denotes its home range cell and the individual's family group (HerdID).                                                                                                                                                  ;;
;;                                                                                                                                                                                                                                                           ;;
;;  Movement: (in select model versions) Adult female movement is restrcited to the 'activity zone' which equals their home range cell (staying strategy). Movement is similar similar for individuals of the sounder (social unit organised around          ;;
;;            adult females and their offspring). Subadult females split in specified weeks of the year from their herds searching for an unoccupied home range cell (see Procedure "HerdSplit").                                                            ;;
;;            Subadult males disperse from natal family groups in subgroups at an age the males begin to seek mates (~5 months), joining other groups in the surrounding area.                                                                               ;;
;;            Adult and elderly males tend to range solitary outside of their activity zone (ranging strategy).                                                                                                                                              ;;
;;                                                                                                                                                                                                                                                           ;;
;;  Landscape: The model landscape is represented by a grid of 4 km^2 square cells which each encompass a wild boar family group's home range (Leaper et al. 1999).                                                                                          ;;
;;             Each habitat cell is characterised by a breeding capacity (Quality), denoting habitat quality by the number of female boars that are allowed to breed, representing density regulation in the model.                                          ;;
;;             The landscape can be either static or dynamic with autocorrelated noise (red noise) or random noise (white noise).                                                                                                                            ;;
;;                                                                                                                                                                                                                                                           ;;
;;                                                                                                                                                                                                                                                           ;;
;;                                                                                                                                                                                                                                                           ;;
;;                                                                                                                                                                                                                                                           ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  STATE VARIABLES                                                                                                                                                                                                                                          ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; NAME                            TYPE/RANGE     DEFAULT        DESCRIPTION

globals
[
  ;Landscape                  ;   [string]       -              way of landscape generation or import
  Habitat                     ;   [set]          -              habitat cells (cells with Quality > 0)
  Matrix                      ;   [set]          -              matrix cells (cells with Quality = 0)
  ;MeanQuality                ;   [#]            9              maximum value for habitat quality to define range of habitat quality (quality between 0 and MeanQuality)
  ;YearRelease                ;   [#]            6              year in which an randomly individual becomes infected (range of WeekRelease)
  WeekRelease                 ;   [#]            261-312        random week in YearRelease of virus release following an expontential distribution with lambda = 1.5
  ;Herd%                      ;   [%]            100            proportion of habitat patches occupied by herds
  ;ReleaseFactor              ;   [#, >=1]       3              multiplies breeding capacity of each occupied patch to determine number of boars
  MaxAge                      ;   [#]            572 weeks      maximum age in weeks
  ;Longevity                  ;   [#]            11 years       maximum age in years (MaxAge / 52)
  ;AgeBlur                    ;   [#]            6 weeks        stochastic deviation from age class thresholds in weeks (in both directions (age class threshold +/- AgeBlur) following a normal distribution
  SurvivalProb                ;   [%]            0.65           annual survival rate of adults and subadults (mean = 0.65, min = 0.4)
  SurvivalProbPig             ;   [%]            0.5            annual survival rate of piglets (mean = 0.5, min = 0.1)
  currHerdID                  ;   [#]            1              current maximum herd ID
  ;CaseFatality               ;   [%]            0.8            general risk for getting lethally infected by a disease transmission
  ;BetaWithin                 ;   [#, 0-1]       0.0208         transmission coefficient within the herd
  ;BetaBetween                ;   [#, 0-1]       0.00208        transmission coefficient between herds
  ;BetaMove                   ;   [#, 0-1]       -              transmission coefficient during roaming movements
  ;PreNatInf                  ;   [#, 0-1]       0.5            probability of pre-natal infection
  ;FertRedInf                 ;   [#, 0-1]       0.625          reduction in fertility of infected females
  ;mue                        ;   [#]            3 weeks        mean of exponential survival time distribution for lethally infected individuals
  ;T_trans                    ;   [#]            1 week         infectious period of transiently infected individuals
  ;T_immu                     ;   [#]            12 weeks       maximum duration of immunity by maternal antibodies (Depner et al. 2000)
  ;T_anti                     ;   [#]            15 weeks       maximum persistence of maternal antibodies (Depner et al. 2000)
  ;T_latent                   ;   [#]            4 weeks        duration of latent transmission
  ProbRepro                   ;   [#]                           monthly reproduction probability
  week                        ;   [#, 1-52]      0              current week of the year
  ;SimulatedYears             ;   [#]            100            number of simulated years during one run
  ;ProbFem                    ;   [%]            0.5            sex ratio (0: only males; 1: only females)
  NewLeth                     ;   [#]            0              number of individuals which got lethally infected during this timestep
  NewTrans                    ;   [#]            0              number of individuals which got transient infected during this timestep
  ;q                          ;   [#, 0-1]       -              probability to move straight to target cells
  Tag                         ;   [who]          -              random chosen roaming male
  TrackCol                    ;   [#]            -              pen color of tracked roaming individuals
  MNR                         ;   [#]            12             mean of log-normal distribution for normal-distance roaming males
  ;FemDispEvent               ;   [#]            -              counting natal dispersal events of females
  ;MaleDispEvent              ;   [#]            -              counting natal dispersal events of males
  NmaxList                    ;   [list]         []             input data: sums of the number of patches adjacent to each patch in a circle of the same area as the infected/infectious area (see Report_fragmentation)
  ReproList                   ;   [list]         []             input data: cumulative reproduction probabilities for each month
  MeanLethPeriod              ;   [#]            -              mean infectious period of lethal infected individuals in weeks (equals mue)
  MaxLethPeriod               ;   [#]            -              maximum infectious period of lethal infected individuals in weeks (equals mue)
  rep                         ;   [#]            -              number of repetitions in case of running BehaviourSpace experiments
  ;seed                       ;   [#]            -              seed used if seed-setup = "ON"
  DONE                        ;   [boolean]      0              1 if stop conditions are TRUE
  week_overlap                ;                                 overlap
  max_id                      ;                                 id helper
  rednoise_list               ;   [list]         []             external input of red noise values

  ;; OUTPUT
  Time                        ;   [min]          0              counting minutes since starting setup
  InfDist                     ;   [cells]        0              maximum spatial distance of disease transmission
  InfX                        ;   [#]            -              coordinate x of the pathogen-introduced cell
  InfY                        ;   [#]            -              coordinate y of the pathogen-introduced cell
  MaxDist                     ;   [#]            -              week when infection reached other side of the tubus
  WeekLast                    ;   [#]            -              week in which the last infected indiviual was noted
  MeanMoveDist                ;   [#]            -              mean weekly movement distance
  outbreak_group              ;   [#]            -              number of individuals in release cell
  outbreak_roaming            ;   [#]            -              number of adult males in release cell
  outbreak_dens5              ;   [#]            -              number of individuals in the release cell and its 5 neighbors
  outbreak_dens12             ;   [#]            -              number of individuals in a radius of 12 cells from the release cell
  F_infected                  ;   [0-1]          -              fragmentation index for infectious cells (patches which contain shedding individuals)
  F_infectious                ;   [0-1]          -              fragmentation index for infected cells (patches which were infected at least one week during the simulation)
  PopStructList               ;   [list]         []             current age class distribution of all individuals in years
  InfStructList               ;   [list]         []             current age class distribution of infected individuals (EpiStat = esLeth or esTrans) in years
  InfTransStructList          ;   [list]         []             current age class distribution of transient infected individuals (EpiStat = esTrans) in years
  InfLethStructList           ;   [list]         []             current age class distribution of lethally infected individuals (EpiStat = esLeth) in years
  CurrInfPeriodList           ;   [list]         []             current infectious periods of lethally infected individuals (EpiStat = esLeth) in weeks
  LethPeriodList              ;   [list]         []             all infectious periods of lethally infected individuals (EpiStat = esLeth) in weeks
  InfPeriodList               ;   [list]         []             all infectious periods of lethally and transient infected individuals (EpiStat = esLeth and esTrans) in weeks
  controll                    ;   [#]            -              controll helper variable
  qual                        ;   [#]            -              controll helper variable
  o_list                      ;   [list]         []             output for infected per lanscape cell
  cells_in_order              ;   [list]         []             all landscape cells ordered by coordinates
  meanStrain                  ;   [#]            -              mean strain of all infected individuals
  BetaWHigh                   ;   [#]            -              virulence / strain variable
  BetaWLow                    ;   [#]            -              virulence / strain variable
  BetaWMid                    ;   [#]            -              virulence / strain variable
  BetaBMid                    ;   [#]            -              virulence / strain variable
  BetaBLow                    ;   [#]            -              virulence / strain variable
  BetaBHigh                   ;   [#]            -              virulence / strain variable
  testCounterAlltime          ;   [#]            -              virulence / strain variable
  sumOfInfectionEvents        ;   [#]            -              virulence / strain variable
  sumOfEvolutionEvents        ;   [#]            -              virulence / strain variable
  BetaW_Middle                ;   [#]            -              virulence / strain variable
  BetaB_Middle                ;   [#]            -              virulence / strain variable
  survivalTimeList            ;   [#]            -              virulence / strain variable
  transmissionFailures        ;   [#]            -              virulence / strain variable
  warningCounter              ;   [#]            -              virulence / strain variable
  fixedExtTimesList           ;   [#]            -              virulence / strain variable
  meanStrainOutput            ;   [#]            -              virulence / strain variable
  transList                   ;   [list]         -              transmission facor list
  strainList                  ;   [list]         -              list of stains
]

turtles-own
[
  Age                         ;   [#, 0-MaxAge]  -              individual age in weeks
  is_female                   ;   [boolean]      -              individual sex: female or male
  HerdID                      ;   [#]            -              herd membership
  IndCaseFatality             ;   [%]            ind.           individual mortality rate (CaseFatality ^ 2 for adults, sqrt CaseFatality for Juveniles)
  AgeGroup                    ;   [string]       -              age group turtle belongs to: piglet, subadult, adult; differs between sexes
  DemStat                     ;   [string]       dsResident     demographic status: dsResident, dsOffspring, dsDispF, dsDispM_adult, dsDispM_subadult
  EpiStat                     ;   [string]       esSusc         epidemiological status: esSusc, esTrans, esLeth, esImm, esImmMat
  is_shedding                 ;   [boolean]      0              1 if the individual is infected (EpiStat = esTrans or esLeth
  TickInf                     ;   [#]            -              week of infection
  TickLeth                    ;   [#]            -              week of death after infection (for esLeth)
  IndBetaMove                 ;   [#, 0-1]       BetaMove/Dist  individual transmission probabilities for roaming individuals based on movement steps of current week (~time spent per cell)
  Family                      ;   [#]            -              contains state variable [who] of the mother
  ChanceBreeding              ;   [#, 0-1]       -              random number for stochastic effects in monthly reproduction for the given year
  TickBirth                   ;   [#]            -              week of birth
  is_matAB                    ;   [boolean]      0              1 if mother had anti bodies (= was immune)
  is_mother                   ;   [boolean]      0              1 if female had offspring in current year
  is_breeder                  ;   [boolean]      0              1 if female is allowed to breed in current year
  MoveDist                    ;   [#]            -              individual weekly movement distance
  NoGroup                     ;                  -
  float_counter               ;   [#]            -              number of boars without valid home cells
  quality_memory              ;   [#]            -              habitat quality of cells
  strain                      ;   [#]            -               strain
  hasActiveStrain             ;   [#]            -               infected
  initialOutbreak             ;   [boolean]      -               initial outbreak
  transientCheck              ;   [#]            -               transient / leathal infection
  strainToGet                 ;   [#]            -               low, med, high virulence
  personalBetaWithin          ;   [#]            -               individual beta
  personalBetaBetwen          ;   [#]            -               individual beta
  markedForInfection          ;   [boolean]      -               determines infection
  strainVirToGet              ;   [#]            -               determines which strain it got its infection from 1 low 2 mid 3 high
  IndBetaMoveHigh             ;   [#]            -               individual beta
  IndBetaMoveLow              ;   [#]            -               individual beta
  counterGotInfected          ;   [#]            -               count variable
  moveBetaList                ;   [list]         -               list of beta for moving indioviduals

]

patches-own
[
  is_occupied                 ;   [boolean]      0              1 if a patch contains a herd of wild boars
  Quality                     ;   [#]            [0-9]          habitat quality in terms of breeding capacity = number of females alllowed to breed
  is_habitat                  ;   [boolean]      0              true if Quality > 0
  Qual_around                 ;   [#]            [0-9]          summed habitat quality of surrounding cells (to decide if movement is possible)
  Capacity                    ;   [#]            [0-40]         overall capacity = Quality * (mean number of 20 boars per cell / mean Quality) = Quality * (20 / 4.5) = 4.4445 boars per quality unit
  NoInfected                  ;   [#]            0              number of infectious indiviuals of this habitat patch
  NoInfectedMove              ;   [#]            0              number of infectious moving indiviuals of this habitat patch
  InfPress                    ;   [#]            0              infection pressure based on infected indivudals in and around this patch
  is_infectious               ;   [boolean]      0              1 if patch is infected at the current step
  is_infected                 ;   [boolean]      0              1 if patch has been infected during simulation
  inf_count                   ;   [#]            0              keeps track of how many times a patch contained one or more infected individuals
  NoBoars                     ;   [#]            0              number of individuals per cell
  NoRepro                     ;   [#]            0              number of reproductive females (subadult + adult) per cell
  NoMale                      ;   [#]            0              number of  adult males per cell
  NoDispF                     ;   [#]            0              number of female dispersers per cell
  DispCells                   ;   [set]          -              cells in dispersal distance of the focal cell
  W_M
  dynamic                     ;   [boolean]      0              1 for dynamic landscape, 0 for static landscape
  counter_d                   ;   [#]            [0-9]          habitat quality / breeding capacity
  qmw                         ;   [#]            [0-9]          mean habitat quality
  dec                         ;   [#]            -              dynamic helper variable
  inc                         ;   [#]            -              dynamic helper variable
  qmax                        ;   [#]            [0-9]          maximum habitat quality
  qmin                        ;   [#]            [0-9]          minimum habitat quality
  sq                          ;   [#]            -              dynamic helper variable
  controll_mean               ;   [#]            -              buffer helper variable
  counter_mean                ;   [#]            -              buffer helper variable
  id                          ;   [#]            -              id
  infection_output_counter    ;   [#]            -              output helper variable
  InfPressLowVir              ;   [#]            -              grouped infection pressure  low
  InfPressHighVir             ;   [#]            -              grouped infection pressure  high
  NoInfectedM                 ;   [#]            -              counters
  NoInfectedLowVir            ;   [#]            -              counters
  NoInfectedHighVir           ;   [#]            -              counters
  NoInfectedAdj               ;   [#]            -              counters
  AllSus                      ;   [#]            -              counters
  tmpToInfectHere             ;   [#]            -              counters
  NoInfectedMoveH             ;   [#]            -              counters
  NoInfectedMoveL             ;   [#]            -              counters
  toInfect                    ;   [#]            -              counters
  ip1                         ;   [#]            -              strain infection pressure
  ip2                         ;   [#]            -              strain infection pressure
  ip3                         ;   [#]            -              strain infection pressure
  ip4                         ;   [#]            -              strain infection pressure
  ip5                         ;   [#]            -              strain infection pressure
  ip6                         ;   [#]            -              strain infection pressure
  ip7                         ;   [#]            -              strain infection pressure
  ip8                         ;   [#]            -              strain infection pressure
  ip9                         ;   [#]            -              strain infection pressure
  ip10                        ;   [#]            -              strain infection pressure
  ip11                        ;   [#]            -              strain infection pressure
  ip12                        ;   [#]            -              strain infection pressure
  ip13                        ;   [#]            -              strain infection pressure
  ip14                        ;   [#]            -              strain infection pressure
  ip15                        ;   [#]            -              strain infection pressure
  ip16                        ;   [#]            -              strain infection pressure
  ip17                        ;   [#]            -              strain infection pressure
  ip18                        ;   [#]            -              strain infection pressure
  ip19                        ;   [#]            -              strain infection pressure
  ip20                        ;   [#]            -              strain infection pressure
  ipMin                       ;   [#]            -              strain infection pressure
  ipMax                       ;   [#]            -              strain infection pressure
  in1                         ;   [#]            -              number of infected per strain
  in2                         ;   [#]            -              number of infected per strain
  in3                         ;   [#]            -              number of infected per strain
  in4                         ;   [#]            -              number of infected per strain
  in5                         ;   [#]            -              number of infected per strain
  in6                         ;   [#]            -              number of infected per strain
  in7                         ;   [#]            -              number of infected per strain
  in8                         ;   [#]            -              number of infected per strain
  in9                         ;   [#]            -              number of infected per strain
  in10                        ;   [#]            -              number of infected per strain
  in11                        ;   [#]            -              number of infected per strain
  in12                        ;   [#]            -              number of infected per strain
  in13                        ;   [#]            -              number of infected per strain
  in14                        ;   [#]            -              number of infected per strain
  in15                        ;   [#]            -              number of infected per strain
  in16                        ;   [#]            -              number of infected per strain
  in17                        ;   [#]            -              number of infected per strain
  in18                        ;   [#]            -              number of infected per strain
  in19                        ;   [#]            -              number of infected per strain
  in20                        ;   [#]            -              number of infected per strain
  inMin                       ;   [#]            -              number of infected per strain
  inMax                       ;   [#]            -              number of infected per strain
  infPressList                ;   [list]         -              list of infection pressures
  virToGetList                ;   [list]         -              list of transmissable strains
  infectionList               ;   [list]         -              list of infections
  initialOutbreakCell         ;   []             -              starting cell for infection
]

extensions
[
  profiler
  palette
]




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  SETUP + GO PROCEDURES                                                                                                                                                                                                                                    ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; SETUP ______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to Setup

  ca
  reset-ticks
  reset-timer

  set rep 150
  set week_overlap 0

  Init_Parameters
  Init_Landscape
  Init_Population
  set o_list []

end

to debug

  Setup                  ;; set up the model
  profiler:start         ;; start profiling
  repeat 200 [ Go ]      ;; run something you want to measure
  profiler:stop          ;; stop profiling
  print profiler:report  ;; view the results
  profiler:reset         ;; clear the data

end

;; GO _________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to Go

  ifelse ticks mod 52 = 0           [ set week 1 ] [ set week week + 1 ]

  if ticks = WeekRelease - 1        [ Disease_Release ]  ;; - 1 to be synchronous*

  if ticks >= WeekRelease - 1       [
                                     Disease_Transmission
                                     infectIndividuals
                                    ]   ;; - 1 to be synchronous*

  if invasion_analysis = true AND ticks >= WeekRelease - 1
  [
    if week = 33 AND ticks >= (WeekRelease - 1 + (invasion_timeDelay * 52))
    [
    invasion
    ]
  ]

                                      Movement_Roaming

  if week = 17                      [ Movement_Dispersal-Males ]

  if week = 29                      [ Movement_Dispersal-Females ]

                                      Host_Reproduction

                                      Host_Death

                                      Host_Ageing

  if ticks >= WeekRelease - 1       [ Disease_Transition ]

                                      dyn

                                      Update

                                      check_floaters

                                      r_move

                                      resetdyn

  testcounter
  infectionCheck

  tick

end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  INITALIZATION PROCEDURES                                                                                                                                                                                                                                 ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; SET PARAMETERS _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Sets all parameters needed for initalisation:
  ;; Current number of habitat patches and maximum number of habitat patches, thresholds for landscape pattern, case fatality, risk of getting lethally infected, week of virus release


to Init_Parameters
  set testCounterAlltime 0
  if seed-setup = "ON" [ random-seed seed ]
  if seed-setup = "BS" [ random-seed behaviorspace-run-number mod rep ]
  set transmissionFailures 0
  set warningCounter 0

  set MaxAge Longevity * 52
  set week 0
  set WeekRelease ((YearRelease - 1) * 52) + release_week + 1   ; YearRelase = 6 -> "during the sixth year" -> reduce number to year 5; + 1 to be between week 1 and 52
  set InfDist 0
  set MaxDist 0
  set sumOfInfectionEvents 0
  set sumOfEvolutionEvents 0
 if uniform_breeding = FALSE
  [
  set ReproList [ 0.06    ;; January
                  0.16    ;; February
                  0.39    ;; March
                  0.73    ;; April
                  0.8     ;; May
                  0.88    ;; June
                  0.94    ;; July
                  0.97    ;; August
                  1.0     ;; September
                  0.0     ;; October
                  0.0     ;; November
                  0.0 ]   ;; December
  ]


  if uniform_breeding = TRUE
  [
   set ReproList [0.08333333    ;; January
                  0.1666667    ;; February
                  0.25    ;; March
                  0.3333333    ;; April
                  0.4166667     ;; May
                  0.5   ;; June
                  0.5833333   ;; July
                  0.6666667    ;; August
                  0.75     ;; September
                  0.8333333     ;; October
                  0.9166667     ;; November
                  1 ]   ;; December
  ]


  ;; INPUT DATA
  file-open "RadiusAll.txt"
  set NmaxList []
  while [ not file-at-end? ] [ set NmaxList lput file-read NmaxList ]
  file-close


  ;; OUTPUT
  set InfStructList []
  set InfTransStructList []
  set InfLethStructList []
  set CurrInfPeriodList []
  set LethPeriodList []
  set InfPeriodList []
  set strainList []


  if experimental_rednoise = true
  [
    file-open "ered_noise.txt"
  set rednoise_list []
  while [ not file-at-end? ] [ set rednoise_list lput file-read rednoise_list ]
  file-close
  ]

  set BetaW_Middle BetaWithin
  set BetaB_Middle BetaWithin / 10 ; / 10 FIX ME
  set BetaWMid BetaWithin
  set BetaWHigh BetaWithinHV
  set BetaWLow BetaWithinLV
  set BetaBMid BetaWMid / 10
  set BetaBHigh BetaWHigh / 10
  set BetaBLow BetaWLow / 10

  set survivalTimeList ;;for non fixed survival times
  [
    16
    16
    8
    4
    2
    1
    0
    -1
    -2
    -4
    -8
    -16
    -32
  ]


  ;;; Fixed survival times selection
  if trade_off_type = "Lin"
  [
    set fixedExtTimesList
    [
      2 ;strain > 10
      4 ;strain  >8 <10
      6 ;strain +5 >6 <8
      8 ;strain +3 >4 <6
      10 ;strain +1 >2 <4
      12 ;strain -1 >0 <2
      14 ;strain 0
      16 ;strain <0 <-2
      18 ;strain -5 <-2 <-4
      20;strain -7 <-4 <-6
      22 ;strain -9 <-6 <-8
      24 ;strain -10 <-8 <-10
      26 ;strain -10 <= -10
    ]
  ]

  if trade_off_type = "Lin2"
  [
    set fixedExtTimesList
    [
      1 ;strain > 10
      2 ;strain  >8 <10
      3 ;strain +5 >6 <8
      4;strain +3 >4 <6
      6;strain +1 >2 <4
      8 ;strain -1 >0 <2
      10 ;strain 0
      11 ;strain <0 <-2
      12 ;strain -5 <-2 <-4
      13;strain -7 <-4 <-6
      14 ;strain -9 <-6 <-8
      15 ;strain -10 <-8 <-10
      16 ;strain -10 <= -10
    ]
  ]


    if trade_off_type = "Non lin A" ; used in current model version
  [
    set fixedExtTimesList
    [
      1 ;strain > 10
      2 ;strain  >8 <10
      2 ;strain +5 >6 <8
      3;strain +3 >4 <6
      3;strain +1 >2 <4
      6 ;strain -1 >0 <2
      8 ;strain 0
      9 ;strain <0 <-2
      10 ;strain -5 <-2 <-4
      11;strain -7 <-4 <-6
      11 ;strain -9 <-6 <-8
      13 ;strain -10 <-8 <-10
      14 ;strain -10 <= -10
    ]
  ]

  if trade_off_type = "Non lin B"
  [
    set fixedExtTimesList
    [
      2
      4
      6
      8
      12
      14
      16
      18
      19
      20
      21
      22
      23
    ]
  ]

  if trade_off_type = "Non lin C"
  [
    set fixedExtTimesList
    [
      1
      3
      6
      9
      12
      15
      18
      19
      20
      21
      22
      23
      24
    ]
  ]

    if trade_off_type = "Non lin D"
  [
    set fixedExtTimesList
    [
      1
      2
      4
      6
      9
      12
      15
      18
      19
      20
      21
      22
      22
    ]
  ]

   if trade_off_type = "Non lin E"
  [
    set fixedExtTimesList
    [
      1
      2
      3
      6
      12
      15
      18
      21
      22
      23
      24
      25
      26
    ]
  ]


     if trade_off_type = "Non lin F"
  [
    set fixedExtTimesList
    [
      2
      2
      3
      3
      4
      5
      10
      12
      14
      16
      20
      22
      26
    ]
  ]

     if trade_off_type = "Non lin G"
  [
    set fixedExtTimesList
    [
      2
      3
      4
      5
      8
      11
      13
      16
      19
      22
      23
      23
      23
    ]
  ]


  if transListSelect = true

  [
   set transList

  [
   0.04
   0.05
   0.06
   0.07
   0.08
   0.09
   0.12
   0.14
   0.18
   0.20
   0.20
   0.20
   0.20
  ]
  ]

   if transListSelect = false
  [
  set transList
  [
   0
  ]
  ]

set meanStrainOutput 0

end



;; LOAD LANDSCAPE _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Habitat patches varying in breeding capacity (Quality) are distributed randomly or homogeneous.
  ;; The darker the color is, the higher is the breeding capacity of the home range.


to Init_Landscape

;; LANDSCAPE ALGORITHM 1 (homogeneous habitat)
  if Landscape = "Generator-Homogeneous" [ ask patches [ set Quality MeanQuality ] ]


;; LANDSCAPE ALGORITHM 2 (random heterogeneous habitat)
  if Landscape = "Generator-Random"
  [
    ask patches with [pycor = min-pycor] [ set Quality MeanQuality ]
    ask patches with [pycor = max-pycor] [ set Quality MeanQuality ]

    ask patches with [pycor != min-pycor OR pycor != max-pycor]
    [
      ifelse Style = "discrete" [ set Quality random (MeanQuality * 2 + 1) ]
                                [ set Quality random-float MeanQuality * 2 ]
     set qmw MeanQuality
  ]

  ]


;; LANDSCAPE ALGORITHM 3 (import file)
  if Landscape = "File-Import-Path"
  [
    let x 0
    let y 0

    ;; file import
    ifelse seed-setup = "BS"
    [
      file-open (word "nlms/" Filename "_" (behaviorspace-run-number mod rep + 1) ".txt")
    ][
      file-open (word "nlms/" Filename ".txt")
    ]

    while [ not file-at-end? ]
    [
      while [ y < max-pycor + 1 ]
      [
        set x 0
        while [ x < max-pxcor + 1 ]
        [
          let next-value file-read ;; iterate over file entries seperated by white-space
          ask patch x y [ set Quality next-value ]
          set x x + 1
        ]
        set y y + 1
    ] ]
    file-close
  ]

  if Landscape = "File-Import-User"
  [
    let x 0
    let y 0

    ;; file import
    file-open user-file

    while [ not file-at-end? ]
    [
      while [ y < max-pycor + 1 ]
      [
        set x 0
        while [ x < max-pxcor + 1 ]
        [
          let next-value file-read ;; iterate over file entries seperated by white-space
          ask patch x y [ set Quality next-value ]
          set x x + 1
        ]
        set y y + 1
    ] ]
    file-close


  ]

  ask patches with [pycor = min-pycor] [ set Quality MeanQuality ]
  ask patches with [pycor = max-pycor] [ set Quality MeanQuality ]

  if Landscape = "File-Import-Path" OR Landscape = "File-Import-User"
  [
    let QualityFactor (MeanQuality) / mean [Quality] of patches with [pycor != min-pycor AND pycor != max-pycor]
    ask patches with [pycor != min-pycor AND pycor != max-pycor]
    [
      set Quality Quality * QualityFactor
      if Style = "discrete" [ set Quality floor Quality ]
  ] ]

  let x 0
  let y 0
  let patch_id 0

  while [ y < max-pycor + 1 ]
    [
      set x 0
      while [ x < max-pxcor + 1 ]
      [
        ask patch x y
        [
          if patch_id = 0
          [
            set id patch_id
            set patch_id patch_id + 1
          ]
          ifelse any? neighbors with [precision quality 2 = [precision quality 2] of myself AND id > 0]
          [
            let min_id min [id] of neighbors with [precision quality 2 = [precision quality 2] of myself AND id > 0]
            set id min_id
            ask neighbors with [precision quality 2 = [precision quality 2] of myself AND id > min_id] [ set id min_id ]
          ]
          [
            set id patch_id set patch_id patch_id + 1
        ] ]
        set x x + 1
      ]
      set y y + 1
  ]



  ask max-one-of patches [id][set max_id id]
  set Habitat patches with [Quality > 0]
  set Matrix patches with [Quality <= 0]
  ask Matrix [ set Quality 0 ]

  let Qual-MW mean [Quality] of patches

  ask Habitat
  [
    set Capacity (Quality * (20 / Qual-MW))   ; Quality * (mean number of 20 boars per cell / mean Quality) = 4.4445 boars per quality unit
    set NoBoars count turtles-here
    set NoRepro (count turtles-here with [is_female = 1 AND is_mother = 0 AND (AgeGroup = "Adult" OR AgeGroup = "Subdult")])
    set DispCells Habitat in-radius DispDist
    set Qual_around sum [Quality] of neighbors
    set is_habitat 1
    set qmw Qual-MW
    set tmpToInfectHere 0
 set InfPress 0
 set InfPressLowVir 0
 set InfPressHighVir 0
 set NoInfectedLowVir 0
 set NoInfectedHighVir 0
 set NoInfectedAdj 0
  ]
  ask Matrix [ set Capacity 0 ]

  ifelse landscapeGscale = true
  [
  ask patches [ set pcolor scale-color grey Quality 0 12 ]
  ]
  [
    ask patches [ set pcolor white ]
  ]

  ask patches [
    set infection_output_counter 0
    set infPressList []
    set virToGetList []
    set infectionList []
    set initialOutbreakCell false
  ]

end



;; CREATE START POPULATION ____________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; One wild boar and its herd is allocated to each habitat patch, group size is initialised as three times breeding capacity (Quality) of each patch.
  ;; Initial age distributions are taken from the results of a 100 years model run conducted by Kramer-Schadt et al. 2009.



to Init_Population

  crt (Herd% / 100 * count Habitat)
  [
     set Age 0
     set DemStat "dsStart"
     set EpiStat "esSusc"
     set is_shedding 0
     set is_breeder 0
     set is_mother 0
     set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
    set strain 999

  ]

  ask turtles [ ht ]

  while [any? turtles with [DemStat = "dsStart"]]
  [
    ask turtles with [HerdID = 0]
    [
      move-to one-of Habitat with [is_occupied = 0]
      set currHerdID currHerdID + 1
      set HerdID currHerdID
      set DemStat "dsResident"
      let offspring (ReleaseFactor * Quality - 1)
      hatch offspring
      if random-float 1 < (offspring - floor(offspring)) [ hatch 1 ]
      set is_occupied 1 ;patches-own
  ] ]

  while [any? turtles with [Age = 0]]
  [
    ask turtles with [Age = 0]
    [
      ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
      let RandomNo random-float 1
      ; generate variation in age of start population
      let Var round random-normal 0 (AgeBlur / 2)
      if Var <= AgeBlur AND Var >= (- AgeBlur)
      [
        ;; intial age distribution based on Kramer-Schadt et al. 2009
        if RandomNo <= 0.38                     [set Age  52 + Var]
        if RandomNo > 0.38 AND RandomNo <= 0.62 [set Age 104 + Var]
        if RandomNo > 0.62 AND RandomNo <= 0.77 [set Age 156 + Var]
        if RandomNo > 0.77 AND RandomNo <= 0.86 [set Age 208 + Var]
        if RandomNo > 0.86 AND RandomNo <= 0.92 [set Age 260 + Var]
        if RandomNo > 0.92 AND RandomNo <= 0.95 [set Age 312 + Var]
        if RandomNo > 0.95 AND RandomNo <= 0.97 [set Age 364 + Var]
        if RandomNo > 0.97 AND RandomNo <= 0.98 [set Age 416 + Var]
        if RandomNo > 0.98 AND RandomNo <= 0.99 [set Age 468 + Var]
        if RandomNo > 0.99 AND RandomNo <= 1    [set Age 520 + Var]
      ]

    ;; age classes males
    if Age <= 21 AND is_female = 0 [ set AgeGroup "Piglet" ]
    if Age > 21 AND Age <= 104 AND is_female = 0 [ set AgeGroup "Subadult" ]
    if Age > 104 AND is_female = 0
    [
      set AgeGroup "Adult"
      set DemStat "dsDispM_adult"
    ]

    ;; age classes females
    if Age <= 34 AND is_female = 1 [ set AgeGroup "Piglet" ]
    if Age > 34 AND Age <= 52 AND is_female = 1 [ set AgeGroup "Subadult" ]
    if Age > 52 AND is_female = 1 [ set AgeGroup "Adult" ]

  ] ]

;; OUTPUT
  set PopStructList []
  ask turtles
  [
    set PopStructList fput Age PopStructList
    set NoGroup count turtles-here
  ]
  set PopStructList map [ ?1 -> ?1 / 52 ] PopStructList

end






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  GO PROCEDURES                                                                                                                                                                                                                                            ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  DISEASE PROCEDURES                                                                                                                                                                                                                                       ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; DISEASE RELEASE ____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; One randomly selected wild boar or selected group (cell) in the middle of the top row gets lethally infected in a random week ("WeekRelease") of a given year (based on input "YearRelease").

to Disease_Release

  infectionCheck

  ask patch 12 25
  [
    set initialOutbreakCell true
    ask n-of 3 neighbors [set initialOutbreakCell true]
  ]

  ask patches with [initialOutbreakCell = true]
  [
    set InfX pxcor
    set InfY pycor
    set is_infectious 1
    set is_infected 1
    set outbreak_group count turtles-here
    set outbreak_roaming count turtles-here with [ DemStat = "dsDispM_adult" ]
    set outbreak_dens5 count turtles in-radius 1.5
    set outbreak_dens12 count turtles in-radius 12
    ;;equal outbreak of all strains

    let susHere count turtles-here
    let thrd susHere ;/ 3

    if Roaming != "OFF"
    [
     let tmpMoThrd outbreak_roaming / 3
      ask at-most-n-of thrd turtles-here with [hasActiveStrain = false AND DemStat = "dsDispM_adult"]
      [
        set initialOutbreak true
        set strainVirToGet ((random 5) + 1) + 7
        Disease_Infection
      ]
      ask at-most-n-of thrd turtles-here with [hasActiveStrain = false AND DemStat = "dsDispM_adult"]
      [
        set initialOutbreak true
        set strainVirToGet ((random 2) + 1) + 5
        Disease_Infection
      ]
      ask at-most-n-of thrd turtles-here with [hasActiveStrain = false AND DemStat = "dsDispM_adult"]
      [
        set initialOutbreak true
        set strainVirToGet ((random 5) + 1)
        Disease_Infection
      ]
    ]

    ask at-most-n-of thrd turtles-here with [hasActiveStrain = false]
      [
        set initialOutbreak true
        set strainVirToGet 5
        Disease_Infection
    ]

  testcounter
  ]
end



;; DISEASE TRANSMISSION _______________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; First, infection pressure is calculated for all habitat cells based on infected indivudals in the herd and the eight herds surrounding the home range cell.
  ;; Afterwards, the model iterates over all individuals and stochastically sets susceptible indivudals to infected if a drawn random number (0-1) is smaller than the infection pressure of the cell.
  ;; The transmission parameter was reversly fitted to the estimated disease spread velocity of approx. 8 km per quarter (Rossi et al. 2010). The transmission probabilities were assigned constant as 0.0208 within a herd and 0.00208 between herds.
  ;; Based on the case fatality rate (CaseFatality), on infection the host is stochastically assigned either as lethally infected (esLeth) or transiently infected (esTrans).
  ;; CaseFatality for lethally infected individuals applies unchanged for yearlings (subadults), is decreased for adults (CaseFatality ^ 2) and increased for piglets (sqrt CaseFatality) to represent age-dependent disease outcomes (Dahle & Liess 1992).
  ;; Lethally infected hosts receive their individual infectious period / survival times (in weeks) which are drawn from an exponential distribution with mean mue (resulting in TickLeth).


to Disease_Transmission

  set NewLeth 0
  set NewTrans 0

  ;; update infection status of cells
  ask Habitat [ set is_infectious 0 ]
  ask Habitat with [any? turtles-here with [is_shedding = 1]]
  [
    set is_infectious 1
    set is_infected 1
  ]

  ask Habitat
  [
    ifelse Roaming = "OFF"
    [

      set AllSus count turtles-here with [EpiStat = "esSusc"]

      set inMin count turtles-here with [is_shedding = 1 AND strain <= -10 AND hasActiveStrain = true]
      set in5 count turtles-here with [is_shedding = 1 AND strain <= -8 AND strain > -10 AND hasActiveStrain = true]
      set in6 count turtles-here with [is_shedding = 1 AND strain <= -6 AND strain > -8 AND hasActiveStrain = true]
      set in7 count turtles-here with [is_shedding = 1 AND strain <= -4 AND strain > -6 AND hasActiveStrain = true]
      set in8 count turtles-here with [is_shedding = 1 AND strain <= -2 AND strain > -4 AND hasActiveStrain = true]
      set in9 count turtles-here with [is_shedding = 1 AND strain < 0 AND strain > -2 AND hasActiveStrain = true]
      set in10 count turtles-here with [is_shedding = 1 AND strain = 0 AND hasActiveStrain = true]
      set in11 count turtles-here with [is_shedding = 1 AND strain > 0 AND strain <= 2 AND hasActiveStrain = true]
      set in12 count turtles-here with [is_shedding = 1 AND strain > 2 AND strain <= 4 AND hasActiveStrain = true]
      set in13 count turtles-here with [is_shedding = 1 AND strain > 4 AND strain <= 6 AND hasActiveStrain = true]
      set in14 count turtles-here with [is_shedding = 1 AND strain > 6 AND strain <= 8 AND hasActiveStrain = true]
      set in15 count turtles-here with [is_shedding = 1 AND strain > 8 AND strain <= 10 AND hasActiveStrain = true]
      set inMax count turtles-here with [is_shedding = 1 AND strain > 10 AND hasActiveStrain = true]

    ]
    [;; roaming ON

      set AllSus count turtles-here with [EpiStat = "esSusc"]
      set inMin count turtles-here with [is_shedding = 1 AND strain <= -10 AND hasActiveStrain = true]
      set in5 count turtles-here with [is_shedding = 1 AND strain <= -8 AND strain > -10 AND hasActiveStrain = true]
      set in6 count turtles-here with [is_shedding = 1 AND strain <= -6 AND strain > -8 AND hasActiveStrain = true]
      set in7 count turtles-here with [is_shedding = 1 AND strain <= -4 AND strain > -6 AND hasActiveStrain = true]
      set in8 count turtles-here with [is_shedding = 1 AND strain <= -2 AND strain > -4 AND hasActiveStrain = true]
      set in9 count turtles-here with [is_shedding = 1 AND strain < 0 AND strain > -2 AND hasActiveStrain = true]
      set in10 count turtles-here with [is_shedding = 1 AND strain = 0 AND hasActiveStrain = true]
      set in11 count turtles-here with [is_shedding = 1 AND strain > 0 AND strain <= 2 AND hasActiveStrain = true]
      set in12 count turtles-here with [is_shedding = 1 AND strain > 2 AND strain <= 4 AND hasActiveStrain = true]
      set in13 count turtles-here with [is_shedding = 1 AND strain > 4 AND strain <= 6 AND hasActiveStrain = true]
      set in14 count turtles-here with [is_shedding = 1 AND strain > 6 AND strain <= 8 AND hasActiveStrain = true]
      set in15 count turtles-here with [is_shedding = 1 AND strain > 8 AND strain <= 10 AND hasActiveStrain = true]
      set inMax count turtles-here with [is_shedding = 1 AND strain > 10 AND hasActiveStrain = true]

    ]
  ]

  ask Habitat
  [
    ifelse Roaming = "OFF"
    [
      set infectionList []
      set infPressList []
      set virToGetList []

      set ip10 (1 - (1 -  BetaW_Middle) ^ in10 * (1 - ( BetaB_Middle)) ^ sum [in10] of neighbors)
      set ip9 (1 - (1 -  (BetaW_Middle - ( 0))) ^ in9 * (1 - (BetaB_Middle - (0 ) )) ^ sum [in9] of neighbors) ;;; 5% difference between each "brackets" beta
      set ip8 (1 - (1 -  (BetaW_Middle - (0 ))) ^ in8 * (1 - (BetaB_Middle - (0 ) )) ^ sum [in8] of neighbors)
      set ip7 (1 - (1 -  (BetaW_Middle - ( 0))) ^ in7 * (1 - (BetaB_Middle - (0 ) )) ^ sum [in7] of neighbors)
      set ip6 (1 - (1 -  (BetaW_Middle - ( 0))) ^ in6 * (1 - (BetaB_Middle - (0 ) )) ^ sum [in6] of neighbors)
      set ip5 (1 - (1 -  (BetaW_Middle - ( 0))) ^ in5 * (1 - (BetaB_Middle - ( 0) )) ^ sum [in5] of neighbors)
      set ipMin (1 - (1 -  (BetaW_Middle - ( 0))) ^ inMin * (1 - (BetaB_Middle - ( 0) )) ^ sum [inMin] of neighbors)

      set ip11 (1 - (1 -  (BetaW_Middle + ( 0))) ^ in11 * (1 - (BetaB_Middle + (0 ) )) ^ sum [in11] of neighbors)
      set ip12 (1 - (1 -  (BetaW_Middle + ( 0))) ^ in12 * (1 - (BetaB_Middle + (0 ) )) ^ sum [in12] of neighbors)
      set ip13 (1 - (1 -  (BetaW_Middle + ( 0))) ^ in13 * (1 - (BetaB_Middle + ( 0) )) ^ sum [in13] of neighbors)
      set ip14 (1 - (1 -  (BetaW_Middle + ( 0))) ^ in14 * (1 - (BetaB_Middle + ( 0) )) ^ sum [in14] of neighbors)
      set ip15 (1 - (1 -  (BetaW_Middle + ( 0))) ^ in15 * (1 - (BetaB_Middle + ( 0) )) ^ sum [in15] of neighbors)
      set ipMax (1 - (1 - (BetaW_Middle + ( 0))) ^ inMax * (1 - (BetaB_Middle + ( 0) )) ^ sum [inMax] of neighbors)

 if ipMin != 0
    [
        set ipMin ipMin + item 0 transList
    ]

      if ip5 != 0
    [
         set ip5 ip5 + item 1 transList
    ]

      if ip6 != 0
    [
        set ip6 ip6 + item 2 transList
    ]

      if ip7 != 0
    [
        set ip7 ip7 + item 3 transList
    ]

      if ip8 != 0
    [
        set ip8 ip8 + item 4 transList
    ]

    if ip9 != 0
    [
        set ip9 ip9 + item 5 transList
    ]
      if ip10 != 0
    [
        set ip10 ip10 + item 6 transList
    ]

    if ip11 != 0
    [
    set ip11 ip11 + item 7 transList
      ]
      if ip11 <= 0 [set ip11 0.0001]

    if ip12 != 0
    [
    set ip12 ip12 + item 8 transList
      ]

    if ip13 != 0
    [
      set ip13 ip13 + item 9 transList
      ]

    if ip14 != 0
    [
    set ip14 ip14 + item 10 transList
      ]

    if ip15 != 0
    [
    set ip15 ip15 + item 11 transList
      ]

    if ipMax != 0
    [
     set ipMax ipMax + item 12 transList
      ]

      set infPressList fput ipMin  infPressList
      set virToGetList fput 1 virToGetList
      set infPressList fput ip5 infPressList
      set virToGetList fput 2 virToGetList
      set infPressList fput ip6 infPressList
      set virToGetList fput 3 virToGetList
      set infPressList fput ip7 infPressList
      set virToGetList fput 4 virToGetList
      set infPressList fput ip8 infPressList
      set virToGetList fput 5 virToGetList
      set infPressList fput ip9 infPressList
      set virToGetList fput 6 virToGetList
      set infPressList fput ip10  infPressList
      set virToGetList fput 7 virToGetList
      set infPressList fput ip11 infPressList
      set virToGetList fput 8 virToGetList
      set infPressList fput ip12 infPressList
      set virToGetList fput 9 virToGetList
      set infPressList fput ip13 infPressList
      set virToGetList fput 10 virToGetList
      set infPressList fput ip14 infPressList
      set virToGetList fput 11 virToGetList
      set infPressList fput ip15 infPressList
      set virToGetList fput 12 virToGetList
      set infPressList fput ipMax infPressList
      set virToGetList fput 13 virToGetList

      set infectionList lput (list infPressList virToGetList) infectionList

    ]   ;c: all individuals are equal -> group's cell and neighbors
    [ ;Roaming ON
      set infectionList []
      set infPressList []
      set virToGetList []

      set ip10 (1 - (1 -  BetaW_Middle) ^ in10 )
      set ip9 (1 - (1 -  (BetaW_Middle - (BetaW_Middle ))) ^ in9)  ;;; 5% difference between each "brackets" beta
      set ip8 (1 - (1 -  (BetaW_Middle - (BetaW_Middle ))) ^ in8)
      set ip7 (1 - (1 -  (BetaW_Middle - (BetaW_Middle ))) ^ in7)
      set ip6 (1 - (1 -  (BetaW_Middle - (BetaW_Middle ))) ^ in6)
      set ip5 (1 - (1 -  (BetaW_Middle - (BetaW_Middle ))) ^ in5)
      set ipMin (1 - (1 -  (BetaW_Middle - (BetaW_Middle ))) ^ inMin)

      set ip11 (1 - (1 -  (BetaW_Middle + (BetaW_Middle ))) ^ in11)
      set ip12 (1 - (1 -  (BetaW_Middle + (BetaW_Middle ))) ^ in12)
      set ip13 (1 - (1 -  (BetaW_Middle + (BetaW_Middle ))) ^ in13)
      set ip14 (1 - (1 -  (BetaW_Middle + (BetaW_Middle ))) ^ in14)
      set ip15 (1 - (1 -  (BetaW_Middle + (BetaW_Middle ))) ^ in15)
      set ipMax (1 - (1 - (BetaW_Middle + (BetaW_Middle ))) ^ inMax)

       if ipMin != 0
    [
      set ipMin ipMin - (transmission_factor + 0.066)
      if ipMin <= 0 [set ipMin 0.0001]
    ]

      set ip10 ip10

    if ip9 != 0
    [
      set ip9 ip9 - (transmission_factor + 0.011)
      if ip9 <= 0 [set ip9 0.0001]
    ]



    if ip8 != 0
    [
      set ip8 ip8 - (transmission_factor + 0.022)
      if ip8 <= 0 [set ip8 0.0001]
    ]

     if ip7 != 0
    [
      set ip7 ip7 - (transmission_factor + 0.033)
      if ip7 <= 0 [set ip7 0.0001]
    ]

      if ip6 != 0
    [
      set ip6 ip6 - (transmission_factor + 0.044)
      if ip6 <= 0 [set ip6 0.0001]
    ]

        if ip5 != 0
    [
      set ip5 ip5 - (transmission_factor + 0.55)
      if ip5 <= 0 [set ip5 0.0001]
    ]

    if ip11 != 0
    [
    set ip11 ip11 + (transmission_factor + 0.011)
    ]

    if ip12 != 0
    [
    set ip12 ip12 + (transmission_factor + 0.022)
    ]

    if ip13 != 0
    [
    set ip13 ip13 + (transmission_factor + 0.033)
    ]

    if ip14 != 0
    [
     set ip14 ip14 + (transmission_factor + 0.044)
    ]

    if ip15 != 0
    [
    set ip15 ip15 + (transmission_factor + 0.055)
    ]

    if ipMax != 0
    [
    set ipMax ipMax + (transmission_factor + 0.066)
    ]


      set infPressList fput ipMin  infPressList
      set virToGetList fput 1 virToGetList
      set infPressList fput ip5 infPressList
      set virToGetList fput 2 virToGetList
      set infPressList fput ip6 infPressList
      set virToGetList fput 3 virToGetList
      set infPressList fput ip7 infPressList
      set virToGetList fput 4 virToGetList
      set infPressList fput ip8 infPressList
      set virToGetList fput 5 virToGetList
      set infPressList fput ip9 infPressList
      set virToGetList fput 6 virToGetList
      set infPressList fput ip10  infPressList
      set virToGetList fput 7 virToGetList
      set infPressList fput ip11 infPressList
      set virToGetList fput 8 virToGetList
      set infPressList fput ip12 infPressList
      set virToGetList fput 9 virToGetList
      set infPressList fput ip13 infPressList
      set virToGetList fput 10 virToGetList
      set infPressList fput ip14 infPressList
      set virToGetList fput 11 virToGetList
      set infPressList fput ip15 infPressList
      set virToGetList fput 12 virToGetList
      set infPressList fput ipMax infPressList
      set virToGetList fput 13 virToGetList

      set infectionList lput (list infPressList virToGetList) infectionList

    ]   ;c: age- and sex-dependent differences -> only group's cell
  ]


  infectionCheck

  ask Habitat [ if is_infectious = 1 [ set inf_count inf_count + 1 ] ]

  ;; OUTPUT
  set CurrInfPeriodList []
  ask turtles with [EpiStat = "esLeth" AND hasActiveStrain = TRUE] [ set CurrInfPeriodList fput (TickLeth - ticks) CurrInfPeriodList ]
  set strainList []
  ask turtles with [EpiStat = "esLeth" AND hasActiveStrain = TRUE] [ set strainList fput strain strainList ]

end



;; DISEASE STATE TRANSITION ___________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Transient shedders (esTrans) are converted to immune (esImm) after T_trans.
  ;; Individuals protected by maternal antibodies (esImmMat) turn suscpetible (esSusc) after reaching the maximum duration of immunity of 12 weeks.


to Disease_Transition

  ask turtles with [EpiStat = "esLeth" AND TickInf + T_trans <= ticks AND initialOutbreak = false AND transientCheck = false AND hasactiveStrain = true] ;AND strain < 0]
  [
  if AgeGroup = "Piglet"   [ set IndCaseFatality sqrt CaseFatality ]
  if  AgeGroup = "Subadult" [ set IndCaseFatality CaseFatality ]
  if  AgeGroup = "Adult"    [ set IndCaseFatality CaseFatality ^ 2 ]



    if random-float 1 < IndCaseFatality
    [
      set EpiStat "esImm"
      set is_shedding 0
      set hasActiveStrain false
      set strain 999
      set markedForInfection false
      set strainVirToGet 0
    ]
    set transientCheck true
  ]

  ;; maternally protected pigs
  ask turtles with [EpiStat = "esImmMat" AND Age >= T_anti]
  [
    set EpiStat "esSusc"
    set is_matAB 0
  ]

  ask turtles with [is_matAB = 1 AND Age > T_anti - 4]   ; 1 month with varying immunity for piglets
  [
    if random-float 1 < 0.5
    [
      set EpiStat "esSusc"
      set is_matAB 0
  ] ]

  ;; OUTPUT
  set InfStructList []
  ask turtles with [is_shedding = 1]  [ set InfStructList fput Age InfStructList ]
  set InfStructList map [ ?1 -> ?1 / 52 ] InfStructList

  set InfTransStructList []
  ask turtles with [EpiStat = "esTrans"]  [ set InfTransStructList fput Age InfTransStructList ]
  set InfTransStructList map [ ?1 -> ?1 / 52 ] InfTransStructList

  set InfLethStructList []
  ask turtles with [EpiStat = "esLeth"]  [ set InfLethStructList fput Age InfLethStructList ]
  set InfLethStructList map [ ?1 -> ?1 / 52 ] InfLethStructList

end



;; INFECTIOUS PROCESS _________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; turtle procedure

  ;; Individual infection based on case fatality



to Disease_Infection

  ;; initial outbreak

  ifelse initialOutbreak = true
  [
    if strainVirToGet = 1
    [
     set strain -10
    ]
      if strainVirToGet = 2
    [
     set strain (random 3) - 10
    ]
    if strainVirToGet = 3
    [
     set strain (random 3) - 8
    ]
    if strainVirToGet = 4
    [
     set strain (random 3) - 6
    ]
      if strainVirToGet = 5
    [
     set strain (random 3) - 4
    ]
    if strainVirToGet = 6
    [
     set strain (random 3) - 2
    ]
    if strainVirToGet = 7
    [
     set strain 0
    ]
      if strainVirToGet = 8
    [
     set strain (random 3)
    ]
    if strainVirToGet = 9
    [
     set strain (random 3) + 2
    ]
    if strainVirToGet = 10
    [
     set strain (random 3) + 4
    ]
      if strainVirToGet = 11
    [
     set strain (random 3) + 6
    ]
    if strainVirToGet = 12
    [
     set strain (random 3) + 8
    ]
    if strainVirToGet = 13
    [
     set strain (10)
    ]

  set HasActiveStrain true

  if EpiStat = "esSusc" AND AgeGroup = "Piglet"   [ set IndCaseFatality sqrt CaseFatality ]
  if EpiStat = "esSusc" AND AgeGroup = "Subadult" [ set IndCaseFatality CaseFatality ]
  if EpiStat = "esSusc" AND AgeGroup = "Adult"    [ set IndCaseFatality CaseFatality ^ 2 ]

  ifelse random-float 1 <= 2 AND hasActiveStrain = true; IndCaseFatality ;force outbreak during disease release
  [

    set EpiStat "esLeth"
    set is_shedding 1
   ; set TickLeth ticks + 1 + floor (- (mue - 0.5) * ln random-float 1)
    set NewLeth NewLeth + 1
    let tempTickLeth 0

    set is_infectious 1  ;patches-own
    set is_infected 1  ;patches-own
    if distancexy InfX InfY > InfDist [ set InfDist distancexy InfX InfY ]  ;patches-own

      if fixedSurvivalTimes = false
    [
    set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)

    let Vir item (strainVirToGet - 1) survivalTimeList

    if strain < -10 [set vir item 0 survivalTimeList]

    set TickLeth round(tempTickLeth + Vir)
    ]

    if fixedSurvivalTimes = true
    [

    set TickLeth get_survival_time strain

    ]

    set LethPeriodList fput (TickLeth - ticks) LethPeriodList  ; list for lethally infected individuals only
    set InfPeriodList fput (TickLeth - ticks) InfPeriodList    ; list for lethally and transient infected individuals
    set  initialOutbreak  false
  ][ ; else
  ]
  ];;; END initial outbreak
  [

  getStrain

  if evolution = true
  [
    evolve
  ]

  if EpiStat = "esSusc" AND AgeGroup = "Piglet"   [ set IndCaseFatality sqrt CaseFatality ]
  if EpiStat = "esSusc" AND AgeGroup = "Subadult" [ set IndCaseFatality CaseFatality ]
  if EpiStat = "esSusc" AND AgeGroup = "Adult"    [ set IndCaseFatality CaseFatality ^ 2 ]

  ifelse random-float 1 < IndCaseFatality ; removed transient infection
  [
    if strain = 999 [user-message (word "ERROR 001 - Disease infection function - strain transmission failure")]
    if strainVirToGet = 0 Or strainVirToGet > 13 [user-message (word "ERROR 002 - Disease infection function - strain virulence allocation failure")]
    set EpiStat "esLeth"
    set is_shedding 1
    let tempTickLeth 0


    ;;; random like survival time allocation
    if fixedSurvivalTimes = false
    [
    set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)

    let Vir item (strainVirToGet - 1) survivalTimeList

    if strain < -10 [set vir item 0 survivalTimeList]

    set TickLeth round(tempTickLeth + Vir)
    ]

    ;;; END random survival times

    ;;;;using the strain/virulence to get the fixed survival time from a list
    if fixedSurvivalTimes = true
    [
    set TickLeth get_survival_time strain
    ]
    ;;;;END fixed survival times


    if TickLeth < 0
      [set TickLeth 2] ; changed from 1 to 2 FIX ME
    set TickLeth TickLeth + ticks
    set NewLeth NewLeth + 1
    set LethPeriodList fput (TickLeth - ticks) LethPeriodList  ; list for lethally infected individuals only
    set InfPeriodList fput (TickLeth - ticks) InfPeriodList    ; list for lethally and transient infected individuals
    set is_infectious 1  ;patches-own
    set is_infected 1  ;patches-own
    if distancexy InfX InfY > InfDist [ set InfDist distancexy InfX InfY ]  ;patches-own
  ][

  ]
  ]

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  HOST PROCEDURES                                                                                                                                                                                                                                          ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; REPRODUCTION _______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Adult and subadult females reproduce only once a year depending on their age class.
  ;; In the first week of the year females are checked wether they are able breed (age + no piglets) and if there is at least one reproducible male within their range (> 2 years).
  ;; Litter size is drawn from a truncated normal distribution and reduced to a constant fraction for infected individuals.
  ;; Depending on the disease state of breeding indiviuals its piglets disease states are adjusted.


to Host_Reproduction

  ask Habitat [ set NoRepro count turtles-here with [is_female = 1 AND is_mother = 0 AND Age > 34] ]
    ;; ... for habitat patches with at least one male to reproduce ...
    ask Habitat
    [
      ;; ... and for all females which can breed
      let num-fem NoRepro
      ;; if number of potential breeders is higher than breeding capacity
      if (num-fem > counter_d)
      [
        ;; … sort breeders descending by age …
        let sort-fem sort-on [(- Age)] turtles-here with [is_female = 1 AND is_mother = 0 AND Age > 34]
        set num-fem floor counter_d
        if length sort-fem > 0
        [
          while [num-fem > 0]
          [
            ;; … and let num-fem of them (difference between capacity and number of breeders) reproduce
            ask first sort-fem
            [
              set is_breeder 1
              set ChanceBreeding random-float 1
            ]
            set sort-fem remove-item 0 sort-fem
            set num-fem (num-fem - 1)
      ] ] ]

      ;; if number of breeders is equal or below than breeding capacity all breeders reproduce
      if (num-fem > 0 AND num-fem <= counter_d)
      [
        ask turtles-here with [is_female = 1 AND AgeGroup = "Adult" AND is_mother = 0]
        [
          set is_breeder 1
          set ChanceBreeding random-float 1
  ] ] ]

  ;; cumulative reproduction probabilities for each month (January to September)
  ifelse week <= 4 [ set ProbRepro item 0 ReproList  Host_prob-Reproduce ][
    ifelse week <= 8 [ set ProbRepro item 1 ReproList  Host_prob-Reproduce ][
      ifelse week <= 13 [ set ProbRepro item 2 ReproList  Host_prob-Reproduce ][
        ifelse week <= 17 [ set ProbRepro item 3 ReproList  Host_prob-Reproduce ][
          ifelse week <= 21 [ set ProbRepro item 4 ReproList  Host_prob-Reproduce ][
            ifelse week <= 26 [ set ProbRepro item 5 ReproList  Host_prob-Reproduce ][
              ifelse week <= 30 [ set ProbRepro item 6 ReproList  Host_prob-Reproduce ][
                ifelse week <= 34 [ set ProbRepro item 7 ReproList  Host_prob-Reproduce ][
                  ifelse week <= 39 [ set ProbRepro item 8 ReproList  Host_prob-Reproduce ][
                    ifelse week <= 43 [ set ProbRepro item 9 ReproList  Host_prob-Reproduce ][
                      ifelse week <= 47 [ set ProbRepro item 10 ReproList  Host_prob-Reproduce ][
                        if week <= 52 [ set ProbRepro item 11 ReproList  Host_prob-Reproduce
  ] ] ] ] ] ] ] ] ] ] ] ]

end


to Host_prob-Reproduce
  ;; observer procedure
    let strTra 99999
    let strVirG 99999

    ask turtles with [is_breeder = 1]
    [
    if hasActiveStrain = true
    [
      set strTra strain
      set strVirG strainVirToGet
    ]
      ;; stochastic effect in probabilities of reproduction
      if ChanceBreeding < ProbRepro
      [
        let RandomNo random-float 1

        ;; stochastic effect in number of piglets
        if (RandomNo > 0.01305859 AND RandomNo <= 0.082208805)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse is_shedding = 1 AND hasActiveStrain = true
          [
            hatch floor (1 * FertRedInf)   ;; round to integer
            [
              set DemStat "dsOffspring"
              set family who
              set float_counter 0
              set hasActiveStrain false
              set initialOutbreak false
              set transientCheck false
              set markedForInfection false
              set strainVirToGet 0
              set counterGotInfected false
              set strain 999
            if strTra != 99999
            [
              ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
              set Age 0
              set AgeGroup "Piglet"
              set is_breeder 0
              set strain strTra
              set strainVirToGet strVirG
              set EpiStat "esLeth"
              set is_shedding 1
              set DemStat "dsResident"
              set hasActiveStrain true
              let tempTickLeth 0
              if fixedSurvivalTimes = false
              [
                set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)
                let Vir item (strainVirToGet - 1) survivalTimeList
                if strain < -10 [set vir item 0 survivalTimeList]
                set TickLeth round(tempTickLeth + Vir)
              ]
              if fixedSurvivalTimes = true
              [
                set TickLeth get_survival_time strain
              ]
              set NewLeth NewLeth + 1
              set LethPeriodList fput (TickLeth - ticks) LethPeriodList
              set InfPeriodList fput (TickLeth - ticks) InfPeriodList
            ]

          ]
        ]
          [
          ;; if the breeder is not infected
            hatch 1
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
        ] ] ]

        if (RandomNo > 0.082208805 AND RandomNo <= 0.24510931)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse is_shedding = 1 AND hasActiveStrain = true
          [
            hatch floor (2 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
           set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
             if strTra != 99999
            [
              ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
              set Age 0
              set AgeGroup "Piglet"
              set is_breeder 0
              set strain strTra
              set strainVirToGet strVirG
              set EpiStat "esLeth"
              set is_shedding 1
              set DemStat "dsResident"
              set hasActiveStrain true
              let tempTickLeth 0
              if fixedSurvivalTimes = false
              [
                set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)
                let Vir item (strainVirToGet - 1) survivalTimeList
                if strain < -10 [set vir item 0 survivalTimeList]
                set TickLeth round(tempTickLeth + Vir)
              ]
              if fixedSurvivalTimes = true
              [
                set TickLeth get_survival_time strain
              ]
              set NewLeth NewLeth + 1
              set LethPeriodList fput (TickLeth - ticks) LethPeriodList
              set InfPeriodList fput (TickLeth - ticks) InfPeriodList
            ]
            ]
          ][
            ;; if the breeder is not infected
            hatch 2
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
        ] ] ]

        if (RandomNo > 0.24510931 AND RandomNo <= 0.495050085)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse is_shedding = 1 AND hasActiveStrain = true
          [
            hatch floor (3 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
             if strTra != 99999
            [
              ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
              set Age 0
              set AgeGroup "Piglet"
              set is_breeder 0
              set strain strTra
              set strainVirToGet strVirG
              set EpiStat "esLeth"
              set is_shedding 1
              set DemStat "dsResident"
              set hasActiveStrain true
              let tempTickLeth 0
              if fixedSurvivalTimes = false
              [
                set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)
                let Vir item (strainVirToGet - 1) survivalTimeList
                if strain < -10 [set vir item 0 survivalTimeList]
                set TickLeth round(tempTickLeth + Vir)
              ]
              if fixedSurvivalTimes = true
              [
                set TickLeth get_survival_time strain
              ]
              set NewLeth NewLeth + 1
              set LethPeriodList fput (TickLeth - ticks) LethPeriodList
              set InfPeriodList fput (TickLeth - ticks) InfPeriodList
            ]
            ]
          ][
            ;; if the breeder is not infected
            hatch 3
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
        ] ] ]

        if (RandomNo > 0.495050085 AND RandomNo <= 0.744990859)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse is_shedding = 1 AND hasActiveStrain = true
          [
            hatch floor (4 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
             if strTra != 99999
            [
              ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
              set Age 0
              set AgeGroup "Piglet"
              set is_breeder 0
              set strain strTra
              set strainVirToGet strVirG
              set EpiStat "esLeth"
              set is_shedding 1
              set DemStat "dsResident"
              set hasActiveStrain true
              let tempTickLeth 0
              if fixedSurvivalTimes = false
              [
                set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)
                let Vir item (strainVirToGet - 1) survivalTimeList
                if strain < -10 [set vir item 0 survivalTimeList]
                set TickLeth round(tempTickLeth + Vir)
              ]
              if fixedSurvivalTimes = true
              [
                set TickLeth get_survival_time strain
              ]
              set NewLeth NewLeth + 1
              set LethPeriodList fput (TickLeth - ticks) LethPeriodList
              set InfPeriodList fput (TickLeth - ticks) InfPeriodList
            ]
          ] ]
          [
            ;; if the breeder is not infected
            hatch 4
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
        ] ] ]

        if (RandomNo > 0.744990859 AND RandomNo <= 0.907891364)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse is_shedding = 1 AND hasActiveStrain = true
          [
            hatch floor (5 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
             if strTra != 99999
            [
              ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
              set Age 0
              set AgeGroup "Piglet"
              set is_breeder 0
              set strain strTra
              set strainVirToGet strVirG
              set EpiStat "esLeth"
              set is_shedding 1
              set DemStat "dsResident"
              set hasActiveStrain true
              let tempTickLeth 0
              if fixedSurvivalTimes = false
              [
                set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)
                let Vir item (strainVirToGet - 1) survivalTimeList
                if strain < -10 [set vir item 0 survivalTimeList]
                set TickLeth round(tempTickLeth + Vir)
              ]
              if fixedSurvivalTimes = true
              [
                set TickLeth get_survival_time strain
              ]
              set NewLeth NewLeth + 1
              set LethPeriodList fput (TickLeth - ticks) LethPeriodList
              set InfPeriodList fput (TickLeth - ticks) InfPeriodList
            ]
            ]
          ][
            ;; if the breeder is not infected
            hatch 5
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
        ] ] ]

        if (RandomNo > 0.907891364 AND RandomNo <= 0.977041579)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse is_shedding = 1 AND hasActiveStrain = true
          [
            hatch floor (6 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
             if strTra != 99999
            [
              ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
              set Age 0
              set AgeGroup "Piglet"
              set is_breeder 0
              set strain strTra
              set strainVirToGet strVirG
              set EpiStat "esLeth"
              set is_shedding 1
              set DemStat "dsResident"
              set hasActiveStrain true
              let tempTickLeth 0
              if fixedSurvivalTimes = false
              [
                set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)
                let Vir item (strainVirToGet - 1) survivalTimeList
                if strain < -10 [set vir item 0 survivalTimeList]
                set TickLeth round(tempTickLeth + Vir)
              ]
              if fixedSurvivalTimes = true
              [
                set TickLeth get_survival_time strain
              ]
              set NewLeth NewLeth + 1
              set LethPeriodList fput (TickLeth - ticks) LethPeriodList
              set InfPeriodList fput (TickLeth - ticks) InfPeriodList
            ]
            ]
          ][
            ;; if the breeder is not infected
            hatch 6
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
        ] ] ]

        if (RandomNo > 0.977041579 AND RandomNo <= 0.996144138)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse is_shedding = 1 AND hasActiveStrain = true
          [
            hatch floor (7 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
             if strTra != 99999
            [
              ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
              set Age 0
              set AgeGroup "Piglet"
              set is_breeder 0
              set strain strTra
              set strainVirToGet strVirG
              set EpiStat "esLeth"
              set is_shedding 1
              set DemStat "dsResident"
              set hasActiveStrain true
              let tempTickLeth 0
              if fixedSurvivalTimes = false
              [
                set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)
                let Vir item (strainVirToGet - 1) survivalTimeList
                if strain < -10 [set vir item 0 survivalTimeList]
                set TickLeth round(tempTickLeth + Vir)
              ]
              if fixedSurvivalTimes = true
              [
                set TickLeth get_survival_time strain
              ]
              set NewLeth NewLeth + 1
              set LethPeriodList fput (TickLeth - ticks) LethPeriodList
              set InfPeriodList fput (TickLeth - ticks) InfPeriodList
            ]
          ] ]
          [
            ;; if the breeder is not infected
            hatch 7
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
        ] ] ]

        if (RandomNo > 0.996144138 AND RandomNo <= 0.999575749)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse is_shedding = 1 AND hasActiveStrain = true
          [
            hatch floor (8 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
             if strTra != 99999
            [
              ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
              set Age 0
              set AgeGroup "Piglet"
              set is_breeder 0
              set strain strTra
              set strainVirToGet strVirG
              set EpiStat "esLeth"
              set is_shedding 1
              set DemStat "dsResident"
              set hasActiveStrain true
              let tempTickLeth 0
              if fixedSurvivalTimes = false
              [
                set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)
                let Vir item (strainVirToGet - 1) survivalTimeList
                if strain < -10 [set vir item 0 survivalTimeList]
                set TickLeth round(tempTickLeth + Vir)
              ]
              if fixedSurvivalTimes = true
              [
                set TickLeth get_survival_time strain
              ]
              set NewLeth NewLeth + 1
              set LethPeriodList fput (TickLeth - ticks) LethPeriodList
              set InfPeriodList fput (TickLeth - ticks) InfPeriodList
            ]
          ] ]
          [
            ;; if the breeder is not infected
            hatch 8
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
        ] ] ]

        if (RandomNo > 0.999575749 AND RandomNo <= 0.99997575)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse is_shedding = 1 AND hasActiveStrain = true
          [
            hatch floor (9 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
             if strTra != 99999
            [
              ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
              set Age 0
              set AgeGroup "Piglet"
              set is_breeder 0
              set strain strTra
              set strainVirToGet strVirG
              set EpiStat "esLeth"
              set is_shedding 1
              set DemStat "dsResident"
              set hasActiveStrain true
              let tempTickLeth 0
              if fixedSurvivalTimes = false
              [
                set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)
                let Vir item (strainVirToGet - 1) survivalTimeList
                if strain < -10 [set vir item 0 survivalTimeList]
                set TickLeth round(tempTickLeth + Vir)
              ]
              if fixedSurvivalTimes = true
              [
                set TickLeth get_survival_time strain
              ]
              set NewLeth NewLeth + 1
              set LethPeriodList fput (TickLeth - ticks) LethPeriodList
              set InfPeriodList fput (TickLeth - ticks) InfPeriodList
            ]
          ] ]
          [
            ;; if the breeder is not infected
            hatch 9
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
        ] ] ]

        if (RandomNo > 0.99997575 AND RandomNo <= 1.0)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse is_shedding = 1 AND hasActiveStrain = true
          [
            hatch floor (10 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
             if strTra != 99999
            [
              ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
              set Age 0
              set AgeGroup "Piglet"
              set is_breeder 0
              set strain strTra
              set strainVirToGet strVirG
              set EpiStat "esLeth"
              set is_shedding 1
              set DemStat "dsResident"
              set hasActiveStrain true
              let tempTickLeth 0
              if fixedSurvivalTimes = false
              [
                set tempTickLeth 1 + floor (- (mue - 0.5) * ln random-float 1)
                let Vir item (strainVirToGet - 1) survivalTimeList
                if strain < -10 [set vir item 0 survivalTimeList]
                set TickLeth round(tempTickLeth + Vir)
              ]
              if fixedSurvivalTimes = true
              [
                set TickLeth get_survival_time strain
              ]
              set NewLeth NewLeth + 1
              set LethPeriodList fput (TickLeth - ticks) LethPeriodList
              set InfPeriodList fput (TickLeth - ticks) InfPeriodList
            ]
          ] ]
          [
            ;; if the breeder is not infected
            hatch 10
            [
              set DemStat "dsOffspring"
              set family who
            set float_counter 0
     set hasActiveStrain false
    set initialOutbreak false
    set transientCheck false
    set markedForInfection false
    set strainVirToGet 0
    set counterGotInfected false
            set strain 999
        ] ] ]

        set TickBirth ticks
        set is_breeder 0
        set is_mother 1
  ]
  set strTra 99999
  set strVirG 99999
  ]


  ask turtles with [DemStat = "dsOffspring"]
  [
    ifelse random-float 1 <= ProbFem [ set is_female 1 ] [ set is_female 0 ]
    set Age 0
    set AgeGroup "Piglet"
    set is_breeder 0

    ;; Susceptibles
    if [EpiStat] of turtle family = "esSusc"
    [
      set EpiStat "esSusc"
      set is_shedding 0
      set DemStat "dsResident"
    ]

    ;; Immunes
    if [EpiStat] of turtle family = "esImm"
    [
      set EpiStat "esImmMat"
      set is_shedding 0
      set DemStat "dsResident"
      set is_matAB 1
      set hasActiveStrain false
    ]
  ]

end

;; DEATH ______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Iterating over the entire population, each individual either stochastically dies due to age class dependent mortality rates by drawing a random number,
  ;;    due to reaching a certain maximum age or due to a lethal infection after a certain ifnection time span (TickLeth).
  ;; Stochastic baseline mortality is age-dependent and are drawn from a Gaussian distribution which was adjusted to annual survival estimates found in the literature and is cutted symmetrically around the mean.
  ;; The stochastic effect resembles 'good' versus 'bad' years for wild boars (e.g. environmental noise).
  ;; Each time step the adjusted age-dependent mortality (MortAge) is applied to all individual.


to Host_Death

  ;; survival probability for subadults and adults
  set SurvivalProb 0.65
  ;; survival probability for piglets
  set SurvivalProbPig 0.5

  ask turtles
  [
    ifelse Age > 34
    [
      if random-float 1 < (1 - (SurvivalProb ^ (1 / 52))) [ die ]
    ][  ;; mean survival probability for yearlings (Gaillard et al. 1987) = 0.65 = mean survival probability for adults (Focardi et al. 1996)
      if random-float 1 < (1 - (SurvivalProbPig ^ (1 / 52))) [ die ]
  ] ]

  ;; mortality due to a lethal infection
  ask turtles with [EpiStat = "esLeth"]
  [
    if TickLeth = ticks [ die ]
  ]


  ;; mortality due to reaching maximum age (more than 10 years)
  ask turtles with [Age > MaxAge] [ die ]

  ask turtles with [DemStat = "floater" AND float_counter >= survival_weeks]
  [
    if (quality_memory / float_counter) <= 4.5
    [
     die
    ]
    if (quality_memory / float_counter) > 4.5
    [
      if random 100 < 90 [die]
      ;[set DemStat  "dsResident"]
    ]
  ]

  ask turtles with [DemStat = "floater_piglet" AND float_counter >= 2] [if random 100 < 75 [die]]
   ask turtles with [DemStat = "floater_piglet" AND float_counter >= 3] [if random 100 < 85 [die]]
   ask turtles with [DemStat = "floater_piglet" AND float_counter >= 4] [die]

  ask Habitat with [not any? turtles-here AND is_occupied = 1] [ set is_occupied 0 ]

  ;;OUTPUT
  set PopStructList []
  ask turtles
  [ set PopStructList fput Age PopStructList ]
  set PopStructList map [ ?1 -> ?1 / 52 ] PopStructList

end

;; AGE ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; For each individual age is incremented one week and disease state transitions are performed.
  ;; Females who gave birth more than 33 weeks before can reproduce again.

to Host_Ageing

  ask turtles
  [
    set Age Age + 1

    ifelse Age = 35
    [
      set AgeGroup "Subadult"
      if is_female = 0 [ set DemStat "dsDispM_subadult" ]
    ]
    [
      ifelse Age = 53 AND is_female = 1
      [ set AgeGroup "Adult" ]
      [
        if Age = 105 AND is_female = 0
        [
          set AgeGroup "Adult"
          if Qual_around > 0 [ set DemStat "dsDispM_adult" ]  ;patches-own
    ] ] ]

    if is_mother = 1 AND (TickBirth + 34) = ticks [ set is_mother 0 ]
  ]

end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  HOST MOVEMENT PROCEDURES                                                                                                                                                                                                                                 ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; HERD SPLITTING _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; During summer, subadult females either remain with their mothers or establish new territories nearby.
  ;; Herd splitting is performed in one week per year only where female subadults (8 months = 34 weeks < Age <= 2 years = 104 weeks) without offspring are able to move between herds.
  ;; All herds to split are extracted matching the conditions of containing a number of females exceeding the cells carrying capacity and containing at least 2 dispersers.
  ;; For each splittable herd an habitat cell without female groups and a carrying capacity above 0 within the dispersal distance of 6 km (3 cells) is selected randomly. All dispersers from the source cell find a new herd on that cell.


to Movement_Dispersal-Females

  let OverCap Habitat with [count turtles-here with [is_female = 1 AND Age > 34] - counter_d >= 2 AND Qual_around > 0]

  ask OverCap
  [
    ;; count females above breeding capacity (Quality)
    let num-overfem ceiling (count turtles-here with [is_female = 1 AND Age > 34] - counter_d)
    ;; ... assign those females a disperser-status...
    let num-disp count turtles-here with [AgeGroup = "Subadult" AND is_female = 1 AND is_breeder = 0]
    ;; ... which are over breeding capacity ...
    if num-disp > num-overfem [ set num-disp num-overfem ]
    ;; ... and if there are at least 2 dispersers ...
    if num-disp >= 2
    [ ask n-of num-disp turtles-here with [AgeGroup = "Subadult" AND is_female = 1 AND is_breeder = 0] [ set DemStat "dsDispF" ] ]
    set NoDispF 0
  ]

  ask OverCap with [any? turtles-here with [DemStat = "dsDispF"]]
  [

    let ListDispF []
    ask turtles-here with [DemStat = "dsDispF"] [ set ListDispF fput self ListDispF ]
    set NoDispF length ListDispF   ; number of dispersing females of this cell
    let tempDispCells DispCells with [not any? turtles-here with [is_female = 1]]
    ifelse count tempDispCells > 0
    [
      let newHabitat one-of tempDispCells
      while [not empty? ListDispF]
      [
        ask first ListDispF
        [
          move-to newhabitat
          let cHerdID HerdID
          set HerdID currHerdID
          set DemStat "esResident"
          if is_shedding = 1
          [
            set is_infected 1 ;patches-own
            if distancexy InfX InfY > InfDist [ set InfDist distancexy InfX InfY ] ;patches-own
          ]
          set ListDispF remove-item 0 ListDispF
        ]
      ]
      ask newhabitat [ set is_occupied 1 ]
      set ListDispF []
    ]
    ;; else = if there are no unoccupied habitat cells within their dispersal distance the dispersers emmigrate/die
    [
      ask turtles-here with [DemStat = "dsDispF"]
      [
        ;; version used in Fernandez et al. 2006: dispersers stay in their maternal cell if no empty habitat cell is available
        set DemStat "esResident"
    ] ]
    set currHerdID currHerdID + 1
  ]

end

;;; MALE DISPERSAL _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Males perform their natal dispersal in an age of 8 months (35 weeks) upwards. All dispersing males form a group of at least 3 subadults and move to a cells below capacity with the highest habitat quality in dispersal distance.

to Movement_Dispersal-Males

  let NatalMales 0

  ask Habitat with [any? turtles-here with [DemStat = "dsDispM_subadult"] AND Qual_around > 0]
  [
    set NatalMales turtles-here with [DemStat = "dsDispM_subadult"]
    if count NatalMales > 2
    [
      let NewHome DispCells with [count turtles-here < Capacity]
      if count NewHome > 0
      [
        set NewHome max-one-of NewHome [counter_d]
        ask NatalMales
        [
          move-to NewHome
          set DemStat "dsResident"
  ] ] ] ]

end

;;; resource response _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to r_move  ;response to changing habitat conditions similar to herd-splitting

  ask Habitat with [count turtles-here  > Capacity ]
  [

    if count turtles-here with [DemStat = "dsResident"] > Capacity
    [
    while [count turtles-here with [DemStat = "dsResident"] > Capacity]
    [

      if any? turtles-here with [is_mother = FALSE OR AgeGroup != "Piglet"] AND count turtles-here  > Capacity
      [
      ask one-of turtles-here with [is_mother = FALSE OR AgeGroup != "Piglet"] [set DemStat "mover" ]
      ]

      if any? turtles-here with [AgeGroup = "Piglet"] AND count turtles-here  > Capacity
      [
      ask one-of turtles-here with [ AgeGroup = "Piglet"][set DemStat "mfd" ]
      ]

      if any? turtles-here with [DemStat = "dsResident"] AND count turtles-here  > Capacity
      [
        ask one-of turtles-here [set DemStat "mover" ]
      ]

    ]
    ]

   let NewHome2 DispCells with [is_occupied = FALSE]

    ask turtles-here with [DemStat = "floater"]
    [
      if count NewHome2 > 0
      [
        let NewHome3 one-of NewHome2

        ask turtles-here with [DemStat = "mover"]
        [
          move-to NewHome3
          set DemStat "dsResident"
        ]
      ]
      if count NewHome2 <= 0
      [
        set quality_memory quality_memory + counter_d
        set DemStat  "floater"
        set float_counter  float_counter + 1
      ]
    ]

    ask turtles-here with [DemStat = "floater_piglet"]
    [
      set quality_memory quality_memory + counter_d
      set float_counter  float_counter + 1
    ]

    ask turtles-here with [DemStat = "mover"]
    [
     if count NewHome2 > 0
      [
        let NewHome3 one-of NewHome2

        ask turtles-here with [DemStat = "mover"]
        [
          move-to NewHome3
          set DemStat "dsResident"
        ]
      ]
      if count NewHome2 <= 0
      [
        set quality_memory counter_d
        set DemStat  "floater"
        set float_counter  float_counter + 1
      ]
    ]

    ask turtles-here with [DemStat = "mfd"]
    [
     set quality_memory counter_d
     set DemStat  "floater_piglet"
     set float_counter  float_counter + 1
    ]
  ]

end

to check_floaters

  ask Habitat with [any? turtles-here with [DemStat = "floater"] OR any? turtles-here with [DemStat = "floater_piglet"]]
  [
    if count turtles-here < Capacity
      [
        ask turtles-here with [DemStat = "floater" OR DemStat = "floater_piglet"]
        [
          set DemStat "dsResident"
          set float_counter 0
          set quality_memory 0
        ]
      ]

  ]

end

;; ROAMING BEHAVIOUR __________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Long-distance movement of males older than 2 years (104 weeks) following the 'ranging strategy'.

  ; The staying strategy consists of remaining within a limited area, the activity zone, and moving short distances at different speeds. This strategy is mostly used by females accompanied by piglets, for activities such as
  ;    resting or social interaction (Morelle et al. 2015).
  ; The ranging strategy, usually employed by males at night, consists of ranging outside the activity zone (Morelle et al. 2015).
  ;
  ; In a rich habitat that provides easy access to food all year and safe resting sites, wild boar restrict their movement to a smaller area and thus reduce energy expenditure (Massei et al. 1997), in accordance with the food-exploitation hypothesis.
  ; In poor nutritional conditions, wild boar move more in search of food and water, consequently increasing their home range (Caley 1997, Massei et al. 1997).
  ;
  ; When population density is high and food availability low, competition increases, and wild boar increase their movement rate to search and compete for food (Bertolotto 2010).
  ;
  ; If ready to reproduce, males travel long distances in search of a sounder of sows, eating little on the way.
  ; Once a sounder has been located, the male drives off all young animals and persistently chases the sows.
  ;
  ; Males begin participating in the rut after 4–5 years.
  ;
  ; male dispersal distances (Morelle et al. 2015):
  ;   - daily range covered = 45% of annual range (1.3 to 2.4 km^2)
  ;   - mean daily distance travelled = 7.2 to 11.4 km
  ;   - daily moved distance = 3 to 4 km with a maximum of 12 km  ->  per week: 21-28 km (10.5-14 cells), max. 84 km (42 cells)

to Movement_Roaming

  if Roaming != "OFF"
  [
    ask patches
    [
      set NoInfectedMove count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain > -5 AND strain < 5]
      set NoInfectedMoveH count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain >= 5 ]
      set NoInfectedMoveL count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain <= -5]
    ]

    if Roaming = "DD"
    [
      ask Habitat
      [
        set W_M (1 - abs((Capacity - count turtles-here - 1) / Capacity))
        if W_M > 1 [ set W_M 1 ]
        if W_M < -1 [ set W_M -1 ]
      ]
    ]

    if Roaming = "FD"
    [
      ask Habitat [ set NoMale count turtles-here with [DemStat = "dsDispM_adult"] ]
      ask patches
      [
        ifelse (any? turtles-here with [DemStat = "dsDispM_adult"] OR NoRepro > 0)
        [ set W_M ((NoRepro - count turtles-here with [DemStat = "dsDispM_adult"]) / (NoRepro + count turtles-here with [DemStat = "dsDispM_adult"])) ]
        [;if no reproductive males or females
          ifelse counter_d >= 0
          [ set W_M -1 ]
          [ set W_M -1.1 ]
    ] ] ]

    if Roaming = "HDD"
    [
      ask Habitat
      [
        ifelse (any? turtles-here)
        [ ;if no turtles
          if Capacity != 0
          [
          set W_M (1 - count turtles-here / Capacity)
          ]
           if W_M > 1 [ set W_M 1 ]
          if W_M < -1 [ set W_M -1 ]
        ][
          ifelse counter_d > 0
          [ set W_M -1 ]
          [ set W_M -1.1 ]
        ]
        ask Matrix [ set W_M -1.1 ]
    ] ]



    ask turtles with [DemStat = "dsDispM_adult" AND Qual_around > 0]
    [
      let Dist 999   ; arbitrary high value
      while [Dist > 42] [ set Dist (Report_weibull 1.3 26) + 1 ]
      ;                         ---CONTROL--- print Dist

      set IndBetaMove (BetaMove / Dist) ;;betamove 0.0224 ADJUSTED for Viru
      set IndBetaMoveHigh ((BetaMoveHV) / Dist)
      set IndBetaMoveLow ((BetaMoveLV) / Dist)
      ;                         ---CONTROL--- print IndBetaMove
     if movement_virulence_tradeoff = true
      [
        let tmpD Dist
        let tmpND dist_vir_tradeoff tmpD
        set Dist tmpND
      ]

      set MoveDist Dist
  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 1 (random walk)
      if Roaming = "RW"
      [
        while [Dist > 0]
        [
          Disease_Infection-on-the-Move
          move-to one-of neighbors with [is_habitat = 1]
          set Dist Dist - 1
      ] ]
  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 2 (correlated random walk)
      if Roaming = "CRW"
      [
        while [Dist > 0]
        [
          Disease_Infection-on-the-Move
          ;; wrapped cauchy distribution (equation from Haefner & Crist 1994, also used by Zollner & Lima 1999, Fletcher 2006)
          rt (2 * (atan (( (1 - q) / (1 + q)) * tan ((random-float 1 - 0.5) * 180)) 1))   ;; produces values between 0-180 and 540-720 if q = 0
          while [ patch-ahead 1 = nobody OR member? patch-ahead 1 Matrix ] [ rt random 180 - 160 ]
          move-to patch-ahead 1
          set Dist Dist - 1
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 3 (habitat-dependent walk)
      if Roaming = "HD"
      [
        while [Dist > 0]
        [
          Disease_Infection-on-the-Move
          if not any? neighbors with [counter_d > 0] [ stop ]
          ifelse random-float 1 < q [ move-to one-of neighbors with-max [floor counter_d] ]
                                    [ move-to one-of neighbors with [is_habitat = 1] ]
          set Dist Dist - 1
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 4 (density-dependent walk)
      if Roaming = "DD"
      [
        while [Dist > 0]
        [
          Disease_Infection-on-the-Move
          let Aim neighbors with [counter_d > 0]
          if count Aim = 0 [ stop ]
          ifelse random-float 1 < q [ move-to max-one-of Aim [W_M] ]
                                    [ move-to one-of Aim ]
          set Dist Dist - 1
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 5 (female-dependent walk)
      if Roaming = "FD"
      [
        while [Dist > 0]
        [
          Disease_Infection-on-the-Move
          let Aim neighbors with [counter_d > 0]
          if count Aim = 0 [ stop ]
          ifelse random-float 1 < q   ; quality of decision
          [
            ifelse any? Aim with [NoRepro > 0]
            [ move-to max-one-of Aim with [NoRepro > 0] [W_M] ]
            [ move-to min-one-of Aim [NoMale] ]
          ]
          [ move-to one-of Aim ]
          set Dist Dist - 1
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 6 (habitat-density-dependent walk)
      if Roaming = "HDD"
      [
        while [Dist > 0]
        [
          Disease_Infection-on-the-Move
          let Aim neighbors with [counter_d > 0]
          if count Aim = 0 [ stop ]
          ifelse random-float 1 < q [ move-to max-one-of Aim [W_M] ]
                                    [ move-to one-of Aim ]
          set Dist Dist - 1
      ] ]
  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 7 (habitat-dependent correlated random walk)
      if Roaming = "HD-CRW"
      [
        while [Dist > 0]
        [

          ;; wrapped cauchy distribution (equation from Haefner & Crist 1994, also used by Zollner & Lima 1999, Fletcher 2006)
          ifelse random-float 1 < q
          [
            ;; CRW for long-range movements only - fix q because local movement arises from habitat-dependent movement
            rt (2 * (atan (( (1 - 0.9) / (1 + 0.9)) * tan ((random-float 1 - 0.5) * 180)) 1))   ;; produces values between 0-180 and 540-720 if q = 0
            let resc_counter 0
            while [ patch-ahead 1 = nobody OR member? patch-ahead 1 Matrix AND resc_counter < 50 ]
            [
              rt random 90 - 45
             set resc_counter resc_counter + 1
            ]
            ifelse patch-ahead 1 = nobody OR member? patch-ahead 1 Matrix
            [
              set Dist 0
            ]
            [
            move-to patch-ahead 1
            ]
          ]
          [
            if not any? neighbors with [counter_d > 0] [ stop ]
            move-to one-of neighbors with-max [floor counter_d]
          ]
          set Dist Dist - 1
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; AFTERWARDS
    ]
    let LethInd turtles with [EpiStat = "esLeth"]
    set CurrInfPeriodList []
    ask LethInd [ set CurrInfPeriodList fput (TickLeth - ticks) CurrInfPeriodList ]
  ]

  if Roaming != "OFF" [ set MeanMoveDist mean [MoveDist] of turtles with [DemStat = "dsDispM_adult"] ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; HELP FUNCTIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; WEIBULL DISTRIBUTION REPORTER ______________________________________________________________________________________________________________________________________________________________________________________________________________________________

to-report Report_weibull [a-scale a-shape]

  let xWei random-float 1
  let yWei (1 / a-scale)
  let zWei ((a-shape * (-1 * ln(1 - xWei)))^ yWei)
  report zWei

end

;; FRAGMENTAITON INDEX by Mitchell & Powell 2004 ______________________________________________________________________________________________________________________________________________________________________________________________________________

to Report_fragmentation

  set F_infected 0
  set F_infectious 0

  let Nmax_infected count patches with [is_infected = 1]
  let Nmax_infectious count patches with [is_infected = 1]

  set Nmax_infected item Nmax_infected NmaxList
  set Nmax_infectious item Nmax_infectious NmaxList

  if Nmax_infected > 0 [ set F_infected (1 - ((sum ([count neighbors with [is_infected = 1]] of patches with [is_infected = 1])) / Nmax_infected )) ]
  if Nmax_infectious > 0 [ set F_infectious (1 - ((sum ([count neighbors with [any? turtles-here with [is_shedding = 1]]] of patches with [any? turtles-here with [is_shedding = 1]])) / Nmax_infectious )) ]

end

;; UPDATE VIEW ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to View

  ;ask turtles [ st ]
  ask turtles with [EpiStat = "esSusc"] [ set color black ]
  ask turtles with [EpiStat = "esImm" OR EpiStat = "esImmMat"] [ set color green ]
  ask turtles with [EpiStat = "esLeth"]
  [
    set color red
    set size 0.75
  ]
  ask turtles with [EpiStat = "esTrans"]
  [
    set color yellow
    set size 0.75
  ]
   ifelse dynmode_visuals = true
  [
    ask patches [ set pcolor scale-color grey counter_d 0 12 ]
  ]

  [
     ifelse landscapeGscale = true
    [
      ask patches [ set pcolor scale-color grey quality 0 12 ]
    ]
    [
      ask patches [ set pcolor white]
    ]
  ]
  if evolution = false
  [
  ask Habitat with [any? turtles-here with [is_shedding = 1]] [ set pcolor scale-color red counter_d 12 0 ]
  ]

  scheme_color_scale_gradient

  if ticks >= WeekRelease - 1 [ ask patch InfX InfY [ set pcolor violet ] ]

end




;; UPDATE _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to Update
  ;; update numbers
  ask turtles
  [
    set NoGroup count turtles-here
    set NoInfected count turtles-here with [is_shedding = 1]
  ]
  ask Habitat [
    set NoBoars count turtles-here
    set tmpToInfectHere 0
  ]

  ;; update infection states patches
  ask Habitat [ set is_infectious 0 ]
  ask Habitat with [NoInfected > 0] [ set is_infectious 1 ]

  ;; stop condition
  if (ticks >= WeekRelease - 1) AND ((not any? turtles with [is_shedding = 1]) OR (ticks = SimulatedYears * 52)) [ set DONE 1 ]
  if DONE = 1 [ stop ]

  ;; OUTPUT
  if (ticks >= WeekRelease AND length LethPeriodList > 0)
  [
    set MeanLethPeriod mean(LethPeriodList)
    set MaxLethPeriod max(LethPeriodList)
  ]
  set WeekLast ticks + 1
  if MaxDist = 0 AND (any? patches with [pycor = 0 AND is_infectious = 1]) [ set MaxDist ticks + 1 ]

  if any? turtles with [hasActiveStrain = true AND strain != "none"]
  [
 set meanStrainOutput mean[strain] of turtles with [hasActiveStrain = true AND strain != "none"]
  ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions for landscape dynmaics
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Landscape buffer _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
to buffer

  if (dynmode = TRUE)
    [

  let counter_mean1  mean [counter_d] of patches
  let controll_mean1 floor mean [Quality] of patches

if (floor controll_mean1 != floor counter_mean1 AND week > 10)[
    while [floor controll_mean1 != floor counter_mean1 AND week > 10]
  [
    if (floor controll_mean1 <= floor counter_mean1)
      [
    ask one-of patches with [inc = true AND counter_d != sq]
    [
    set counter_d counter_d - 1
    set counter_mean1 floor mean [counter_d] of patches
    ]
      ]
        if (floor controll_mean1 >= floor counter_mean1)
      [
        ask one-of patches with [dec = true AND counter_d != sq]
    [
    set counter_d counter_d + 1
    set counter_mean1 floor mean [counter_d] of patches
    ]
      ]
     if (floor controll_mean1 = floor counter_mean1) [stop]

  ]
  ]

    ]

end

;; Landscape dynamic _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to dyn

  ask patches
  [
    let qb mean [counter_d] of patches
    set qmw qb
     if (Quality < round 5 AND dynamic != true  )
    [
      set sq round Quality
      set counter_d round Quality
      set qmin Quality
      set dynamic true
      set inc true
      set dec false
    ]

     if (Quality > round 5 AND dynamic != true )
    [
     set sq round Quality
     set counter_d round Quality
     set qmax Quality
     set dynamic true
     set dec true
     set inc false
    ]

    set Qual_around sum [counter_d] of neighbors

    ;; temporal lag was seperated into individual functions for each level of lag ( 0 - 100 in 25 increments) for ease of access


    if (qmw != 0)
 [
 if (dynmode = true AND noise_type = "Red Noise" AND experimental_rednoise = false AND temporal_lag = 100)
 [
      if (week = 27 )
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
         if (week = 33 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 37 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 42 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 47 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 3 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
          if (week = 7 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 12 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 17 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 22 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 1)[

          set counter_d round Quality
          set Capacity (counter_d * (20 / qmw))
        ]

 ]
    ]




      if (qmw != 0)
 [
 if (dynmode = true AND noise_type = "Red Noise" AND experimental_rednoise = false AND temporal_lag = 0)
 [
      if (week = 5 + dyntimer_up - dyntimer_down + week_overlap AND counter_d < 9)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
         if (week = 10 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 15 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 20 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 25 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 30 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
          if (week = 35 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 40 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 45 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 50 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 1)[

          set counter_d round Quality
          set Capacity (counter_d * (20 / qmw))
        ]

 ]
    ]

    if (qmw != 0)
 [
if (dynmode = true AND noise_type = "Red Noise" AND experimental_rednoise = false AND temporal_lag = 25)
 [
      if (week = 10 )
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
         if (week = 15 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 20 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 25 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 30 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 35 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
          if (week = 40 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 45 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 50 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 5 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 5)[
          set counter_d round Quality
          set Capacity (counter_d * (20 / qmw))
        ]
      ]
 ]


    if (qmw != 0)
 [
if (dynmode = true AND noise_type = "Red Noise" AND experimental_rednoise = false AND temporal_lag = 50)
 [
      if (week = 15 )
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
         if (week = 20 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 25 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 30 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 35 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 40 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
          if (week = 45 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 50 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 5 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 10 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 10)[
          set counter_d round Quality
          set Capacity (counter_d * (20 / qmw))
        ]
      ]
 ]




    if (qmw != 0)
 [
if (dynmode = true AND noise_type = "Red Noise" AND experimental_rednoise = false AND temporal_lag = 75)
 [
      if (week = 20 )
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
         if (week = 25 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 30 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 35 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 40 + dyntimer_up - dyntimer_down AND counter_d < 9 AND counter_d != 0)
      [
        set counter_d counter_d + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 45 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
          if (week = 50 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 5 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 10 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 15 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d counter_d - 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 15)[
          set counter_d round Quality
          set Capacity (counter_d * (20 / qmw))
        ]
      ]
 ]

    if (qmw != 0)
 [
    if (dynmode = true AND noise_type = "White Noise")
    [
    if (week = 5 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d random (9) + 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 10 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d random (9) + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 15 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d random (9) + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 20 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d random (9) + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 25 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d random (9) + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 30 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d random (9) + 1
        set Capacity (counter_d * (20 / qmw))
      ]
       if (week = 35 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d random (9) + 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 40 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d random (9) + 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 45 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d random (9) + 1
        set Capacity (counter_d * (20 / qmw))
      ]
        if (week = 50 + dyntimer_up - dyntimer_down AND counter_d != 0)
      [
        set counter_d random (9) + 1
        set Capacity (counter_d * (20 / qmw))
      ]
    ]
    ]

    if (qmw != 0)
    [
    if (dynmode = true AND noise_type = "Red Noise" AND experimental_rednoise = true)
    [
     let counter_li ticks


        let helper15 item counter_li rednoise_list
        set counter_d ( round helper15)

        if counter_d > 9
        [
          set counter_d 9
        ]
        if counter_d <= 0
        [
          set counter_d 0
        ]

        set Capacity (counter_d * (20 / qmw))

    ]
    ]
  ]

end

;; Landscape reset _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to resetdyn

;; landscape is reset to the starting landscape after each year ( 52 simulation ticks)

  if (dynmode = TRUE)
  [
  ask patches
  [
   if (counter_d >= round 9 AND dynamic = true)
    [
      set inc false
      set dec true

    ]

    if (counter_d <= round 0 AND dynamic = true)
    [
      set inc true
      set dec false

    ]

    if (counter_d = sq AND week > 10 AND dynamic = true)
    [
      if (counter_d <= round 5 AND dynamic = true)
      [
      set inc true
      set dec false
      set counter_d counter_d + 1

      if qmw != 0
        [
      set Capacity (counter_d * (20 / qmw))
        ]
      ]

      if (counter_d > round 5 AND dynamic = true)
      [
      set inc false
      set dec true
      set counter_d counter_d - 1
      set Capacity (counter_d * (20 / qmw))
      ]
    ]
    ]
  ]
end

;; Output spatial distribution _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


to output_list

  set cells_in_order patch-set sort patches

  ask patches
  [
    if (any? turtles-here with [is_shedding = 1]) [set infection_output_counter infection_output_counter + 1]
  ]
  ask patches
  [
    if(DONE = 1 or ticks = 2199)
    [
      set o_list fput (list pxcor";"pycor";" infection_output_counter ";" counter_d) o_list
    ]
  ]
end

;; Evolution

to getStrain
              ;; called within a ask turtle
  let transferStrain 9999 ;for control since strain will be between -50 and 50
  let r 1
  set sumOfInfectionEvents sumOfInfectionEvents + 1

  if evolution = false
  [
    set strain 0
  ]

  if evolution = true
  [
    ;strainVirToGet ; determines which strain it got its infection from 1 low 2 mid 3 high
    if strainTransmissionType = "exp"
    [
      if strainVirToGet = 13
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 10 AND strain <= 50 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 10 AND strain <= 50    ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 10 AND strain <= 50   ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 10 AND strain <= 50    ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure13"
            ]
          ]
        ]
      ]

      if strainVirToGet = 12
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 8 AND strain <= 10 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 8 AND strain <= 10    ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 8 AND strain <= 10   ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 8 AND strain <= 10    ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure12"
            ]
          ]
        ]
      ]

      if strainVirToGet = 11
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 6 AND strain <= 8 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 6 AND strain <= 8    ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 6 AND strain <= 8   ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 6 AND strain <= 8    ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure11"
            ]
          ]
        ]
      ]

      if strainVirToGet = 10
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 4 AND strain <= 6 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 4 AND strain <= 6   ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 4 AND strain <= 6  ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 4 AND strain <= 6   ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure10"
            ]
          ]
        ]
      ]

      if strainVirToGet = 9
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 2 AND strain <= 4 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 2 AND strain <= 4  ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 2 AND strain <= 4  ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 2 AND strain <= 4  ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure9"
            ]
          ]
        ]
      ]

      if strainVirToGet = 8
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 0 AND strain <= 2 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 0 AND strain <= 2  ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 0 AND strain <= 2  ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > 0 AND strain <= 2  ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure8"
            ]
          ]
        ]
      ]

      if strainVirToGet = 7
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain = -0 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain = 0  ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain = -0  ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain = -0  ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure7"
            ]
          ]
        ]
      ]


      if strainVirToGet = 6
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > -2 AND strain < 0 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > -2 AND strain < 0  ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > -2 AND strain < 0  ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain > -2 AND strain < 0  ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure6"
            ]
          ]
        ]
      ]

      if strainVirToGet = 5
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -2 AND strain >= -4 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -2 AND strain >= -4 ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -2 AND strain >= -4 ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -2 AND strain >= -4 ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure5"
            ]
          ]
        ]
      ]

      if strainVirToGet = 4
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -4 AND strain >= -6 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -4 AND strain >= -6 ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -4 AND strain >= -6 ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -4 AND strain >= -6 ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure4"
            ]
          ]
        ]
      ]


      if strainVirToGet = 3
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -6 AND strain >= -8 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -6 AND strain >= -8 ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -6 AND strain >= -8 ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -6 AND strain >= -8 ]
              [
                set transferStrain strain
              ]
            ]
            [
              ;print "failure3"
            ]
          ]
        ]
      ]

      if strainVirToGet = 2
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -8 AND strain >= -10 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -8 AND strain >= -10 ]
          [
            set transferStrain strain
            ;print (word "strainset1:"  strain) ;; DEBUG remove me
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -8 AND strain >= -10 ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -8 AND strain >= -10 ]
              [
                set transferStrain strain
                ; print (word "strainset2:"  transferStrain) ;;DEBUG
              ]
            ]
            [
              ;print "failure2"
            ]
          ]
        ]
      ]



      if strainVirToGet = 1
      [
        ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -10 ]
        [
          ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -10 ]
          [
            set transferStrain strain
          ]
        ]
        [
          ask neighbors
          [
            ifelse any? other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -10 ]
            [
              ask one-of other turtles-here with [hasActiveStrain = true AND strain != 999 AND strain < -10 ]
              [
                set transferStrain strain
              ]
            ]
            [
             ; print "failure1"
            ]
          ]
        ]
      ]

      if strainVirToGet = 0 OR strainVirToGet > 13 OR strainVirToGet < -100 ;; FAILSAVE FIX ME
      [
        print "FAILSAFE"
        set strain (random 20 )- 10
      ]



      ;;;Failsave in case of mortalilty or movement removing individuals from the pool

      if transferStrain = 9999 AND strainVirToGet > 0 AND strainVirToGet < 14
      [
        ;print "Failsave"
        ifelse strainVirToGet = 1 [set transferstrain -10 set hasActiveStrain true][
          ifelse strainVirToGet = 2 [set transferstrain ((10 - (random 3)) * -1) set hasActiveStrain true][
            ifelse strainVirToGet = 3 [set transferstrain ((8 - (random 3)) * -1) set hasActiveStrain true][
              ifelse strainVirToGet = 4 [set transferstrain ((6 - (random 3)) * -1)set hasActiveStrain true][
                ifelse strainVirToGet = 5 [set transferstrain ((4 - (random 3)) * -1)set hasActiveStrain true][
                  ifelse strainVirToGet = 6 [set transferstrain ((2 - (random 3)) * -1)set hasActiveStrain true][
                    ifelse strainVirToGet = 7 [set transferstrain (0)set hasActiveStrain true][
                      ifelse strainVirToGet = 8 [set transferstrain (( (random 3)) )set hasActiveStrain true][
                        ifelse strainVirToGet = 9 [set transferstrain ((2 + (random 3)))set hasActiveStrain true][
                          ifelse strainVirToGet = 10 [set transferstrain ((4 + (random 3)))set hasActiveStrain true][
                            ifelse strainVirToGet = 11 [set transferstrain ((6 + (random 3)))set hasActiveStrain true][
                              ifelse strainVirToGet = 12 [set transferstrain ((8 + (random 3)))set hasActiveStrain true][
                                if strainVirToGet = 13 [set transferstrain (10)set hasActiveStrain true]
        ]]]]]]]]]]]]
      ]


      ifelse transferStrain != 9999
      [
        set strain transferStrain
        set hasActiveStrain true
      ]
      [
        ;user-message (word "ERROR 003 - transfer strain out of bounds - getStrain_function")
        set warningCounter warningCounter + 1
        print (word "Warning -" warningCounter "- strain transmission failure")
        set transmissionFailures transmissionFailures + 1
        print (word "Total number of transmission failures: " transmissionFailures)
        if transmissionFailures > 10 [
          user-message (word "ERROR 003 - Total number of transmission failures above acceptable threshold - getStrain_function")
        ]
      ]

      if strain < -50
      [set strain -50]
      if strain != 999 AND strain != 9999
      [
        if strain > 50
        [
          set strain 50

        ]
      ]
    ]

    if strainTransmissionType = "random"
    [
      if any? turtles-here with [hasActiveStrain = true]
      [
        ask one-of turtles with [hasActiveStrain = true]
        [
          set transferStrain strain
        ]
      ]

      set strain transferStrain
      set hasActiveStrain true
    ]
  ]

end


to evolve
let rng4 random-float 1
let rng5 random 1000000

  if rng4 < mutationRate
  [
;print "I evolved" ;DEBUG
    let tempStrain strain
    set sumOfEvolutionEvents sumOfEvolutionEvents + 1
    set strain random-normal tempStrain evolutionFactor
    ;; softcap for testing
    if strain < -15
    [set strain -10]
    if strain > 15
    [set strain 10]

  ]



end

to virulenceCheck
; removed

end



to Disease_Infection-on-the-Move
; not impemented in this version
  if is_shedding = 1 AND hasActiveStrain = true
  [

  let TmpMoveBetaList []
  let TmpVirList []

  let ipM10 (1 - (1 -  IndBetaMove) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain = 0 AND hasActiveStrain = true] ))
  let ipM9 (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain < 0 AND strain > -2 AND hasActiveStrain = true]))  ;;;  difference between each "brackets" beta
  let ipM8 (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain <= -2 AND strain > -4 AND hasActiveStrain = true]))
  let ipM7 (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain <= -4 AND strain > -6 AND hasActiveStrain = true]))
  let ipM6 (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain <= -6 AND strain > -8 AND hasActiveStrain = true]))
  let ipM5 (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain <= -8 AND strain > -10 AND hasActiveStrain = true]))
  let ipMMin (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain <= -10 AND hasActiveStrain = true]))
  let ipM11 (1 - (1 -  (IndBetaMove + (IndBetaMove * 0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain > 0 AND strain <= 2 AND hasActiveStrain = true]))
  let ipM12 (1 - (1 -  (IndBetaMove + (IndBetaMove * 0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain > 2 AND strain <= 4 AND hasActiveStrain = true]))
  let ipM13 (1 - (1 -  (IndBetaMove + (IndBetaMove * 0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain > 4 AND strain <= 6 AND hasActiveStrain = true]))
  let ipM14 (1 - (1 -  (IndBetaMove + (IndBetaMove * 0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain > 6 AND strain <= 8 AND hasActiveStrain = true]))
  let ipM15 (1 - (1 -  (IndBetaMove + (IndBetaMove * 0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain > 8 AND strain <= 10 AND hasActiveStrain = true]))
  let ipMMax (1 - (1 - (IndBetaMove + (IndBetaMove * 0))) ^ (count turtles-here with [DemStat = "dsDispM_adult" AND is_shedding = 1 AND strain > 10 AND hasActiveStrain = true])) ;ITEM 0

    if ipMMin != 0
    [
      set ipMMin ipMMin - (transmission_factor + 0.066 - lowVirOffset)
      if ipMMin <= 0 [set ipMMin 0.0001]
    ]

    set ipM10 ipM10

    if ipM9 != 0
    [
      set ipM9 ipM9 - (transmission_factor + 0.011 - lowVirOffset)
      if ipM9 <= 0 [set ipM9 0.0001]
    ]

    if ipM8 != 0
    [
      set ipM8 ipM8 - (transmission_factor + 0.022 - lowVirOffset)
      if ipM8 <= 0 [set ipM8 0.0001]
    ]

     if ipM7 != 0
    [
      set ipM7 ipM7 - (transmission_factor + 0.033 - lowVirOffset)
      if ipM7 <= 0 [set ipM7 0.0001]
    ]

      if ipM6 != 0
    [
      set ipM6 ipM6 - (transmission_factor + 0.044 - lowVirOffset)
      if ipM6 <= 0 [set ipM6 0.0001]
    ]

        if ipM5 != 0
    [
      set ipM5 ipM5 - (transmission_factor + 0.055 - lowVirOffset)
      if ipM5 <= 0 [set ipM5 0.0001]
    ]

    if ipM11 != 0
    [
    set ipM11 ipM11 + 0.011
    ]

    if ipM12 != 0
    [
    set ipM12 ipM12 + 0.022
    ]

    if ipM13 != 0
    [
    set ipM13 ipM13 + 0.033
    ]

    if ipM14 != 0
    [
     set ipM14 ipM14 + 0.044
    ]

    if ipM15 != 0
    [
    set ipM15 ipM15 + 0.055
    ]

    if ipMMax != 0
    [
    set ipMMax ipMMax + 0.066
    ]

  set TmpMoveBetaList lput ipMMin TmpMoveBetaList ;ITEM 0
  set TmpMoveBetaList lput ipM5 TmpMoveBetaList
  set TmpMoveBetaList lput ipM6 TmpMoveBetaList
  set TmpMoveBetaList lput ipM7 TmpMoveBetaList
  set TmpMoveBetaList lput ipM8 TmpMoveBetaList
  set TmpMoveBetaList lput ipM9 TmpMoveBetaList
  set TmpMoveBetaList lput ipM10 TmpMoveBetaList
  set TmpMoveBetaList lput ipM11 TmpMoveBetaList
  set TmpMoveBetaList lput ipM12 TmpMoveBetaList
  set TmpMoveBetaList lput ipM13 TmpMoveBetaList
  set TmpMoveBetaList lput ipM14 TmpMoveBetaList
  set TmpMoveBetaList lput ipM15 TmpMoveBetaList
  set TmpMoveBetaList lput ipMMax TmpMoveBetaList ; ITEM 15


    let tm item (strainVirToGet - 1) TmpMoveBetaList
    let tmi round [(tm * AllSus)] of patch-here
    let tStrain strainVirToGet


        ask at-most-n-of tmi other turtles-here with [markedForInfection = true]
          [
            if random-float 1 < tm
            [
              ; print "test" ;; DEBUG remove me
              set strainVirToGet tStrain ; determines which strain it got its infection from
              Disease_Infection

              ifelse strain != 999
              [
                set markedForInfection false
              ]
          [print "ERROR 004"]
            ]
          ]
  ]

  if EpiStat = "esSusc" AND any? other turtles-here with [is_shedding = 1 AND hasActiveStrain = true] ;(NoInfectedM > 0 OR NoInfectedMove > 0 OR NoInfectedHighVir > 0 OR NoInfectedLowVir > 0 OR NoInfectedMoveH > 0 OR NoInfectedMoveL > 0)
  [

  let ip10M (1 - (1 -  IndBetaMove) ^ (count turtles-here with [is_shedding = 1 AND strain = 0 AND hasActiveStrain = true] ))
  let ip9M (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.05))) ^ (count turtles-here with [is_shedding = 1 AND strain < 0 AND strain > -2 AND hasActiveStrain = true]))  ;;; 5% difference between each "brackets" beta
  let ip8M (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.10))) ^ (count turtles-here with [is_shedding = 1 AND strain <= -2 AND strain > -4 AND hasActiveStrain = true]))
  let ip7M (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.15))) ^ (count turtles-here with [is_shedding = 1 AND strain <= -4 AND strain > -6 AND hasActiveStrain = true]))
  let ip6M (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.20))) ^ (count turtles-here with [is_shedding = 1 AND strain <= -6 AND strain > -8 AND hasActiveStrain = true]))
  let ip5M (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.25))) ^ (count turtles-here with [is_shedding = 1 AND strain <= -8 AND strain > -10 AND hasActiveStrain = true]))
  let ipMinM (1 - (1 -  (IndBetaMove - (IndBetaMove * 0.30))) ^ (count turtles-here with [is_shedding = 1 AND strain <= -10 AND hasActiveStrain = true]))
  let ip11M (1 - (1 -  (IndBetaMove + (IndBetaMove * 0.05))) ^ (count turtles-here with [is_shedding = 1 AND strain > 0 AND strain <= 2 AND hasActiveStrain = true]))
  let ip12M (1 - (1 -  (IndBetaMove + (IndBetaMove * 0.10))) ^ (count turtles-here with [is_shedding = 1 AND strain > 2 AND strain <= 4 AND hasActiveStrain = true]))
  let ip13M (1 - (1 -  (IndBetaMove + (IndBetaMove * 0.15))) ^ (count turtles-here with [is_shedding = 1 AND strain > 4 AND strain <= 6 AND hasActiveStrain = true]))
  let ip14M (1 - (1 -  (IndBetaMove + (IndBetaMove * 0.20))) ^ (count turtles-here with [is_shedding = 1 AND strain > 6 AND strain <= 8 AND hasActiveStrain = true]))
  let ip15M (1 - (1 -  (IndBetaMove + (IndBetaMove * 0.25))) ^ (count turtles-here with [is_shedding = 1 AND strain > 8 AND strain <= 10 AND hasActiveStrain = true]))
  let ipMaxM (1 - (1 - (IndBetaMove + (IndBetaMove * 0.30))) ^ (count turtles-here with [is_shedding = 1 AND strain > 10 AND hasActiveStrain = true])) ;ITEM 0

    ifelse count other turtles-here with [is_shedding = 1 AND hasActiveStrain = true] > 0
    [

    let rng12 random-float 1

    ifelse rng12 <= [ipMinM] of patch-here AND rng12 > 0 [set strainVirToGet 1 Disease_Infection][
      ifelse rng12 <= [ip5M] of patch-here[set strainVirToGet 2 Disease_Infection][
        ifelse rng12 <= [ip6M] of patch-here[set strainVirToGet 3 Disease_Infection][
          ifelse rng12 <= [ip7M] of patch-here[set strainVirToGet 4 Disease_Infection][
            ifelse rng12 <= [ip8M] of patch-here[set strainVirToGet 5 Disease_Infection][
              ifelse rng12 <= [ip9M] of patch-here[set strainVirToGet 6 Disease_Infection][
                ifelse rng12 <= [ip10M] of patch-here[set strainVirToGet 7 Disease_Infection][
                  ifelse rng12 <= [ip11M] of patch-here[set strainVirToGet 8 Disease_Infection][
                    ifelse rng12 <= [ip12M] of patch-here[set strainVirToGet 9 Disease_Infection][
                      ifelse rng12 <= [ip13M] of patch-here[set strainVirToGet 10 Disease_Infection][
                        ifelse rng12 <= [ip14M] of patch-here[set strainVirToGet 11 Disease_Infection][
                          ifelse rng12 <= [ip15M] of patch-here[set strainVirToGet 12 Disease_Infection][
                            if rng12 <= [ipMaxM] of patch-here[set strainVirToGet 13 Disease_Infection]]]]]]]]]]]]]

    ]
    [
        ; not implemented in this version
    ]
  ]

end

to infectionCheck

  ask turtles with [EPiStat = "esSusc" ]
  [
    set markedForInfection true

  ]

end

to infectIndividuals

  ask habitat
  [
    set toInfect  round ((sum infPressList) * AllSus)
    let n 0
    if toInfect > 0
    [
      while [ n < 12 ]
      [  ;;; number of entrys in the list
        let ti item n infPressList
        set ti round (ti * AllSus)
        let si item n infPressList

        if any? turtles-here with [markedForInfection = true] ; check to not ask 0
        [
          ask at-most-n-of ti turtles-here with [markedForInfection = true]
          [
            if random-float 1 < si
            [
              set strainVirToGet item n virToGetList ; determines which strain it got its infection from
              Disease_Infection

              if strain != 999
              [
                set markedForInfection false
                ask patch-here
                [
                  set tmpToInfectHere tmpToInfectHere - 1
                ]
              ]
            ]
          ]
        ]
        set n n + 1
      ]
    ]
  ]


end



to invasion

  ask one-of habitat
  [
    ask turtles-here with [hasActiveStrain = false AND  EpiStat  = "esSusc"]
    [
      set strainVirToGet 12
      set strain 12
      set hasActiveStrain true
      Disease_Infection
    ]
  ]

end

to testcounter
  ask turtles with [hasActiveStrain = true AND strain > 500]
  [
    print strain

  ]

end

to-report at-most-n-of [ n agentset ]
  ifelse count agentset > n [
    report n-of n agentset
  ] [
    report agentset
  ]
end

to-report closestValue [valueToLookFor listToBrowse]
  let valueToLookForList n-values 13 [valueToLookFor] ; 13 catagories
  let diffList map abs(map -  valueToLookForList listToBrowse)
  let closestValueR item (position (min diffList) diffList) listToBrowse
  report (closestValueR)
end

to-report ListWithoutZero [inputList]
  let zList [0]
  let ipL inputList
  let result filter [ x -> not member? x zList ] ipL
  report result
end

to scheme_color_scale_gradient
  ask habitat  with [any? turtles-here with [is_shedding = 1 AND hasActiveStrain = true  ]]
  [
    set pcolor palette:scale-gradient palette:scheme-colors "Divergent" "Spectral" 9 max[strain * -1]  of turtles-here with [hasActiveStrain = true] -10 10
  ]
end



to-report dist_vir_tradeoff [Mdist]
  ; inside ask turtle
  let dist Mdist
  if hasActiveStrain = true
  [
    if strain > 5
    [
      let tmpS round (strain - 5)
      if tmpS > 0
      [
       let reductionPercent tmpS / 10

        set dist Mdist - (Mdist * reductionPercent)
        if dist <= 0
        [
          set dist 1
        ]
      ]
    ]
  ]
  report dist

end

to-report get_survival_time [vir]

  let surTime 11 ;;default value of "base"-strain

  ifelse vir < -10 [set surTime item 12 fixedExtTimesList][
    ifelse vir <= -9 [set surTime item 11 fixedExtTimesList][
      ifelse vir <= -8 [set surTime item 11 fixedExtTimesList][
        ifelse vir <= -7 [set surTime item 10 fixedExtTimesList][
          ifelse vir <= -6 [set surTime item 10 fixedExtTimesList][
            ifelse vir <= -5 [set surTime item 9 fixedExtTimesList][
              ifelse vir <= -4 [set surTime item 9 fixedExtTimesList][
                ifelse vir <= -3 [set surTime item 8 fixedExtTimesList][
                  ifelse vir <= -2 [set surTime item 8 fixedExtTimesList][
                    ifelse vir <= -1 [set surTime item 7 fixedExtTimesList][
                      ifelse vir <= 0 [set surTime item 7 fixedExtTimesList][
                        ifelse vir > 0 AND vir < 1 [set surTime item 5 fixedExtTimesList][
                          ifelse vir >= 1 AND vir < 2 [set surTime item 5 fixedExtTimesList][
                            ifelse vir >= 2 AND vir < 3 [set surTime item 4 fixedExtTimesList][
                              ifelse vir >= 3 AND vir < 4 [set surTime item 4 fixedExtTimesList][
                                ifelse vir >= 4 AND vir < 5 [set surTime item 3 fixedExtTimesList][
                                  ifelse vir >= 5 AND vir < 6 [set surTime item 3 fixedExtTimesList][
                                    ifelse vir >= 6 AND vir < 7 [set surTime item 2 fixedExtTimesList][
                                      ifelse vir >= 7 AND vir < 8 [set surTime item 2 fixedExtTimesList][
                                        ifelse vir >= 8 AND vir < 9 [set surTime item 1 fixedExtTimesList][
                                          ifelse vir >= 9 AND vir < 10 [set surTime item 1 fixedExtTimesList][
                                            ifelse vir >= 10 [set surTime item 0 fixedExtTimesList][
                                              if vir > 10 [set surTime item 0 fixedExtTimesList]
  ]]]]]]]]]]]]]]]]]]]]]]

  report surTime
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;; END                                                                                                                                                                                                                                                       ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@#$#@#$#@
GRAPHICS-WINDOW
916
11
1291
755
-1
-1
14.7
1
10
1
1
1
0
0
0
1
0
24
0
49
1
1
1
ticks
30.0

PLOT
8
154
353
354
Population during outbreak
Time [weeks]
Number of Individuals
0.0
10.0
0.0
100.0
true
true
"" ""
PENS
"Population" 1.0 0 -16777216 true "" "if  any? turtles with [hasActiveStrain = true] [plot count turtles]"
"mean Virulence" 1.0 0 -5298144 true "" "if  any? turtles with [hasActiveStrain = true] [plot (mean[strain] of turtles with [hasActiveStrain = true AND strain != \"none\"] * 3000)]"

PLOT
8
357
353
554
Epidemiological Status (All)
Time [weeks]
Number of Individuals
0.0
10.0
0.0
100.0
true
true
"" "if not any? turtles [ stop ]"
PENS
"Immune" 1.0 0 -14439633 true "" "plot count turtles with [EpiStat = \"esImm\"]"
"Immune Mat" 1.0 0 -11085214 true "" "plot count turtles with [EpiStat = \"esImmMat\"]"
"Lethally Inf" 1.0 0 -8053223 true "" "plot count turtles with [EpiStat = \"esLeth\"]"
"Transient Inf" 1.0 0 -817084 true "" "plot count turtles with [EpiStat = \"esTrans\"]"

PLOT
688
11
911
190
Control
NIL
NIL
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

TEXTBOX
692
26
905
185
i
11
139.9
0

BUTTON
695
29
750
62
Setup
Setup\n;View\nask turtles [ ht ]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
752
29
826
62
Run Week
Go\nView
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
828
29
902
62
Run Year
Go\nView\nif ticks mod 52 = 0 [stop]
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
811
64
902
97
Run Until End
if DONE = 0 [ Go View ]
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
695
64
808
97
Run Until Infection
while [ ticks + 1 < (WeekRelease) ] [ Go ]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
1650
17
1873
153
Landscape
NIL
NIL
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

TEXTBOX
1657
34
1868
148
i
11
139.9
0

PLOT
1650
157
1873
338
Host Population
NIL
NIL
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

TEXTBOX
1654
174
1868
334
i
11
139.9
0

SLIDER
1657
178
1759
211
Herd%
Herd%
0
100
100.0
1
1
NIL
HORIZONTAL

SLIDER
1762
178
1865
211
ReleaseFactor
ReleaseFactor
0
10
4.5
0.5
1
NIL
HORIZONTAL

SLIDER
1657
214
1759
247
Longevity
Longevity
0
50
11.0
1
1
yrs
HORIZONTAL

SLIDER
1762
214
1865
247
AgeBlur
AgeBlur
0
26
6.0
1
1
weeks
HORIZONTAL

PLOT
1650
340
1873
564
Disease Dynamics
NIL
NIL
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

TEXTBOX
1656
358
1869
559
i
11
139.9
0

INPUTBOX
1723
496
1799
556
BetaBetween
0.00208
1
0
Number

SLIDER
1657
360
1758
393
YearRelease
YearRelease
1
10
2.0
1
1
NIL
HORIZONTAL

SLIDER
1657
397
1758
430
PreNatInf
PreNatInf
0
1
0.5
0.05
1
NIL
HORIZONTAL

SLIDER
1762
360
1865
393
FertRedInf
FertRedInf
0.025
1
0.5
0.025
1
NIL
HORIZONTAL

INPUTBOX
1762
433
1812
493
T_trans
1.0
1
0
Number

INPUTBOX
1815
433
1865
493
T_anti
12.0
1
0
Number

INPUTBOX
695
100
784
160
SimulatedYears
100.0
1
0
Number

SLIDER
1762
396
1865
429
CaseFatality
CaseFatality
0
1
0.5
0.05
1
NIL
HORIZONTAL

MONITOR
916
753
966
798
Year
ceiling (ticks / 52)
0
1
11

PLOT
1876
357
2202
554
Age Classes of Infected Ind.
Age Class [years]
Numberof Infected Ind.
0.0
11.0
0.0
100.0
true
true
"" "if not any? turtles with [EpiStat = \"esLeth\" or EpiStat = \"esTrans\"]\n[\n  clear-plot\n  stop\n]"
PENS
"pen-0" 1.0 1 -14737633 true "" ""
"pen-1" 1.0 1 -5298144 true "" ""
"pen-2" 1.0 1 -817084 true "" ""
"pen-3" 1.0 0 -7500403 true "" ""
"oen-4" 1.0 0 -2674135 true "" ""
"Total Inf" 1.0 1 -14737633 true "" "histogram InfStructList"
"Lethally Inf" 1.0 1 -5298144 true "" "histogram InfLethStructList"
"Transient Inf" 1.0 1 -817084 true "" "histogram InfTransStructList"

SLIDER
1762
250
1865
283
DispDist
DispDist
0
25
3.0
1
1
NIL
HORIZONTAL

SLIDER
1657
250
1759
283
ProbFem
ProbFem
0
1
0.5
0.01
1
NIL
HORIZONTAL

TEXTBOX
1833
376
1983
394
NIL
11
0.0
1

MONITOR
2121
410
2195
451
New Transient
NewTrans
17
1
10

MONITOR
577
10
681
55
Run Time [min]
Time
17
1
11

PLOT
1876
557
2202
756
Survival Time of Lethally Infected Ind.
Duration Infectious Period [weeks]
Number of Infected
0.0
50.0
0.0
10.0
true
false
"" "if not any? turtles with [EpiStat = \"esLeth\"]\n[\n  clear-plot\n  stop\n]"
PENS
"default" 1.0 1 -14737633 true "" "histogram CurrInfPeriodList"

MONITOR
2141
578
2192
619
Maximum
max LethPeriodList
17
1
10

PLOT
11
560
356
759
Epidemiological Status (Infected Ind.)
Time [weeks]
Number of Infected Ind.
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Lethally Inf" 1.0 0 -8053223 true "" "plot count turtles with [EpiStat = \"esLeth\"]"
"Transient Inf" 1.0 0 -817084 true "" "plot count turtles with [EpiStat = \"esTrans\"]"

MONITOR
274
692
351
733
Transient Inf
count turtles with [EpiStat = \"esTrans\"]
17
1
10

MONITOR
274
652
351
693
Lethally Inf
count turtles with [EpiStat = \"esLeth\"]
17
1
10

MONITOR
274
613
351
654
Immunes
count turtles with [EpiStat = \"esImm\" or EpiStat = \"esImmMat\"]
17
1
10

MONITOR
274
574
351
615
Susceptibles
count turtles with [EpiStat = \"esSusc\"]
17
1
10

PLOT
1876
157
2202
357
Age Classes
Age Class [years]
Number of Individuals
0.0
11.0
0.0
10.0
true
false
"" "if not any? turtles\n[\n  clear-plot\n  stop\n]"
PENS
"default" 1.0 1 -16777216 true "" "histogram PopStructList"

MONITOR
2077
255
2138
296
Males
count turtles with [is_female = 0]
17
1
10

MONITOR
2077
215
2138
256
Females
count turtles with [is_female = 1]
17
1
10

MONITOR
2017
175
2078
216
Population
count turtles
17
1
10

MONITOR
2134
255
2191
296
Adults
count turtles with [AgeGroup = \"Adult\"]
17
1
10

MONITOR
2134
215
2191
256
Subadults
count turtles with [AgeGroup = \"Subadult\"]
17
1
10

MONITOR
2121
371
2195
412
New Lethal
NewLeth
17
1
10

MONITOR
965
753
1015
798
Week
week
17
1
11

OUTPUT
2202
154
2342
238
11

MONITOR
577
58
681
103
Ticks per second
ticks / (Time * 60)
1
1
11

CHOOSER
1660
37
1771
82
Landscape
Landscape
"Generator-Homogeneous" "Generator-Random" "File-Import-Path" "File-Import-User"
2

SLIDER
1774
85
1866
118
MeanQuality
MeanQuality
1
15
4.5
0.5
1
NIL
HORIZONTAL

INPUTBOX
1660
85
1772
145
Filename
patch_m_1
1
0
String

BUTTON
787
148
902
181
Default Values
;; CONTROL\nset SimulatedYears 50\nset seed-setup \"ON\"\nset seed 1\n\n;; LANDSCAPE\nset Landscape \"File-Import-Path\"\nset Filename \"patch_m_1\"\nset MeanQuality 4.5\nset Style \"continuous\"\nset survival_weeks 10\nset temporal_lag 0\nset dynmode TRUE\n\n;; HOST POPULATION\nset Herd% 100\nset ReleaseFactor 4.5\nset Longevity 11\nset AgeBlur 6\nset ProbFem 0.5\nset DispDist 3\nset uniform_breeding FALSE\n\n;; DISEASE DYNAMICS\nset YearRelease 2\nset FertRedInf 0.625\nset PreNatInf 0.5\nset CaseFatality 0.5\nset BetaWithin 0.0208\nset BetaBetween 0.00208\nset BetaMove 0.0224\nset mue 6\nset T_trans 1\nset T_anti 12\n\n;; MOVEMENT\nset q 0.5\nset Roaming \"CRW\"
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
2077
175
2138
216
Roamers
count turtles with [DemStat = \"dsDispM_adult\"]
17
1
10

SLIDER
1763
286
1866
319
q
q
0
1.0
0.5
0.05
1
NIL
HORIZONTAL

CHOOSER
1658
286
1760
331
Roaming
Roaming
"OFF" "RW" "CRW" "HD" "DD" "FD" "HDD" "HD-CRW"
0

CHOOSER
787
100
902
145
seed-setup
seed-setup
"OFF" "ON" "BS"
0

SLIDER
1657
433
1758
466
mue
mue
1
15
8.0
1
1
NIL
HORIZONTAL

MONITOR
2134
175
2191
216
Piglets
count turtles with [AgeGroup = \"Piglet\"]
17
1
10

CHOOSER
1773
37
1866
82
Style
Style
"continuous" "discrete"
0

INPUTBOX
734
120
784
180
seed
55.0
1
0
Number

SLIDER
1657
465
1760
498
release_week
release_week
0
52
20.0
1
1
NIL
HORIZONTAL

MONITOR
352
152
466
197
Mean Habitat Quality
mean [counter_d] of patches
17
1
11

SWITCH
1301
37
1439
70
dynmode
dynmode
0
1
-1000

SWITCH
1302
185
1441
218
uniform_breeding
uniform_breeding
1
1
-1000

PLOT
7
10
215
151
Habitat quality
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean[counter_d] of patches"

SLIDER
2023
50
2195
83
dyntimer_down
dyntimer_down
0
5
0.0
1
1
weeks
HORIZONTAL

SLIDER
1301
117
1440
150
survival_weeks
survival_weeks
0
25
10.0
1
1
weeks
HORIZONTAL

CHOOSER
1301
71
1439
116
noise_type
noise_type
"Red Noise" "White Noise"
0

SWITCH
2023
17
2194
50
experimental_rednoise
experimental_rednoise
1
1
-1000

MONITOR
2017
216
2074
261
Breeder
count turtles with [is_female = 1 AND is_mother = 0 AND Age > 34]
17
1
11

SLIDER
1301
151
1440
184
temporal_lag
temporal_lag
0
100
0.0
25
1
%
HORIZONTAL

SLIDER
2023
82
2195
115
dyntimer_up
dyntimer_up
0
5
0.0
1
1
weeks
HORIZONTAL

TEXTBOX
1315
15
1465
33
Landscape dynamic
12
0.0
0

TEXTBOX
1451
37
1601
65
ON-  dynamic landscape\nOFF- static landscape
11
0.0
1

TEXTBOX
1451
122
1654
150
maximum survival time of individuals without access to resources
11
0.0
1

TEXTBOX
1451
156
1643
198
temporal shift between resource peak and seasonal reproduction peak
11
0.0
1

TEXTBOX
1451
187
1625
229
switch from seasonal reproduction peak to year-round reproduction
11
0.0
1

TEXTBOX
2203
29
2353
99
Debugging\ndefaults: \nexperimental red noise - OFF\ndyntimer down - 0\ndyntimer up - 0
11
0.0
1

PLOT
410
190
687
370
Mean virulence over all individuals
Time [weeks]
Mean virulence
0.0
10.0
-10.0
10.0
true
false
"" ""
PENS
"pen_20" 1.0 0 -5825686 true "" "plot mean[strain] of turtles with [hasActiveStrain = true AND strain != \"none\"]"

TEXTBOX
1305
238
1455
256
Virulence evolution\n
14
0.0
1

CHOOSER
1302
385
1494
430
strainTransmissionType
strainTransmissionType
"percentile" "highest infection pressure" "random" "exp"
3

SLIDER
1304
639
1476
672
BetaMove
BetaMove
0
0.5
0.0541
0.0001
1
NIL
HORIZONTAL

SLIDER
1302
436
1474
469
BetaWithin
BetaWithin
0
0.3
0.0208
0.0001
1
NIL
HORIZONTAL

SLIDER
1650
566
1822
599
BetaWithinHV
BetaWithinHV
0
0.1
0.0439
0.0001
1
NIL
HORIZONTAL

SLIDER
1650
600
1822
633
BetaWithinLV
BetaWithinLV
0
0.1
0.0216
0.0001
1
NIL
HORIZONTAL

SLIDER
1651
633
1823
666
BetaMoveHV
BetaMoveHV
0
0.1
0.0363
0.0001
1
NIL
HORIZONTAL

SLIDER
1651
666
1823
699
BetaMoveLV
BetaMoveLV
0
0.1
0.0236
0.0001
1
NIL
HORIZONTAL

MONITOR
2415
10
2533
55
strain transfer fails
count turtles with [strainVirToGet = 4]
17
1
11

INPUTBOX
1302
325
1389
385
mutationRate
0.01
1
0
Number

INPUTBOX
1388
325
1543
385
virulenceEffectOnMortality
1.0
1
0
Number

SLIDER
1302
293
1474
326
evolutionFactor
evolutionFactor
0
5
3.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
1480
295
1630
323
SD of the normal \ndistribution
11
0.0
1

SWITCH
1302
261
1408
294
evolution
evolution
0
1
-1000

PLOT
687
546
911
723
Low strains
Time [weeks]
Number of individuals
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"-2 - -4" 1.0 0 -15302303 true "" "plot count turtles with [hasActiveStrain = true AND strain < -2 AND strain > -4 AND strain != \"none\"]"
"-4 - -6" 1.0 0 -12345184 true "" "plot count turtles with [hasActiveStrain = true AND strain < -4 AND  strain > -6 AND strain != \"none\"]"
"-6 - -8" 1.0 0 -11033397 true "" "plot count turtles with [hasActiveStrain = true AND strain < -6 AND strain > -8 AND strain != \"none\"]"
"< -8" 1.0 0 -14454117 true "" "plot count turtles with [hasActiveStrain = true AND strain < -8 AND strain != \"none\"]"

PLOT
687
370
911
545
Medium strains
Time [weeks]
Number of individuals
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"0 - 2" 1.0 0 -13840069 true "" "plot count turtles with [hasActiveStrain = true AND strain > 0 AND strain < 2 AND strain != 999]"
"0 - -2" 1.0 0 -11221820 true "" "plot count turtles with [hasActiveStrain = true AND strain < 0 AND strain > -2 AND strain != \"none\"]"

MONITOR
852
662
909
707
N_LV
count turtles with [strain < 0 AND HasactiveStrain = true AND strain != \"none\"]
17
1
11

MONITOR
854
475
911
520
N_MV
count turtles with [strain < 2 AND strain > -2 AND HasactiveStrain = true AND strain != 999\n]
17
1
11

TEXTBOX
1480
436
1630
464
mean infectiousness (neiborhood infection)
11
0.0
1

TEXTBOX
1483
642
1633
670
mean infectiousness (movement)
11
0.0
1

SWITCH
1450
72
1593
105
dynmode_visuals
dynmode_visuals
1
1
-1000

SWITCH
1303
508
1454
541
invasion_analysis
invasion_analysis
1
1
-1000

SWITCH
1304
605
1527
638
movement_virulence_tradeoff
movement_virulence_tradeoff
1
1
-1000

SLIDER
1651
700
1823
733
transmission_factor
transmission_factor
0
1
0.0
0.001
1
NIL
HORIZONTAL

INPUTBOX
1651
732
1729
792
lowVirOffset
0.0
1
0
Number

SWITCH
1302
469
1453
502
fixedSurvivalTimes
fixedSurvivalTimes
0
1
-1000

INPUTBOX
1303
540
1412
600
invasion_timeDelay
5.0
1
0
Number

SWITCH
1729
734
1869
767
LinearTradeOff
LinearTradeOff
1
1
-1000

CHOOSER
1460
482
1598
527
trade_off_type
trade_off_type
"Lin" "Lin2" "Non lin A" "Non lin B" "Non lin C" "Non lin D" "Non lin E" "Non lin F" "Non lin G"
2

PLOT
410
371
687
547
Strains of infected individuals
Strain
Number of Infected
-10.0
10.0
0.0
10.0
true
false
"" ""
PENS
"pen-0" 1.0 1 -7500403 true "" "histogram strainList"

MONITOR
629
206
679
251
Mean
mean[strain] of turtles with [hasActiveStrain = true AND strain != \"none\"]
2
1
11

PLOT
410
546
688
723
Max strain
Time [weeks]
N
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot  count turtles with [hasActiveStrain = true AND strain >= 10 AND strain != 999] "

PLOT
686
193
913
370
High strains
Time [weeks]
Individuals
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"2 - 4" 1.0 0 -612749 true "" "plot count turtles with [hasActiveStrain = true AND strain > 2 AND strain < 4 AND strain != 999]"
"4 - 6" 1.0 0 -10146808 true "" "plot count turtles with [hasActiveStrain = true AND strain > 4 AND strain < 6 AND strain != 999]"
"6 - 8" 1.0 0 -1604481 true "" "plot  count turtles with [hasActiveStrain = true AND strain >= 6 AND strain < 8 AND strain != 999] "
"8 - 10" 1.0 0 -5298144 true "" "plot  count turtles with [hasActiveStrain = true AND strain >= 8 AND strain < 11 AND strain != 999] "

MONITOR
855
307
912
352
N_HV
count turtles with [strain >= 2 AND HasactiveStrain = true AND strain != 999]
17
1
11

MONITOR
624
562
681
607
N
count turtles with [hasActiveStrain = true AND strain >= 10 AND strain != 999]
17
1
11

BUTTON
1292
720
1358
753
NIL
debug
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
215
10
547
151
Population
Time [weeks]
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count turtles"
"pen-1" 1.0 0 -7500403 true "" "plot mean[counter_d * 3000] of patches"

SWITCH
1874
52
2021
85
landscapeGscale
landscapeGscale
0
1
-1000

SWITCH
1874
18
2010
51
transListSelect
transListSelect
0
1
-1000

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Manuscript" repetitions="150" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="624"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>count turtles</metric>
    <metric>count turtles with [EpiStat = "esSusc"]</metric>
    <metric>count turtles with [EpiStat = "esTrans"]</metric>
    <metric>count turtles with [EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &lt;= 34]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &gt; 104]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esLeth"]</metric>
    <metric>NewTrans</metric>
    <metric>NewLeth</metric>
    <metric>count turtles with [is_female = 0]</metric>
    <metric>count turtles with [is_female = 1]</metric>
    <metric>MaxLethPeriod</metric>
    <metric>MeanLethPeriod</metric>
    <metric>count patches with [is_infected = 1]</metric>
    <metric>count patches with [is_infectious = 1]</metric>
    <metric>outbreak_group</metric>
    <metric>outbreak_roaming</metric>
    <metric>outbreak_dens5</metric>
    <metric>outbreak_dens12</metric>
    <metric>F_infected</metric>
    <metric>F_infectious</metric>
    <metric>InfDist</metric>
    <metric>MaxDist</metric>
    <metric>MeanMoveDist</metric>
    <metric>WeekRelease</metric>
    <metric>WeekLast</metric>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;hom&quot;"/>
      <value value="&quot;rand&quot;"/>
      <value value="&quot;patch_s&quot;"/>
      <value value="&quot;patch_m&quot;"/>
      <value value="&quot;patch_l&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="CaseFatality" first="0" step="0.25" last="1"/>
    <enumeratedValueSet variable="mue">
      <value value="3"/>
      <value value="6"/>
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;RW&quot;"/>
      <value value="&quot;CRW&quot;"/>
      <value value="&quot;HD&quot;"/>
      <value value="&quot;FD&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;BS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import-Path&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanQuality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Style">
      <value value="&quot;continuous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Herd%">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReleaseFactor">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AgeBlur">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProbFem">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DispDist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="YearRelease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PreNatInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaBetween">
      <value value="0.00208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaMove">
      <value value="0.0224"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_test" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="624"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>count turtles</metric>
    <metric>count turtles with [EpiStat = "esSusc"]</metric>
    <metric>count turtles with [EpiStat = "esTrans"]</metric>
    <metric>count turtles with [EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &lt;= 34]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &gt; 104]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esLeth"]</metric>
    <metric>NewTrans</metric>
    <metric>NewLeth</metric>
    <metric>count turtles with [is_female = 0]</metric>
    <metric>count turtles with [is_female = 1]</metric>
    <metric>MaxLethPeriod</metric>
    <metric>MeanLethPeriod</metric>
    <metric>count patches with [is_infected = 1]</metric>
    <metric>count patches with [is_infectious = 1]</metric>
    <metric>outbreak_group</metric>
    <metric>outbreak_roaming</metric>
    <metric>outbreak_dens5</metric>
    <metric>outbreak_dens12</metric>
    <metric>F_infected</metric>
    <metric>F_infectious</metric>
    <metric>InfDist</metric>
    <metric>MaxDist</metric>
    <metric>MeanMoveDist</metric>
    <metric>WeekRelease</metric>
    <metric>WeekLast</metric>
    <enumeratedValueSet variable="mue">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="12"/>
    </enumeratedValueSet>
    <steppedValueSet variable="release_week" first="5" step="2" last="50"/>
  </experiment>
  <experiment name="dyn_nodyn_run" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2601"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>count turtles</metric>
    <metric>count turtles with [EpiStat = "esSusc"]</metric>
    <metric>count turtles with [EpiStat = "esTrans"]</metric>
    <metric>count turtles with [EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &lt;= 34]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &gt; 104]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esLeth"]</metric>
    <metric>NewTrans</metric>
    <metric>NewLeth</metric>
    <metric>count turtles with [is_female = 0]</metric>
    <metric>count turtles with [is_female = 1]</metric>
    <metric>MaxLethPeriod</metric>
    <metric>MeanLethPeriod</metric>
    <metric>count patches with [is_infected = 1]</metric>
    <metric>count patches with [is_infectious = 1]</metric>
    <metric>outbreak_group</metric>
    <metric>outbreak_roaming</metric>
    <metric>outbreak_dens5</metric>
    <metric>outbreak_dens12</metric>
    <metric>F_infected</metric>
    <metric>F_infectious</metric>
    <metric>InfDist</metric>
    <metric>MaxDist</metric>
    <metric>MeanMoveDist</metric>
    <metric>WeekRelease</metric>
    <metric>WeekLast</metric>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;CRW&quot;"/>
      <value value="&quot;HD-CRW&quot;"/>
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="noise_type">
      <value value="&quot;Red Noise&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="6"/>
      <value value="8"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="survival_weeks">
      <value value="5"/>
      <value value="10"/>
      <value value="15"/>
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="dyn_nodyn_ru_uniform" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="5200"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>count turtles</metric>
    <metric>count turtles with [EpiStat = "esSusc"]</metric>
    <metric>count turtles with [EpiStat = "esTrans"]</metric>
    <metric>count turtles with [EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &lt;= 34]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &gt; 104]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esLeth"]</metric>
    <metric>NewTrans</metric>
    <metric>NewLeth</metric>
    <metric>count turtles with [is_female = 0]</metric>
    <metric>count turtles with [is_female = 1]</metric>
    <metric>MaxLethPeriod</metric>
    <metric>MeanLethPeriod</metric>
    <metric>count patches with [is_infected = 1]</metric>
    <metric>count patches with [is_infectious = 1]</metric>
    <metric>outbreak_group</metric>
    <metric>outbreak_roaming</metric>
    <metric>outbreak_dens5</metric>
    <metric>outbreak_dens12</metric>
    <metric>F_infected</metric>
    <metric>F_infectious</metric>
    <metric>InfDist</metric>
    <metric>MaxDist</metric>
    <metric>MeanMoveDist</metric>
    <metric>WeekRelease</metric>
    <metric>WeekLast</metric>
    <metric>mean[counter_d] of patches</metric>
    <metric>ceiling (ticks / 52)</metric>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;CRW&quot;"/>
      <value value="&quot;HD-CRW&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="white_noise">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_week">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="6"/>
      <value value="8"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2600"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>o_list</metric>
    <enumeratedValueSet variable="noise_type">
      <value value="&quot;Red Noise&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Style">
      <value value="&quot;continuous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AgeBlur">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="survival_weeks">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CaseFatality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experimental_rednoise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.575"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;Generator-Random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaBetween">
      <value value="0.00208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaMove">
      <value value="0.0224"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_down">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Herd%">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReleaseFactor">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;CRW&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_week">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_up">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="YearRelease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PreNatInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanQuality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProbFem">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;patch_m_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DispDist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="hdcrw_cluster" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2600"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>o_list</metric>
    <enumeratedValueSet variable="noise_type">
      <value value="&quot;Red Noise&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Style">
      <value value="&quot;continuous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AgeBlur">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="survival_weeks">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CaseFatality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experimental_rednoise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.575"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import-Path&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaBetween">
      <value value="0.00208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaMove">
      <value value="0.0224"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_down">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Herd%">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReleaseFactor">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;HD-CRW&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_week">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_up">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="YearRelease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PreNatInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanQuality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProbFem">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;patch_m_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DispDist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="full_mismatch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="climate_change_factor">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="climate_slider_run" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2600"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>count turtles</metric>
    <metric>count turtles with [EpiStat = "esSusc"]</metric>
    <metric>count turtles with [EpiStat = "esTrans"]</metric>
    <metric>count turtles with [EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &lt;= 34]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &gt; 104]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esSusc"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [DemStat = "dsDispM_adult" AND EpiStat = "esLeth"]</metric>
    <metric>NewTrans</metric>
    <metric>NewLeth</metric>
    <metric>count turtles with [is_female = 0]</metric>
    <metric>count turtles with [is_female = 1]</metric>
    <metric>MaxLethPeriod</metric>
    <metric>MeanLethPeriod</metric>
    <metric>count patches with [is_infected = 1]</metric>
    <metric>count patches with [is_infectious = 1]</metric>
    <metric>outbreak_group</metric>
    <metric>outbreak_roaming</metric>
    <metric>outbreak_dens5</metric>
    <metric>outbreak_dens12</metric>
    <metric>F_infected</metric>
    <metric>F_infectious</metric>
    <metric>InfDist</metric>
    <metric>MaxDist</metric>
    <metric>MeanMoveDist</metric>
    <metric>WeekRelease</metric>
    <metric>WeekLast</metric>
    <metric>mean[counter_d] of patches</metric>
    <metric>ceiling (ticks / 52)</metric>
    <enumeratedValueSet variable="noise_type">
      <value value="&quot;Red Noise&quot;"/>
      <value value="&quot;White Noise&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Style">
      <value value="&quot;continuous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AgeBlur">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CaseFatality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="survival_weeks">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experimental_rednoise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="full_mismatch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.575"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import-Path&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaBetween">
      <value value="0.00208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaMove">
      <value value="0.0224"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_down">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="climate_change_factor">
      <value value="0"/>
      <value value="25"/>
      <value value="50"/>
      <value value="75"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReleaseFactor">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Herd%">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;HD-CRW&quot;"/>
      <value value="&quot;CRW&quot;"/>
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_week">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_up">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="YearRelease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PreNatInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanQuality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProbFem">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DispDist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;patch_m_1&quot;"/>
      <value value="&quot;patch_s_1&quot;"/>
      <value value="&quot;patch_l_1&quot;"/>
      <value value="&quot;patch_m_2&quot;"/>
      <value value="&quot;patch_s_2&quot;"/>
      <value value="&quot;patch_l_2&quot;"/>
      <value value="&quot;rand_1&quot;"/>
      <value value="&quot;rand_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="hdcrw_cluster_white_noise" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2600"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>o_list</metric>
    <enumeratedValueSet variable="noise_type">
      <value value="&quot;White Noise&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Style">
      <value value="&quot;continuous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AgeBlur">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="survival_weeks">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CaseFatality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experimental_rednoise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.575"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import-Path&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaBetween">
      <value value="0.00208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaMove">
      <value value="0.0224"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_down">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Herd%">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReleaseFactor">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;HD-CRW&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_week">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_up">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="YearRelease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PreNatInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanQuality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProbFem">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;patch_m_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DispDist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="full_mismatch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="climate_change_factor">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="scluster" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2600"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>o_list</metric>
    <enumeratedValueSet variable="noise_type">
      <value value="&quot;Red Noise&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Style">
      <value value="&quot;continuous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AgeBlur">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CaseFatality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="survival_weeks">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experimental_rednoise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="full_mismatch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.575"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import-Path&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaBetween">
      <value value="0.00208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaMove">
      <value value="0.0224"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_down">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Herd%">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="climate_change_factor">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReleaseFactor">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;HD-CRW&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_week">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_up">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="YearRelease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PreNatInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanQuality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProbFem">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;patch_m_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DispDist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2600"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>o_list</metric>
    <enumeratedValueSet variable="noise_type">
      <value value="&quot;Red Noise&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Style">
      <value value="&quot;continuous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;ON&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AgeBlur">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CaseFatality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="survival_weeks">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experimental_rednoise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="full_mismatch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.575"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;Generator-Random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaBetween">
      <value value="0.00208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaMove">
      <value value="0.0224"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_down">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Herd%">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="climate_change_factor">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReleaseFactor">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;HD-CRW&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_week">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_up">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="YearRelease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PreNatInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanQuality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProbFem">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;patch_l_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DispDist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_crw_rework" repetitions="25" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2601"/>
    <exitCondition>Done = 1</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="noise_type">
      <value value="&quot;Red Noise&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Style">
      <value value="&quot;continuous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AgeBlur">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CaseFatality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="survival_weeks">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experimental_rednoise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="full_mismatch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.575"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import-Path&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaBetween">
      <value value="0.00208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaMove">
      <value value="0.0224"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_down">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Herd%">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="climate_change_factor">
      <value value="0"/>
      <value value="25"/>
      <value value="50"/>
      <value value="75"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReleaseFactor">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="6"/>
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;CRW&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_week">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_up">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="YearRelease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PreNatInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanQuality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProbFem">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;patch_s_1&quot;"/>
      <value value="&quot;patch_l_1&quot;"/>
      <value value="&quot;patch_m_1&quot;"/>
      <value value="&quot;rand_2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DispDist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_occ" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2601"/>
    <exitCondition>Done = 1</exitCondition>
    <metric>o_list</metric>
    <enumeratedValueSet variable="noise_type">
      <value value="&quot;Red Noise&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Style">
      <value value="&quot;continuous&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SimulatedYears">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AgeBlur">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="survival_weeks">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CaseFatality">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="experimental_rednoise">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="full_mismatch">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.575"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import-Path&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaMove">
      <value value="0.0224"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaBetween">
      <value value="0.00208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_down">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="climate_change_factor">
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Herd%">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ReleaseFactor">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Roaming">
      <value value="&quot;HD-CRW&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_week">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dyntimer_up">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="YearRelease">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PreNatInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanQuality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ProbFem">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;patch_m_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="DispDist">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Beta_test" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1300"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>WeekLast</metric>
    <steppedValueSet variable="BetaWithin" first="0.01" step="0.001" last="0.05"/>
    <enumeratedValueSet variable="transmission_factor">
      <value value="0.01"/>
      <value value="0.02"/>
      <value value="0.03"/>
      <value value="0.04"/>
      <value value="0.05"/>
      <value value="0.06"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Beta_test_2" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1300"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>WeekLast</metric>
    <steppedValueSet variable="BetaWithin" first="0.03" step="0.001" last="0.09"/>
  </experiment>
  <experiment name="tr_factor_test" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1300"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>WeekLast</metric>
    <steppedValueSet variable="transmission_factor" first="0.01" step="0.001" last="0.05"/>
  </experiment>
  <experiment name="Viru_evo_testRun_1" repetitions="20" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="5201"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>count turtles</metric>
    <metric>meanStrainOutput</metric>
    <metric>count turtles with [EpiStat = "esSusc"]</metric>
    <metric>count turtles with [EpiStat = "esLeth"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 0 AND strain &lt; 2 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 2 AND strain &lt; 4 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 4 AND strain &lt; 5 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 5 AND strain &lt; 6 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 6 AND strain &lt; 8 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 8 AND strain &lt; 10 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 10 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; 0 AND strain &gt; -2 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -2 AND strain &gt; -4 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -4 AND  strain &gt; -6 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -6 AND strain &gt; -8 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -8 AND strain != "none"]</metric>
    <metric>WeekLast</metric>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;hom_1&quot;"/>
      <value value="&quot;rand_1&quot;"/>
      <value value="&quot;rand_3&quot;"/>
      <value value="&quot;rand_6&quot;"/>
      <value value="&quot;patch_s_1&quot;"/>
      <value value="&quot;patch_s_3&quot;"/>
      <value value="&quot;patch_s_6&quot;"/>
      <value value="&quot;patch_m_1&quot;"/>
      <value value="&quot;patch_m_3&quot;"/>
      <value value="&quot;patch_m_6&quot;"/>
      <value value="&quot;patch_l_1&quot;"/>
      <value value="&quot;patch_l_2&quot;"/>
      <value value="&quot;patch_l_6&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="LinearTradeOff">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Viru_evo_testRun_multiTradeOff" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="5201"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>count turtles</metric>
    <metric>meanStrainOutput</metric>
    <metric>count turtles with [EpiStat = "esSusc"]</metric>
    <metric>count turtles with [EpiStat = "esLeth"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 0 AND strain &lt; 2 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 2 AND strain &lt; 4 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 4 AND strain &lt; 5 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 5 AND strain &lt; 6 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 6 AND strain &lt; 8 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 8 AND strain &lt; 10 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 10 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; 0 AND strain &gt; -2 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -2 AND strain &gt; -4 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -4 AND  strain &gt; -6 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -6 AND strain &gt; -8 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -8 AND strain != "none"]</metric>
    <metric>WeekLast</metric>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;hom_1&quot;"/>
      <value value="&quot;rand_1&quot;"/>
      <value value="&quot;patch_s_1&quot;"/>
      <value value="&quot;patch_m_1&quot;"/>
      <value value="&quot;patch_l_6&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="trade_off_type">
      <value value="&quot;Lin&quot;"/>
      <value value="&quot;Non lin A&quot;"/>
      <value value="&quot;Non lin B&quot;"/>
      <value value="&quot;Non lin C&quot;"/>
      <value value="&quot;Non lin D&quot;"/>
      <value value="&quot;Non lin E&quot;"/>
      <value value="&quot;Non lin F&quot;"/>
      <value value="&quot;Non lin G&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evolutionFactor">
      <value value="2.5"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import-Path&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="vir_evo_first_sims" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="5200"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>count turtles</metric>
    <metric>meanStrainOutput</metric>
    <metric>count turtles with [EpiStat = "esSusc"]</metric>
    <metric>count turtles with [EpiStat = "esLeth"]</metric>
    <metric>count turtles with [EpiStat = "esImm"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 0 AND strain &lt; 2 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 2 AND strain &lt; 4 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 4 AND strain &lt; 5 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 5 AND strain &lt; 6 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 6 AND strain &lt; 8 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 8 AND strain &lt; 10 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 10 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; 0 AND strain &gt; -2 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -2 AND strain &gt; -4 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -4 AND  strain &gt; -6 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -6 AND strain &gt; -8 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -8 AND strain != "none"]</metric>
    <metric>WeekLast</metric>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;hom_1&quot;"/>
      <value value="&quot;patch_m_1&quot;"/>
      <value value="&quot;patch_l_2&quot;"/>
      <value value="&quot;patch_s_1&quot;"/>
      <value value="&quot;rand_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="trade_off_type">
      <value value="&quot;Lin&quot;"/>
      <value value="&quot;Lin2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evolutionFactor">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed">
      <value value="55"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import-Path&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temporal_lag">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transListSelect">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="vir_evo_first_sims_shrt" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="5200"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>count turtles</metric>
    <metric>meanStrainOutput</metric>
    <metric>count turtles with [EpiStat = "esSusc"]</metric>
    <metric>count turtles with [EpiStat = "esLeth"]</metric>
    <metric>count turtles with [EpiStat = "esImm"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 0 AND strain &lt; 2 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 2 AND strain &lt; 4 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt; 4 AND strain &lt; 5 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 5 AND strain &lt; 6 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 6 AND strain &lt; 8 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 8 AND strain &lt; 10 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &gt;= 10 AND strain != 999]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; 0 AND strain &gt; -2 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -2 AND strain &gt; -4 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -4 AND  strain &gt; -6 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -6 AND strain &gt; -8 AND strain != "none"]</metric>
    <metric>count turtles with [hasActiveStrain = true AND strain &lt; -8 AND strain != "none"]</metric>
    <metric>WeekLast</metric>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;patch_m_1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="trade_off_type">
      <value value="&quot;Lin&quot;"/>
      <value value="&quot;Lin2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform_breeding">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evolutionFactor">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-setup">
      <value value="&quot;OFF&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed">
      <value value="55"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="FertRedInf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import-Path&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="BetaWithin">
      <value value="0.0208"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dynmode">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temporal_lag">
      <value value="0"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transListSelect">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
