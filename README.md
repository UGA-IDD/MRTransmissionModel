# MRTransmissionModel

This is an measles and rubella transmission model. It is an age-structured discrete time stochastic transmission model that incorporates both epidemiological and demographic transitions, building on framework introduced by [1,2]. The model structure described here was originally presented in [3] and was used to describe rubella dynamics in [4]. We structured the population into five epidemiological stages (maternally immune ‘M', susceptible ‘S', infected ‘I', recovered ‘R', and vaccinated ‘V', taken to indicate the effectively vaccinated). The key feature of the model is a matrix that at every time-step defines transitions from every possible epidemiological stage and age class combination to every other possible epidemiological stage and age class combination. 

[1] Klepac P, Caswell H. The stage-structured epidemic: linking disease and demography with a multi-state matrix approach model. Theor Ecol 2011;4:301–19. https://doi.org/10.1007/s12080-010-0079-8.

[2] Klepac P, Pomeroy LW, Bjornstad ON, Kuiken T, Osterhaus ADME, Rijks JM. Stage-structured transmission of phocine distemper virus in the Dutch 2002 outbreak. Proceedings of the Royal Society B-Biological Sciences 2009;276:2469–76. https://doi.org/10.1098/rspb.2009.0175.

[3] Metcalf CJE, Lessler J, Klepac P, Morice A, Grenfell BT, Bjornstad ON. Structured models of infectious disease: Inference with discrete data. Theor 
Popul Biol 2012;82:275–82. https://doi.org/DOI 10.1016/j.tpb.2011.12.001.

[4] Metcalf CJE, Lessler J, Klepac P, Cutts F, Grenfell BT. Impact of birth rate, seasonality and transmission rate on minimum levels of coverage needed for rubella vaccination. Epidemiol Infect 2012;140:2290–301. https://doi.org/Doi 10.1017/S0950268812000131.

Contributors: Amy Winter, Jess Metcalf, Justin Lessler, Petra Klepac, Lauren Hoskovec, Indrajit Ghosh, Matt Ferrari


