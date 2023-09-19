# Master helper file of ODE functions ## MULTIPLE TURNOVER ##

# CTS is the used model
# MAX is monophasic

exp.params = function(p){
  # Exp transform k constants among parameters
  which.k = names(p) != "Fa"
  p[which.k] = exp(p[which.k])
  # print(p)
  return(p)
}

# ODE system
# Given state and parameters, gives you the d/dt of that time point
ODEsys.MAX = function(t, state, parameters) {
  parameters = exp.params(parameters)  # transform back to linear form

  kon = parameters["kon"]
  kcat = parameters["kcat"]
  kph2 = 0
  koffP = parameters["koffP"]
  Fa = parameters["Fa"]
  A = state["A"]; I = state["I"]; R = state["R"]
  AR = state["AR"]; IR = state["IR"]; P = state["P"]
  AP = state["AP"]; IP = state["IP"]

  dA  = - kon*R*A                                 + koffP*AP - kon*A*P
  dI  = 0
  dR  = - kon*R*A
  dAR = + kon*R*A - kcat*AR
  dIR = 0
  dP  =                                           + koffP*AP - kon*A*P
  dAP =           + kcat*AR                       - koffP*AP + kon*A*P
  dIP = 0
  list(c(dA, dI, dR, dAR, dIR, dP, dAP, dIP))
}

ODEsys.CTR = function(t, state, parameters) {
  parameters = exp.params(parameters)  # transform back to linear form

  kon = parameters["kon"]
  kcat = parameters["kcat"]
  kph2 = parameters["kph2"]
  koffP = parameters["koffP"]
  Fa = parameters["Fa"]
  A = state["A"]; I = state["I"]; R = state["R"]
  AR = state["AR"]; IR = state["IR"]; P = state["P"]
  AP = state["AP"]; IP = state["IP"]

  dA  = - kon*R*A                                 + koffP*AP            - kon*A*P
  dI  =                     - kon*R*I + kph2*IR              + koffP*IP           - kon*I*P
  dR  = - kon*R*A           - kon*R*I + kph2*IR
  dAR = + kon*R*A - kcat*AR
  dIR =                     + kon*R*I - kph2*IR
  dP  =                                           + koffP*AP + koffP*IP - kon*A*P - kon*I*P
  dAP =           + kcat*AR                       - koffP*AP            + kon*A*P
  dIP =                                                      - koffP*IP           + kon*I*P
  list(c(dA, dI, dR, dAR, dIR, dP, dAP, dIP))
}

ODEsys.CTS = function(t, state, parameters) {
  parameters = exp.params(parameters)  # transform back to linear form

  kon = parameters["kon"]
  kcat = parameters["kcat"]
  kph2 = parameters["kph2"]
  koffP = parameters["koffP"]
  Fa = parameters["Fa"]
  A = state["A"]; I = state["I"]; R = state["R"]
  AR = state["AR"]; IR = state["IR"]; P = state["P"]
  AP = state["AP"]; IP = state["IP"]

  dA  = - kon*R*A                                 + koffP*AP            - kon*A*P
  dI  =                     - kon*R*I                        + koffP*IP           - kon*I*P
  dR  = - kon*R*A           - kon*R*I
  dAR = + kon*R*A - kcat*AR
  dIR =                     + kon*R*I - kph2*IR
  dP  =                                           + koffP*AP + koffP*IP - kon*A*P - kon*I*P
  dAP =           + kcat*AR                       - koffP*AP            + kon*A*P
  dIP =                               + kph2*IR              - koffP*IP           + kon*I*P
  list(c(dA, dI, dR, dAR, dIR, dP, dAP, dIP))
}

ODEsys.UNL = function(t, state, parameters) {
  parameters = exp.params(parameters)  # transform back to linear form

  kon = parameters["kon"]
  kcat = parameters["kcat"]
  kph2 = parameters["kph2"]
  koffP = parameters["koffP"]
  Fa = parameters["Fa"]
  A = state["A"]; I = state["I"]; R = state["R"]
  AR = state["AR"]; IR = state["IR"]; P = state["P"]
  AP = state["AP"]; IP = state["IP"]

  kcomp = parameters["kcat"] * (1-Fa)/Fa

  dA  = - kon*R*A                                  + koffP*AP - kon*A*P
  dI  =                                + kph2*IR
  dR  = - kon*R*A                      + kph2*IR
  dAR = + kon*R*A - kcat*AR - kcomp*AR
  dIR =                     + kcomp*AR - kph2*IR
  dP  =                                            + koffP*AP - kon*A*P
  dAP =           + kcat*AR                        - koffP*AP + kon*A*P
  dIP = 0
  list(c(dA, dI, dR, dAR, dIR, dP, dAP, dIP))
}

ODEsys.CFI = function(t, state, parameters) {
  parameters = exp.params(parameters)  # transform back to linear form

  kon = parameters["kon"]
  kcat = parameters["kcat"]
  kph2 = parameters["kph2"]
  koffP = parameters["koffP"]
  Fa = parameters["Fa"]
  A = state["A"]; I = state["I"]; R = state["R"]
  AR = state["AR"]; IR = state["IR"]; P = state["P"]
  AP = state["AP"]; IP = state["IP"]

  kcomp = parameters["kcat"] * (1-Fa)/Fa

  dA  = - kon*R*A                      + kph2*IR   + koffP*AP - kon*A*P
  dI  = 0
  dR  = - kon*R*A                      + kph2*IR
  dAR = + kon*R*A - kcat*AR - kcomp*AR
  dIR =                     + kcomp*AR - kph2*IR
  dP  =                                            + koffP*AP - kon*A*P
  dAP =           + kcat*AR                        - koffP*AP + kon*A*P
  dIP = 0
  list(c(dA, dI, dR, dAR, dIR, dP, dAP, dIP))
}

ODEsys.CFR = function(t, state, parameters) {
  parameters = exp.params(parameters)  # transform back to linear form

  kon = parameters["kon"]
  kcat = parameters["kcat"]
  kph2 = parameters["kph2"]
  koffP = parameters["koffP"]
  Fa = parameters["Fa"]
  A = state["A"]; I = state["I"]; R = state["R"]
  AR = state["AR"]; IR = state["IR"]; P = state["P"]
  AP = state["AP"]; IP = state["IP"]

  kcomp = parameters["kcat"] * (1-Fa)/Fa

  dA  = - kon*R*A                                  + koffP*AP - kon*A*P
  dI  = 0
  dR  = - kon*R*A
  dAR = + kon*R*A - kcat*AR - kcomp*AR + kph2*IR
  dIR =                     + kcomp*AR - kph2*IR
  dP  =                                            + koffP*AP - kon*A*P
  dAP =           + kcat*AR                        - koffP*AP + kon*A*P
  dIP = 0
  list(c(dA, dI, dR, dAR, dIR, dP, dAP, dIP))
}

ODEsys.CFS = function(t, state, parameters) {
  parameters = exp.params(parameters)  # transform back to linear form

  kon = parameters["kon"]
  kcat = parameters["kcat"]
  kph2 = parameters["kph2"]
  koffP = parameters["koffP"]
  Fa = parameters["Fa"]
  A = state["A"]; I = state["I"]; R = state["R"]
  AR = state["AR"]; IR = state["IR"]; P = state["P"]
  AP = state["AP"]; IP = state["IP"]

  kcomp = parameters["kcat"] * (1-Fa)/Fa

  dA  = - kon*R*A                                  + koffP*AP - kon*A*P
  dI  = 0
  dR  = - kon*R*A
  dAR = + kon*R*A - kcat*AR - kcomp*AR
  dIR =                     + kcomp*AR - kph2*IR
  dP  =                                            + koffP*AP - kon*A*P
  dAP =           + kcat*AR            + kph2*IR   - koffP*AP + kon*A*P
  dIP = 0
  list(c(dA, dI, dR, dAR, dIR, dP, dAP, dIP))
}

ODEsys.BNR = function(t, state, parameters) {
  parameters = exp.params(parameters)  # transform back to linear form

  kon = parameters["kon"]
  kcat = parameters["kcat"]
  kph2 = parameters["kph2"]
  koffP = parameters["koffP"]
  Fa = parameters["Fa"]
  A = state["A"]; I = state["I"]; R = state["R"]
  AR = state["AR"]; IR = state["IR"]; P = state["P"]
  AP = state["AP"]; IP = state["IP"]

  kcomp = parameters["kon"] * (1-Fa)/Fa

  dA  = - kon*R*A           - kcomp*R*A + kph2*IR   + koffP*AP - kon*A*P
  dI  = 0
  dR  = - kon*R*A           - kcomp*R*A + kph2*IR
  dAR = + kon*R*A - kcat*AR
  dIR =                     + kcomp*R*A - kph2*IR
  dP  =                                             + koffP*AP - kon*A*P
  dAP =           + kcat*AR                         - koffP*AP + kon*A*P
  dIP = 0
  list(c(dA, dI, dR, dAR, dIR, dP, dAP, dIP))
}

ODEsys.BNS = function(t, state, parameters) {
  parameters = exp.params(parameters)  # transform back to linear form

  kon = parameters["kon"]
  kcat = parameters["kcat"]
  kph2 = parameters["kph2"]
  koffP = parameters["koffP"]
  Fa = parameters["Fa"]
  A = state["A"]; I = state["I"]; R = state["R"]
  AR = state["AR"]; IR = state["IR"]; P = state["P"]
  AP = state["AP"]; IP = state["IP"]

  kcomp = parameters["kon"] * (1-Fa)/Fa

  dA  = - kon*R*A           - kcomp*R*A             + koffP*AP - kon*A*P
  dI  = 0
  dR  = - kon*R*A           - kcomp*R*A
  dAR = + kon*R*A - kcat*AR
  dIR =                     + kcomp*R*A - kph2*IR
  dP  =                                             + koffP*AP - kon*A*P
  dAP =           + kcat*AR             + kph2*IR   - koffP*AP + kon*A*P
  dIP = 0
  list(c(dA, dI, dR, dAR, dIR, dP, dAP, dIP))
}


ODElib = function(modelmode){
  return(switch(modelmode,
                MAX = ODEsys.MAX,
                CTR = ODEsys.CTR,
                CTS = ODEsys.CTS,
                UNL = ODEsys.UNL,
                CFI = ODEsys.CFI,
                CFR = ODEsys.CFR,
                CFS = ODEsys.CFS,
                BNS = ODEsys.BNS,
                BNR = ODEsys.BNR
  ))
}

ODEstateinit = function(modelmode, A_T, R_T, parameters){
  # Generate initial states for ODE to run

  # Get Fa parameter
  Fa = parameters["Fa"]
  Fi = 1-Fa

  state = switch(modelmode,
                 #       A       I       R       AR IR P  AP IP
                 MAX = c(A_T,    0,      R_T*Fa, 0, 0, 0, 0, 0),
                 CTR = c(A_T*Fa, A_T*Fi, R_T,    0, 0, 0, 0, 0),
                 CTS = c(A_T*Fa, A_T*Fi, R_T,    0, 0, 0, 0, 0),
                 UNL = c(A_T,    0,      R_T,    0, 0, 0, 0, 0),
                 CFI = c(A_T,    0,      R_T,    0, 0, 0, 0, 0),
                 CFR = c(A_T,    0,      R_T,    0, 0, 0, 0, 0),
                 CFS = c(A_T,    0,      R_T,    0, 0, 0, 0, 0),
                 BNS = c(A_T,    0,      R_T,    0, 0, 0, 0, 0),
                 BNR = c(A_T,    0,      R_T,    0, 0, 0, 0, 0)
  )

  names(state) = c("A","I","R","AR","IR","P","AP","IP")
  return(state)
}

