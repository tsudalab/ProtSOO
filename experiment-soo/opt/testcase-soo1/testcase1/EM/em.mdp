integrator  = steep     ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0    ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size

nsteps      = 1000       ; Maximum number of (minimization) steps to perform
; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist     = 20        ; Frequency to update the neighbor list and long range forces
ns_type     = grid      ; Method to determine neighbor list (simple, grid)
rlist       = 1.0       ; Cut-off for making neighbor list (short range forces)
coulombtype = cut-off   ; Treatment of long range electrostatic interactions
rcoulomb    = 1.0       ; Short-range electrostatic cut-off
rvdw        = 1.0       ; Short-range Van der Waals cut-off
pbc         = xy        ; Periodic Boundary Conditions (yes/no)
ewald_geometry          = 3dc

; cutoff-scheme for gromacs version 5 of above
cutoff-scheme = group

; Output frequency and precision for xtc file
nstxout                 = 1
nstvout                 = 100
nstfout                 = 0
nstenergy               =  1 ; step
nstxtcout               =  1 ; step
;; The following parameters will be updated by ProtPOS program
;energygrps              = Protein PTS
;energygrp_excl          = PTS PTS
energygrps              = Protein PTS
energygrp_excl          = PTS PTS

; GENERALIZED BORN ELECTROSTATICS
implicit_solvent        = GBSA
; Algorithm for calculating Born radii
gb_algorithm            = OBC
; Frequency of calculating the Born radii inside rlist
nstgbradii              = 1
; Cutoff for Born radii calculation; the contribution from atoms
; between rlist and rgbradii is updated every nstlist steps
rgbradii                = 1
; Dielectric coefficient of the implicit solvent
gb_epsilon_solvent      = 78.3
; Salt concentration in M for Generalized Born models
gb_saltconc             = 0
; Scaling factors used in the OBC GB model. Default values are OBC(II)
gb_obc_alpha            = 1
gb_obc_beta             = 0.8
gb_obc_gamma            = 4.85
gb_dielectric_offset    = 0.009
sa_algorithm            = Ace-approximation
; Surface tension (kJ/mol/nm^2) for the SA (nonpolar surface) part of GBSA
; The value -1 will set default value for Still/HCT/OBC GB-models.
sa_surface_tension      = -1



; Non-equilibrium MD stuff
acc-grps                =
accelerate              =

;; The following parameters will be updated by ProtPOS program
;freezegrps              = PTS
freezegrps              = PTS
freezedim               = Y Y Y
cos-acceleration        = 0
deform
