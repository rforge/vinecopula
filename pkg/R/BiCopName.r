BiCopName <- function(family, short = TRUE) {
    if (is.logical(short) == FALSE) 
        stop("'short' has to be a logical variable.")
    
    fam <- NA
    
    if (is.numeric(family)) {
        # Zahl zu Name
        if (short == TRUE) {
            # kurzer Name
            if (family == 0) 
                fam <- "I"
            if (family == 1) 
                fam <- "N"
            if (family == 2) 
                fam <- "t"
            if (family == 3) 
                fam <- "C"
            if (family == 4) 
                fam <- "G"
            if (family == 5) 
                fam <- "F"
            if (family == 6) 
                fam <- "J"
            if (family == 7) 
                fam <- "BB1"
            if (family == 8) 
                fam <- "BB6"
            if (family == 9) 
                fam <- "BB7"
            if (family == 10) 
                fam <- "BB8"
            if (family == 13) 
                fam <- "SC"
            if (family == 14) 
                fam <- "SG"
            if (family == 16) 
                fam <- "SJ"
            if (family == 17) 
                fam <- "SBB1"
            if (family == 18) 
                fam <- "SBB6"
            if (family == 19) 
                fam <- "SBB7"
            if (family == 20) 
                fam <- "SBB8"
            if (family == 23) 
                fam <- "C90"
            if (family == 24) 
                fam <- "G90"
            if (family == 26) 
                fam <- "J90"
            if (family == 27) 
                fam <- "BB1_90"
            if (family == 28) 
                fam <- "BB6_90"
            if (family == 29) 
                fam <- "BB7_90"
            if (family == 30) 
                fam <- "BB8_90"
            if (family == 33) 
                fam <- "C270"
            if (family == 34) 
                fam <- "G270"
            if (family == 36) 
                fam <- "J270"
            if (family == 37) 
                fam <- "BB1_270"
            if (family == 38) 
                fam <- "BB6_270"
            if (family == 39) 
                fam <- "BB7_270"
            if (family == 40) 
                fam <- "BB8_270"
            if (family == 100) 
                fam <- "NP"  #changed Mathias
            if (family == 41) 
                fam <- "1-par AS"
            if (family == 51) 
                fam <- "1-par AS180"
            if (family == 61) 
                fam <- "1-par AS90"
            if (family == 71) 
                fam <- "1-par AS270"
            if (family == 104) 
                fam <- "Tawn"
            if (family == 114) 
                fam <- "Tawn180"
            if (family == 124) 
                fam <- "Tawn90"
            if (family == 134) 
                fam <- "Tawn270"
            if (family == 204) 
                fam <- "Tawn2"
            if (family == 214) 
                fam <- "Tawn2_180"
            if (family == 224) 
                fam <- "Tawn2_90"
            if (family == 234) 
                fam <- "Tawn2_270"
        } else {
            # langer Name
            if (family == 0) 
                fam <- "Independence"
            if (family == 1) 
                fam <- "Gaussian"
            if (family == 2) 
                fam <- "t"
            if (family == 3) 
                fam <- "Clayton"
            if (family == 4) 
                fam <- "Gumbel"
            if (family == 5) 
                fam <- "Frank"
            if (family == 6) 
                fam <- "Joe"
            if (family == 7) 
                fam <- "Clayton-Gumbel"
            if (family == 8) 
                fam <- "Joe-Gumbel"
            if (family == 9) 
                fam <- "Joe-Clayton"
            if (family == 10) 
                fam <- "Joe-Frank"
            if (family == 13) 
                fam <- "Survival Clayton"
            if (family == 14) 
                fam <- "Survival Gumbel"
            if (family == 16) 
                fam <- "Survival Joe"
            if (family == 17) 
                fam <- "Survival Clayton-Gumbel"
            if (family == 18) 
                fam <- "Survival Joe-Gumbel"
            if (family == 19) 
                fam <- "Survival Joe-Clayton"
            if (family == 20) 
                fam <- "Survival Joe-Frank"
            if (family == 23) 
                fam <- "Rotated Clayton 90 degrees"
            if (family == 24) 
                fam <- "Rotated Gumbel 90 degrees"
            if (family == 26) 
                fam <- "Rotated Joe 90 degrees"
            if (family == 27) 
                fam <- "Rotated Clayton-Gumbel 90 degrees"
            if (family == 28) 
                fam <- "Rotated Joe-Gumbel 90 degrees"
            if (family == 29) 
                fam <- "Rotated Joe-Clayton 90 degrees"
            if (family == 30) 
                fam <- "Rotated Joe-Frank 90 degrees"
            if (family == 33) 
                fam <- "Rotated Clayton 270 degrees"
            if (family == 34) 
                fam <- "Rotated Gumbel 270 degrees"
            if (family == 36) 
                fam <- "Rotated Joe 270 degrees"
            if (family == 37) 
                fam <- "Rotated Clayton-Gumbel 270 degrees"
            if (family == 38) 
                fam <- "Rotated Joe-Gumbel 270 degrees"
            if (family == 39) 
                fam <- "Rotated Joe-Clayton 270 degrees"
            if (family == 40) 
                fam <- "Rotated Joe-Frank 270 degrees"
            if (family == 100) 
                fam <- "Nonparametric"  #changed Mathias
            if (family == 41) 
                fam <- "1-parametric asymmetric"
            if (family == 51) 
                fam <- "Rotated 1-parametric asymmetric 180 degree"
            if (family == 61) 
                fam <- "Rotated 1-parametric asymmetric 90 degree"
            if (family == 71) 
                fam <- "Rotated 1-parametric asymmetric 270 degree"
            if (family == 104) 
                fam <- "Tawn  type 1"
            if (family == 114) 
                fam <- "Rotated Tawn type 1 180 degrees"
            if (family == 124) 
                fam <- "Rotated Tawn type 1 90 degrees"
            if (family == 134) 
                fam <- "Rotated Tawn type 1 270 degrees"
            if (family == 204) 
                fam <- "Tawn  type 2"
            if (family == 214) 
                fam <- "Rotated Tawn type 2 180 degrees"
            if (family == 224) 
                fam <- "Rotated Tawn type 2 90 degrees"
            if (family == 234) 
                fam <- "Rotated Tawn type 2 270 degrees"
        }
    } else {
        # Name zu Zahl
        if (family == "I" || family == "Independence") 
            fam <- 0
        if (family == "N" || family == "Gaussian") 
            fam <- 1
        if (family == "t") 
            fam <- 2
        if (family == "C" || family == "Clayton") 
            fam <- 3
        if (family == "G" || family == "Gumbel") 
            fam <- 4
        if (family == "F" || family == "Frank") 
            fam <- 5
        if (family == "J" || family == "Joe") 
            fam <- 6
        if (family == "BB1" || family == "Clayton-Gumbel") 
            fam <- 7
        if (family == "BB6" || family == "Joe-Gumbel") 
            fam <- 8
        if (family == "BB7" || family == "Joe-Clayton") 
            fam <- 9
        if (family == "SC" || family == "Survival Clayton") 
            fam <- 13
        if (family == "SG" || family == "Survival Gumbel") 
            fam <- 14
        if (family == "SJ" || family == "Survival Joe") 
            fam <- 16
        if (family == "SBB1" || family == "Survival Clayton-Gumbel") 
            fam <- 17
        if (family == "SBB6" || family == "Survival Joe-Gumbel") 
            fam <- 18
        if (family == "SBB7" || family == "Survival Joe-Clayton") 
            fam <- 19
        if (family == "SBB8" || family == "Survival Joe-Frank") 
            fam <- 20
        if (family == "C90" || family == "Rotated Clayton 90 degrees") 
            fam <- 23
        if (family == "G90" || family == "Rotated Gumbel 90 degrees") 
            fam <- 24
        if (family == "J90" || family == "Rotated Joe 90 degrees") 
            fam <- 26
        if (family == "BB1_90" || family == "Rotated Clayton-Gumbel 90 degrees") 
            fam <- 27
        if (family == "BB6_90" || family == "Rotated Joe-Gumbel 90 degrees") 
            fam <- 28
        if (family == "BB7_90" || family == "Rotated Joe-Clayton 90 degrees") 
            fam <- 29
        if (family == "BB8_90" || family == "Rotated Joe-Frank 90 degrees") 
            fam <- 30
        if (family == "C270" || family == "Rotated Clayton 270 degrees") 
            fam <- 33
        if (family == "G270" || family == "Rotated Gumbel 270 degrees") 
            fam <- 34
        if (family == "J270" || family == "Rotated Joe 270 degrees") 
            fam <- 36
        if (family == "BB1_270" || family == "Rotated Clayton-Gumbel 270 degrees") 
            fam <- 37
        if (family == "BB6_270" || family == "Rotated Joe-Gumbel 270 degrees") 
            fam <- 38
        if (family == "BB7_270" || family == "Rotated Joe-Clayton 270 degrees") 
            fam <- 39
        if (family == "BB8_270" || family == "Rotated Joe-Frank 270 degrees") 
            fam <- 40
        if (family == "NP" || family == "Nonparametric") 
            fam <- 100  #changed Mathias
        if (family == "1-par AS" || family == "1-parametric asymmetric") 
            fam <- 41
        if (family == "1-par AS180" || family == "Rotated 1-parametric asymmetric 180 degree") 
            fam <- 51
        if (family == "1-par AS90" || family == "Rotated 1-parametric asymmetric 90 degree") 
            fam <- 61
        if (family == "1-par AS270" || family == "Rotated 1-parametric asymmetric 270 degree") 
            fam <- 71
        if (family == "Tawn" || family == "Tawn type 1") 
            fam <- 104
        if (family == "Tawn180" || family == "Rotated Tawn type 1 180 degrees") 
            fam <- 114
        if (family == "Tawn90" || family == "Rotated Tawn type 1 90 degrees") 
            fam <- 124
        if (family == "Tawn270" || family == "Rotated Tawn type 1 270 degrees") 
            fam <- 134
        if (family == "Tawn2" || family == "Tawn type 2") 
            fam <- 204
        if (family == "Tawn2_180" || family == "Rotated Tawn type 2 180 degrees") 
            fam <- 214
        if (family == "Tawn2_90" || family == "Rotated Tawn type 2 90 degrees") 
            fam <- 224
        if (family == "Tawn2_270" || family == "Rotated Tawn type 2 270 degrees") 
            fam <- 234
    }
    if (is.na(fam)) 
        stop("Family not implemented.")
    return(fam)
}