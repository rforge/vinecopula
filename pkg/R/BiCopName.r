BiCopName <- function(family, short=TRUE)
{
if(is.logical(short)==FALSE) stop("'short' has to be a logical variable.")

if(is.numeric(family))	# Zahl zu Name
{
	if(short==TRUE)		# kurzer Name
	{
		if(family==0) fam="I"
		else if(family==1) fam="N"
		else if(family==2) fam="t"
		else if(family==3) fam="C"
		else if(family==4) fam="G"
		else if(family==5) fam="F"
		else if(family==6) fam="J"
		else if(family==7) fam="BB1"
		else if(family==8) fam="BB6"
		else if(family==9) fam="BB7"
		else if(family==10) fam="BB8"
		else if(family==13) fam="SC"
		else if(family==14) fam="SG"
		else if(family==16) fam="SJ"
		else if(family==17) fam="SBB1"
		else if(family==18) fam="SBB6"
		else if(family==19) fam="SBB7"
		else if(family==20) fam="SBB8"
		else if(family==23) fam="C90"
		else if(family==24) fam="G90"
		else if(family==26) fam="J90"
		else if(family==27) fam="BB1_90"
		else if(family==28) fam="BB6_90"
		else if(family==29) fam="BB7_90"
		else if(family==30) fam="BB8_90"
		else if(family==33) fam="C270"
		else if(family==34) fam="G270"
		else if(family==36) fam="J270"
		else if(family==37) fam="BB1_270"
		else if(family==38) fam="BB6_270"
		else if(family==39) fam="BB7_270"
		else if(family==40) fam="BB8_270"
		else if(family==100) fam="NP"    #changed Mathias
		else if(family==41) fam="1-par AS"
		else if(family==51) fam="1-par AS180"
		else if(family==61) fam="1-par AS90"
		else if(family==71) fam="1-par AS270"
		else if(family==104) fam="Tawn"
		else if(family==114) fam="Tawn180"
		else if(family==124) fam="Tawn90"
		else if(family==134) fam="Tawn270"
		else if(family==204) fam="Tawn2"
		else if(family==214) fam="Tawn2_180"
		else if(family==224) fam="Tawn2_90"
		else if(family==234) fam="Tawn2_270"
		else stop("Family not implemented.")
	}
	else			# langer Name
	{
		if(family==0) fam="Independence"
		else if(family==1) fam="Gaussian"
		else if(family==2) fam="t"
		else if(family==3) fam="Clayton"
		else if(family==4) fam="Gumbel"
		else if(family==5) fam="Frank"
		else if(family==6) fam="Joe"
		else if(family==7) fam="Clayton-Gumbel"
		else if(family==8) fam="Joe-Gumbel"
		else if(family==9) fam="Joe-Clayton"
		else if(family==10) fam="Joe-Frank"
		else if(family==13) fam="Survival Clayton"
		else if(family==14) fam="Survival Gumbel"
		else if(family==16) fam="Survival Joe"
		else if(family==17) fam="Survival Clayton-Gumbel"
		else if(family==18) fam="Survival Joe-Gumbel"
		else if(family==19) fam="Survival Joe-Clayton"
		else if(family==20) fam="Survival Joe-Frank"
		else if(family==23) fam="Rotated Clayton 90 degrees"
		else if(family==24) fam="Rotated Gumbel 90 degrees"
		else if(family==26) fam="Rotated Joe 90 degrees"
		else if(family==27) fam="Rotated Clayton-Gumbel 90 degrees"
		else if(family==28) fam="Rotated Joe-Gumbel 90 degrees"
		else if(family==29) fam="Rotated Joe-Clayton 90 degrees"
		else if(family==30) fam="Rotated Joe-Frank 90 degrees"
		else if(family==33) fam="Rotated Clayton 270 degrees"
		else if(family==34) fam="Rotated Gumbel 270 degrees"
		else if(family==36) fam="Rotated Joe 270 degrees"
		else if(family==37) fam="Rotated Clayton-Gumbel 270 degrees"
		else if(family==38) fam="Rotated Joe-Gumbel 270 degrees"
		else if(family==39) fam="Rotated Joe-Clayton 270 degrees"
		else if(family==40) fam="Rotated Joe-Frank 270 degrees"
		else if(family==100) fam="Nonparametric"  #changed Mathias
		else if(family==41) fam="1-parametric asymmetric"
		else if(family==51) fam="Rotated 1-parametric asymmetric 180 degree"
		else if(family==61) fam="Rotated 1-parametric asymmetric 90 degree"
		else if(family==71) fam="Rotated 1-parametric asymmetric 270 degree"
		else if(family==104) fam="Tawn"
		else if(family==114) fam="Rotated Tawn 180 degrees"
		else if(family==124) fam="Rotated Tawn 90 degrees"
		else if(family==134) fam="Rotated Tawn 270 degrees"
		else if(family==204) fam="Tawn2"
		else if(family==214) fam="Rotated Tawn2 180 degrees"
		else if(family==224) fam="Rotated Tawn2 90 degrees"
		else if(family==234) fam="Rotated Tawn2 270 degrees"
		else stop("Family not implemented.")
	}
}
else	# Name zu Zahl
{
	if(family=="I" || family=="Independence") fam=0
	else if(family=="N" || family=="Gaussian") fam=1
	else if(family=="t") fam=2
	else if(family=="C" || family=="Clayton") fam=3
	else if(family=="G" || family=="Gumbel") fam=4
	else if(family=="F" || family=="Frank") fam=5
	else if(family=="J" || family=="Joe") fam=6
	else if(family=="BB1" || family=="Clayton-Gumbel") fam=7
	else if(family=="BB6" || family=="Joe-Gumbel") fam=8
	else if(family=="BB7" || family=="Joe-Clayton") fam=9
	else if(family=="SC" || family=="Survival Clayton") fam=13
	else if(family=="SG" || family=="Survival Gumbel") fam=14
	else if(family=="SJ" || family=="Survival Joe") fam=16
	else if(family=="SBB1" || family=="Survival Clayton-Gumbel") fam=17
	else if(family=="SBB6" || family=="Survival Joe-Gumbel") fam=18
	else if(family=="SBB7" || family=="Survival Joe-Clayton") fam=19
	else if(family=="SBB8" || family=="Survival Joe-Frank") fam=20
	else if(family=="C90" || family=="Rotated Clayton 90 degrees") fam=23
	else if(family=="G90" || family=="Rotated Gumbel 90 degrees") fam=24
	else if(family=="J90" || family=="Rotated Joe 90 degrees") fam=26
	else if(family=="BB1_90" || family=="Rotated Clayton-Gumbel 90 degrees") fam=27
	else if(family=="BB6_90" || family=="Rotated Joe-Gumbel 90 degrees") fam=28
	else if(family=="BB7_90" || family=="Rotated Joe-Clayton 90 degrees") fam=29
	else if(family=="BB8_90" || family=="Rotated Joe-Frank 90 degrees") fam=30
	else if(family=="C270" || family=="Rotated Clayton 270 degrees") fam=33
	else if(family=="G270" || family=="Rotated Gumbel 270 degrees") fam=34
	else if(family=="J270" || family=="Rotated Joe 270 degrees") fam=36
	else if(family=="BB1_270" || family=="Rotated Clayton-Gumbel 270 degrees") fam=37
	else if(family=="BB6_270" || family=="Rotated Joe-Gumbel 270 degrees") fam=38
	else if(family=="BB7_270" || family=="Rotated Joe-Clayton 270 degrees") fam=39
	else if(family=="BB8_270" || family=="Rotated Joe-Frank 270 degrees") fam=40
	else if(family=="NP" || family=="Nonparametric") fam=100 #changed Mathias
	else if(family=="1-par AS" || family=="1-parametric asymmetric") fam=41
	else if(family=="1-par AS180" || family=="Rotated 1-parametric asymmetric 180 degree") fam=51
	else if(family=="1-par AS90" || family=="Rotated 1-parametric asymmetric 90 degree") fam=61
	else if(family=="1-par AS270" || family=="Rotated 1-parametric asymmetric 270 degree") fam=71
	else if(family=="Tawn") fam=104
	else if(family=="Tawn180" || family=="Rotated Tawn 180 degrees") fam=114
	else if(family=="Tawn90" || family=="Rotated Tawn 90 degrees") fam=124
	else if(family=="Tawn270" || family=="Rotated Tawn 270 degrees") fam=134
	else if(family=="Tawn2") fam=204
	else if(family=="Tawn2_180" || family=="Rotated Tawn2 180 degrees") fam=214
	else if(family=="Tawn2_90" || family=="Rotated Tawn2 90 degrees") fam=224
	else if(family=="Tawn2_270" || family=="Rotated Tawn2 270 degrees") fam=234
	else stop("Family not implemented.")
}

return(fam)
}