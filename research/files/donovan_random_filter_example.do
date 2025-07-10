* Pierce Donovan, pierce.donovan@unr.edu
* Random Filter Identification Strategy: Simulated Example
* February 2024

* Just some code to replicate the simulation results from the introduction
* of my Random Filter lecture. Replicates the four tables that demonstrate
* identification of the ATT amid selection into treatment. I also added 
* some bootstrap code at the bottom that grabs the standard error in the
* base simulation. It can easily be applied to the other simulations.

****************************************************

* base simulation

clear
set seed 3.14
set obs 10000
gen U = rbinomial(1,.5) // financial savvy
gen T = rbinomial(1,.25+.5*U) // sign up for account
gen M = rbinomial(1,.67*T) // usable account
replace M = 0 if M==. // fix 0 probability case
gen Y = rnormal(500+500*U+300*M,25) // yearly savings

* naive and correct specifications
qui reg Y T
est sto naive

qui reg Y T U
est sto benchmark

* random filter gives right answer
qui reg Y M T
est sto second

qui reg M T
est sto first

* full regression shows T picks up only the confounding variation
qui reg Y M T U
est sto full

* put everything in one table
est tab naive benchmark first second full, b se

****************************************************

* multiple causal channels

clear
set seed 3.14
set obs 10000
gen U = rbinomial(1,.5) // financial savvy
gen T = rbinomial(1,.25+.5*U) // sign up for account
gen M = rbinomial(1,.67*T) // usable account
replace M = 0 if M==. // fix 0 probability case
gen Y = rnormal(500+500*U+300*M+30*T,25) // yearly savings

* naive and correct specifications
qui reg Y T
est sto naive

qui reg Y T U
est sto benchmark

* random filter gives right answer
qui reg Y M T
est sto second

qui reg M T
est sto first

* full regression shows both effects separated
qui reg Y M T U
est sto full

* put everything in one table
est tab naive benchmark first second full, b se

****************************************************

* conditional exogeneity

clear
set seed 3.14
set obs 10000
gen U = rbinomial(1,.5) // financial savvy
gen T = rbinomial(1,.25+.5*U) // sign up for account
gen F = rbinomial(1,.25+.5*U) // knows how to use
gen M = rbinomial(1,(.33+.67*K)*T) // usable account
replace M = 0 if M==. // fix 0 probability case
replace M = 1 if T==1 & K==1 // fix 1 probability case
gen Y = rnormal(500+500*U+300*M,25) // yearly savings

* naive and correct specifications
qui reg Y T
est sto naive

qui reg Y T U
est sto benchmark

* unconditional filter, the doubly-biased one
qui reg Y M T
est sto second_unc

qui reg M T
est sto first_unc

* conditional filter, the correct one
qui reg Y M F T
est sto second_con

qui reg M T F
est sto first_con

* put everything in one table
est tab naive benchmark first_unc second_unc first_con second_con, b se

****************************************************

* compare to IV

clear
set seed 3.14
set obs 10000
gen U = rbinomial(1,.5) // financial savvy
gen Z = rbinomial(1,.5) // random assignment instrument
gen C = rbinomial(1,.5) // complier dummy
gen T = (C==0)*rbinomial(1,.25+0.5*U) + (C==1)*Z // sign up for account
gen M = rbinomial(1,.67*T) // usable account
replace M = 0 if M==. // fix 0 probability case
gen Y = rnormal(500+500*U+(450-300*C)*M,25) // yearly savings
* last line creates different impact for compliers vs whole population

* benchmarking
qui reg Y T
est sto naive

qui reg Y T##C U
est sto benchmark

qui reg Y T U if C==1
est sto complier

* RF for ATT
qui reg Y M T
est sto second

qui reg M T
est sto first

* IV for LATE
qui reg Y Z
est sto iv_two

qui reg T Z
est sto iv_one

* IV comparison
est tab naive benchmark complier first second iv*, b se

****************************************************

* bootstrap standard error for base simulation

clear
set seed 3.14
set obs 10000
gen U = rbinomial(1,.5) // financial savvy
gen T = rbinomial(1,.25+.5*U) // sign up for account
gen M = rbinomial(1,.67*T) // usable account
replace M = 0 if M==. // fix 0 probability case
gen Y = rnormal(500+500*U+300*M,25) // yearly savings

* first, the program that the bootstrap references
* refer to this routine as "surboot", stores output in r()
capture program drop surboot
program surboot, rclass
	
	*compatibility
	version 14
	
	* SUR
	sureg (M T) (Y M T)
	
	* extract the two coefficients of interest
	local firststage = _b[M: T]
	local secondstage = _b[Y: M]
	
	* program ouputs the multiplication of the two
    return scalar beta_rf = `firststage'*`secondstage'
	
end

* do the bootstrap with this function
bootstrap r(beta_rf), reps(1000) seed(89509): surboot

