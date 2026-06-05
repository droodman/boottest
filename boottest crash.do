webuse laborsup, clear
ppmlhdfe kids fem_inc male_educ, a(fem_work) d
predict double mu_hat, mu
gen double z_star = ln(mu_hat) + (kids - mu_hat) / mu_hat
reghdfe z_star fem_inc male_educ [aw=mu_hat], absorb(fem_work)

// net install boottest, replace from(https://raw.github.com/droodman/boottest/v4.5.2)
// run "c:\ado\plus\b\boottest 4.5.2.ado"
// boottest fem_inc, noci

// net install boottest, replace from(https://raw.github.com/droodman/boottest/v4.5.3)
run "c:\ado\plus\b\boottest 4.5.3.ado"
boottest fem_inc, noci
