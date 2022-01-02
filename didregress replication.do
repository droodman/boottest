use https://www.stata-press.com/data/r17/hospdd, clear

qui didregress (satis) (procedure), group(hospital) time(month) // wildbootstrap(rseed(123) errorweight(webb))
boottest, nogr seed(123)
qui areg satis procedure i.month, a(hospital) cluster(hospital)
boottest procedure, nogr seed(123)

qui didregress (satis) (procedure), group(hospital) time(month) nogteffects
boottest, nogr seed(123)
qui reg satis procedure, cluster(hospital)
boottest procedure, nogr seed(123)

generate hightrt = procedure==1 & (frequency==3 | frequency==4)
label define trt 0 "Untreated" 1 "Treated"
label values hightrt trt
didregress (satis) (hightrt), group(hospital frequency) time(month) // wildbootstrap(rseed(123) errorweight(webb))
boottest, nogr seed(123)
qui areg satis hightrt month#frequency month#hospital frequency#hospital, cluster(hospital) a(hospital)
boottest hightrt, nogr seed(123)

didregress (satis) (hightrt), group(hospital frequency) time(month) nogteffects  // no change in estimate but more entries of e(b) non-zero
boottest, nogr seed(123)
qui reg satis hightrt month#frequency month#hospital frequency#hospital, cluster(hospital)
boottest hightrt, nogr seed(123)

didregress (satis) (hightrt), group(hospital frequency) time(month) nointeract
boottest, nogr seed(123)
qui areg satis hightrt i.month i.frequency, cluster(hospital) a(hospital) 
boottest hightrt, nogr seed(123)

use https://www.stata-press.com/data/r17/patents, clear
qui xtdidregress (uspatents fpatents) (gotpatent), group(classid) time(year) // wildbootstrap(rseed(123) errorweight(webb) reps(99))
boottest, nogr seed(123) reps(99)
qui xtreg uspatents gotpatent fpatents i.year, fe cluster(classid)
boottest gotpatent, nogr seed(123) reps(99)
