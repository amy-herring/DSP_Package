

# Check if input is valid ------------------------------------------------------
#
# no copies of id in baseline
# no copies of id+cyc in cycle
# no copies of id+cyc+day in daily
#
# baseline/cycle are NULL/mat/arr/df
# daily is mat/arr/df
#
# X, Y TRUE/FALSE or 0/1 or yes/no
#
# (baseline == NULL) <==> (varInclNames$baseline == NULL and similarly for cycle
# varInclNames all match a name in corresp dataset and similarly for other *Name
# no multiple pregnancies
# no multiple id/cyc in cycle or id/cyc/day in day
# if preg in daily then preg consistent throughout cycle
#
# cycle and cycleDay are numeric
# no duplicates in fwDays