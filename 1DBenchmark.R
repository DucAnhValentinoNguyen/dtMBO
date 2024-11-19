library(mlrMBO)
library(lhs)
library(DiceKriging)

obj.fun = makeSingleObjectiveFunction(
  name = "my_black_box",
  fn = function(x)
    (x[2] - 0.1 * x[1] ** 2 + x[1] - 6) ** 2 + cos(x[1]),
  par.set = makeParamSet(
    makeNumericParam("x1", lower = -5, upper = 10),
    makeNumericParam("x2", lower = 0, upper = 15)
  )
)

# sym um x = 0
obj.fun = makeAckleyFunction(2) # (0,0,0)
# sym um x = 0
obj.fun = makeGriewankFunction(1) # (0,0,0)
obj.fun = makeAlpine01Function(1) # (0,0,0)
# sym um x = 5
obj.fun = makeDeflectedCorrugatedSpringFunction(1) #(5, 5,-1)
obj.fun = makeSchwefelFunction(1) # (420.9687, 420.9687, -418.9829)
obj.fun = makeRosenbrockFunction(1) # (1,1,0), only take dim >= 2
obj.fun

ggplot2::autoplot(obj.fun, render.levels = TRUE, show.optimum = TRUE)

# obj.fun = makeBraninFunction()
# obj.fun = makeSixHumpCamelFunction()
# obj.fun = makeEggholderFunction()


set.seed(123)
# It is recommended to use a Latin Hypercube Design by calling generateDesign() and
# passing the number of desired points. If no design is given by the user, mlrMBO
# will generate a maximin Latin Hypercube Design of size 4 times the number of the
# black - box function’s parameters.

# Create initial random Latin Hypercube Design of 10 points
# design = generateDesign(n = 5L * 2L, getParamSet(fn), fun = randomLHS)
design = generateDesign(5 * getNumberOfParameters(obj.fun),
                        getParamSet(obj.fun),
                        fun = lhs::maximinLHS)
design

# For surrogate regression, Kriging (makeLearner("regr.km")) and random forests
# (makeLearner("regr.randomForest"))

# If no regression method is supplied
# by the user, the fallback is a Kriging model with a Matern-3/2 kernel and the
# “GENetic Optimization Using Derivatives” (genoud) fitting algorithm in a fully
# numeric setting, and a random forest with jackknife variance estimation otherwise.

# if expected improvement or LCB is chosen as the infill criterion, the surrogate
# either has to provide an uncertainty estimator, or has to be combined with a bagging
# approach using the makeBaggingWrapper() in mlr.


# Specify kriging model with standard error estimation
surrogate = makeLearner("regr.km", predict.type = "se",
                        covtype = "matern3_2")
surrogate

# Set general controls
ctrl = makeMBOControl()
ctrl = setMBOControlTermination(ctrl, iters = 20 * getNumberOfParameters(obj.fun))

# If the infill optimization is unspecified, mlrMBO uses LCB as infill criterion
# with λ = 1 in a fully numeric setting and λ = 2 if at least one discrete
# parameter is present.

# LCB
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritCB())

# EI
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())

# MaxMinCB
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritMaxMinCB())

# MaxMaxCB
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritMaxMaxCB())

# MinMaxCB
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritMinMaxCB())



# Hurwicz
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritHurwicz())

# HodgesLehmann
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritHodgesLehmann())

# Bayes Decision
ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritBayesDecision())


# To optimize the criterion, focus search with nrestarts = 3, niters = 5 and
# npoints = 1000 is used by default.

run = exampleRun(obj.fun,
                 design = design,
                 learner = surrogate,
                 control = ctrl)
plotExampleRun(run, pause = FALSE)


# start mbo
obj.fun
ctrl

result = mbo(obj.fun , design = design, control = ctrl)
print(result)
plot(result)
resultAckley <-  result
resultAckley