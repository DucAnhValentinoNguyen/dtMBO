# Einzelne Ergebnisse der synthetischen Funktionen
# Experimente mit synthetischen Funktionen
# Bei den Benchmark Experimenten mit synthetischen Funktionen wurden sowohl
# Funktionen mit stetigem Parameterraum, als auch mit gemischtem Parameterraum
# ausgewählt. Die Auswirkung der Dimension auf die Methoden ließ sich mit den
# Funktionen mit stetigem Parameterraum untersuchen.
library(mlrMBO)
library(batchtools)
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)

reg = makeExperimentRegistry("benchmarkInfill")

# addProblem(name = "sphere",
#            data = # smooth, unimodal, sym um x = 0
#              makeSphereFunction(2))

# addProblem(name = "1Dgriewank",
#            data = makeGriewankFunction(2)) # (0,0,0)) # (0,0,0))


#(5, 5,-1)
addProblem(name = "2DDeflectedCorrugatedSpring",
           data = makeDeflectedCorrugatedSpringFunction(2))

# schwefel range too wide, less than rosenbrock though
addProblem(name = "2DSchwefel",
           data = makeSchwefelFunction(2)) # (420.9687, 420.9687, -418.9829))

# rosenborck range too wide
addProblem(name = "2DRosenbrock",
           data = # sym um x = 0
             # multimodal, continuous, differentiable
             makeRosenbrockFunction(2)) # (0,0,0))





# addProblem(name = "4DSphere",
#            data = makeSphereFunction(4))

#(5, 5,-1)
addProblem(name = "4DDeflectedCorrugatedSpring",
           data = makeDeflectedCorrugatedSpringFunction(4))

# schwefel range too wide, less than rosenbrock though
# (420.9687, 420.9687, -418.9829)
addProblem(name = "4DSchwefel",
           data = makeSchwefelFunction(4))

# Global optimum objective value of -3.3224 at
# x1       x2       x3       x4       x5     x6
# 1 0.20169 0.150011 0.476874 0.275332 0.311652 0.6573
addProblem(name = "4DHartmann",
           data = makeHartmannFunction(4))




#(5, 5,-1)
addProblem(name = "6DDeflectedCorrugatedSpring",
           data = makeDeflectedCorrugatedSpringFunction(6))

addProblem(name = "6DSchwefel",
           data = makeSchwefelFunction(6))

addProblem(name = "6DRosenbrock",
           data = makeRosenbrockFunction(6))

addProblem(name = "6DHartmann",
           data = makeHartmannFunction(6))



mbo_optimization <-
  function(job, data, instance, infill_crit, ...) {
    # Define the design based on the parameter set of the instance
    design = generateDesign(5 * getNumberOfParameters(instance),
                            getParamSet(instance),
                            fun = lhs::maximinLHS)

    # Create the surrogate model using kriging with Matern 3/2 covariance
    surrogate = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")

    # Set up MBO control with termination after 200 iterations
    ctrl = makeMBOControl()
    ctrl = setMBOControlTermination(ctrl, iters = 200)

    # Set the specified infill criterion
    ctrl = setMBOControlInfill(ctrl, crit = infill_crit)

    # Run the optimization process
    result = mbo(instance, design = design, control = ctrl)
    run.time = result$final.opt.state$time.used
    best.y = result$y
    return(list(best.y, run.time))
  }


# For UCB
addAlgorithm(
  name = "ucb",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritCB(), ...)
)

# For EI
addAlgorithm(
  name = "ei",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritEI(), ...)
)

# For Maximin CB
addAlgorithm(
  name = "maximin",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritMaximinCB(), ...)
)

addAlgorithm(
  name = "minimax",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritMinimaxCB(), ...)
)


# addAlgorithm(name = "hurwicz", fun = function(job, data, instance, ...)
#   mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritHurwicz(), ...)
# )

addAlgorithm(
  name = "hurwiczAlpha0.7",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritHurwicz(0.7), ...)
)

# addAlgorithm(name = "hodgesLehmann", fun = function(job, data, instance, ...)
#   mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritHodgesLehmann(), ...)
# )

addAlgorithm(
  name = "hodgesLehmannFactor0.8",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritHodgesLehmann(0.8), ...)
)



# addAlgorithm(
#   name = "ucb",
#   fun = function(job, data, instance, ...) {
#     design = generateDesign(5 * getNumberOfParameters(instance),
#                             getParamSet(instance),
#                             fun = lhs::maximinLHS)
#
#     # Create the surrogate model using kriging
#     surrogate = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")
#
#     # Set up MBO control with LCB and termination after 200 iterations
#     ctrl = makeMBOControl()
#     ctrl = setMBOControlTermination(ctrl, iters = 200)
#     ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritCB())
#
#     # Run the optimization process
#     result = mbo(instance, design = design, control = ctrl)
#     run.time = result$final.opt.state$time.used
#     best.y = result$y
#     return(list(best.y, run.time))
#   }
# )
#
#
# addAlgorithm(
#   name = "ei",
#   fun = function(job, data, instance, ...) {
#     design = generateDesign(5 * getNumberOfParameters(instance),
#                             getParamSet(instance),
#                             fun = lhs::maximinLHS)
# 
#     # Create the surrogate model using kriging
#     surrogate = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")
# 
#     # Set up MBO control with LCB and termination after 200 iterations
#     ctrl = makeMBOControl()
#     ctrl = setMBOControlTermination(ctrl, iters = 200)
#     ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())
# 
#     # Run the optimization process
#     result = mbo(instance, design = design, control = ctrl)
#     run.time = result$final.opt.state$time.used
#     best.y = result$y
#     return(list(best.y, run.time))
#   }
# )
# 
# addAlgorithm(
#   name = "maximin",
#   fun = function(job, data, instance, ...) {
#     design = generateDesign(5 * getNumberOfParameters(instance),
#                             getParamSet(instance),
#                             fun = lhs::maximinLHS)
# 
#     # Create the surrogate model using kriging
#     surrogate = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")
# 
#     # Set up MBO control with LCB and termination after 200 iterations
#     ctrl = makeMBOControl()
#     ctrl = setMBOControlTermination(ctrl, iters = 200)
#     ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritMaximinCB())
# 
#     # Run the optimization process
#     result = mbo(instance, design = design, control = ctrl)
#     run.time = result$final.opt.state$time.used
#     best.y = result$y
#     return(list(best.y, run.time))
#   }
# )
# 
# 
# addAlgorithm(
#   name = "minimax",
#   fun = function(job, data, instance, ...) {
#     design = generateDesign(5 * getNumberOfParameters(instance),
#                             getParamSet(instance),
#                             fun = lhs::maximinLHS)
# 
#     # Create the surrogate model using kriging
#     surrogate = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")
# 
#     # Set up MBO control with LCB and termination after 200 iterations
#     ctrl = makeMBOControl()
#     ctrl = setMBOControlTermination(ctrl, iters = 200)
#     ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritMinimaxCB())
# 
#     # Run the optimization process
#     result = mbo(instance, design = design, control = ctrl)
#     run.time = result$final.opt.state$time.used
#     best.y = result$y
#     return(list(best.y, run.time))
#   }
# )


# addAlgorithm(
#   name = "hurwicz",
#   fun = function(job, data, instance, ...) {
#     design = generateDesign(5 * getNumberOfParameters(instance),
#                             getParamSet(instance),
#                             fun = lhs::maximinLHS)
#
#     # Create the surrogate model using kriging
#     surrogate = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")
#
#     # Set up MBO control with LCB and termination after 200 iterations
#     ctrl = makeMBOControl()
#     ctrl = setMBOControlTermination(ctrl, iters = 200)
#     ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritHurwicz())
#
#     # Run the optimization process
#     result = mbo(instance, design = design, control = ctrl)
#     run.time = result$final.opt.state$time.used
#     best.y = result$y
#     return(list(best.y, run.time))
#   }
# )

# addAlgorithm(
#   name = "hurwiczalpha0.7",
#   fun = function(job, data, instance, ...) {
#     design = generateDesign(5 * getNumberOfParameters(instance),
#                             getParamSet(instance),
#                             fun = lhs::maximinLHS)
# 
#     # Create the surrogate model using kriging
#     surrogate = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")
# 
#     # Set up MBO control with LCB and termination after 200 iterations
#     ctrl = makeMBOControl()
#     ctrl = setMBOControlTermination(ctrl, iters = 200)
#     ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritHurwicz(0.7))
# 
#     # Run the optimization process
#     result = mbo(instance, design = design, control = ctrl)
#     run.time = result$final.opt.state$time.used
#     best.y = result$y
#     return(list(best.y, run.time))
#   }
# )

# addAlgorithm(
#   name = "hodgesLehmann",
#   fun = function(job, data, instance, ...) {
#     design = generateDesign(5 * getNumberOfParameters(instance),
#                             getParamSet(instance),
#                             fun = lhs::maximinLHS)
#
#     # Create the surrogate model using kriging
#     surrogate = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")
#
#     # Set up MBO control with LCB and termination after 200 iterations
#     ctrl = makeMBOControl()
#     ctrl = setMBOControlTermination(ctrl, iters = 200)
#     ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritHodgesLehmann())
#
#     # Run the optimization process
#     result = mbo(instance, design = design, control = ctrl)
#     run.time = result$final.opt.state$time.used
#     best.y = result$y
#     return(list(best.y, run.time))
#   }
# )

# addAlgorithm(
#   name = "hodgesLehmannFactor0.7",
#   fun = function(job, data, instance, ...) {
#     design = generateDesign(5 * getNumberOfParameters(instance),
#                             getParamSet(instance),
#                             fun = lhs::maximinLHS)
# 
#     # Create the surrogate model using kriging
#     surrogate = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")
# 
#     # Set up MBO control with LCB and termination after 200 iterations
#     ctrl = makeMBOControl()
#     ctrl = setMBOControlTermination(ctrl, iters = 200)
#     ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritHodgesLehmann(0.7))
# 
#     # Run the optimization process
#     result = mbo(instance, design = design, control = ctrl)
#     run.time = result$final.opt.state$time.used
#     best.y = result$y
#     return(list(best.y, run.time))
#   }
# )
# 



# Overview over defined experiments
reg$problems
reg$algorithms

# add experiments and submit
ids = addExperiments(reg = reg, repls = 1)
summarizeExperiments(reg = reg)
# summarizeExperiments(reg = tmp, by = c("problem", "algorithm", "n"))
# ids = findExperiments(prob.pars = (n == 50), reg = tmp)
# print(unwrap(getJobPars(ids, reg = tmp)))
# Submit jobs
submitJobs()
waitForJobs()

# # Reduce the results of algorithm ucb
# ids.ucb = findExperiments(algo.name = "ucb", reg = reg)
# reduceResults(
#   ids.ucb,
#   fun = function(aggr, res, ...)
#     c(aggr, res),
#   reg = reg
# )
# # # Join info table with all results and calculate mean of results
# # # grouped by n and algorithm
# pars = unwrap(getJobPars(ids, reg = reg))
#
# results = unwrap(reduceResultsDataTable(
#   ids,
#   fun = function(res)
#     list(res = res),
#   reg = reg
# ))

res = ijoin(unwrap(getJobPars(ids, reg = reg)),
            reduceResultsDataTable(
              reg = reg,
              fun = function(x)
                list(res = x)
            ))
res <- res |> as.data.frame()

# Step 1: Extract values from the nested list
res_values <- lapply(res$result, function(x) {
  best_y <- x$res[[1]]  # Extract res[[1]]
  running_time <-
    as.numeric(x$res[[2]])  # Extract res[[2]] and convert to numeric
  return(c(best_y, running_time))
})

# Step 2: Convert the list to a data frame
res_df <- do.call(rbind, res_values)  # Combine all rows
colnames(res_df) <- c("best.y", "running.time")  # Name the columns
res <- cbind(res, res_df)
res$result = NULL





plot <- ggplot(res, aes(x = algorithm, y = best.y)) +
  geom_boxplot(aes(fill = algorithm)) +
  facet_wrap(~ problem, ncol = 2, scales = "free") +
  labs(title = "Best objective value (on y axis) found by respective algorithms on respective test
function",
       x = "algorithms",
       y = "best y") + theme_bw()

plot
ggsave(path = "Plots",
       filename = "Benchmark1.png")

avg_running_time <- res %>%
  group_by(problem, algorithm) %>%
  summarise(avg_time = mean(running.time))

# Print the result
avg_running_time



# Add a rank column based on the best.y values
res <- res %>%
  mutate(rank = dense_rank(best.y))



avg_rank <- res %>%
  group_by(algorithm) %>%
  summarise(avg_rank = mean(rank))

avg_rank

res

res <- res %>%
  left_join(avg_running_time, by = c("problem", "algorithm")) |>
  left_join(avg_rank, by = c("problem", "algorithm"))


# Print the final result
res

#################################


# 
# tmp = makeExperimentRegistry(file.dir = NA, make.default = FALSE)
# # add first problem
# fun = function(job, data, n, mean, sd, ...)
#   rnorm(n, mean = mean, sd = sd)
# addProblem("rnorm", fun = fun, reg = tmp)
# # add second problem
# fun = function(job, data, n, lambda, ...)
#   rexp(n, rate = lambda)
# addProblem("rexp", fun = fun, reg = tmp)
# # add first algorithm
# fun = function(instance, method, ...)
#   if (method == "mean")
#     mean(instance)
# else
#   median(instance)
# addAlgorithm("average", fun = fun, reg = tmp)
# # add second algorithm
# fun = function(instance, ...)
#   sd(instance)
# addAlgorithm("deviation", fun = fun, reg = tmp)
# # define problem and algorithm designs
# library(data.table)
# prob.designs = algo.designs = list()
# prob.designs$rnorm = CJ(n = 100, mean = -1:1, sd = 1:5)
# prob.designs$rexp = data.table(n = 100, lambda = 1:5)
# algo.designs$average = data.table(method = c("mean", "median"))
# algo.designs$deviation = data.table()
# # add experiments and submit
# addExperiments(prob.designs, algo.designs, reg = tmp)
# # check what has been created
# summarizeExperiments(reg = tmp)
# unwrap(getJobPars(reg = tmp))



#################################
# reg = makeExperimentRegistry(file.dir = NA, make.default = FALSE)
# # Define one problem, two algorithms and add them with some parameters:
# addProblem(
#   reg = reg,
#   "p1",
#   fun = function(job, data, n, mean, sd, ...)
#     rnorm(n, mean = mean, sd = sd)
# )
# addAlgorithm(
#   reg = reg,
#   "a1",
#   fun = function(job, data, instance, ...)
#     mean(instance)
# )
# addAlgorithm(
#   reg = reg,
#   "a2",
#   fun = function(job, data, instance, ...)
#     median(instance)
# )
# ids = addExperiments(reg = reg, list(p1 = data.table::CJ(
#   n = c(50, 100),
#   mean = -2:2,
#   sd = 1:4
# )))
# # Overview over defined experiments
# tmp$problems
# tmp$algorithms
# summarizeExperiments(reg = tmp)
# summarizeExperiments(reg = tmp, by = c("problem", "algorithm", "n"))
# ids = findExperiments(prob.pars = (n == 50), reg = tmp)
# print(unwrap(getJobPars(ids, reg = tmp)))
# # Submit jobs
# submitJobs(reg = tmp)
# waitForJobs(reg = tmp)
# # Reduce the results of algorithm a1
# ids.mean = findExperiments(algo.name = "a1", reg = tmp)
# reduceResults(
#   ids.mean,
#   fun = function(aggr, res, ...)
#     c(aggr, res),
#   reg = tmp
# )
# # Join info table with all results and calculate mean of results
# # grouped by n and algorithm
# ids = findDone(reg = tmp)
# pars = unwrap(getJobPars(ids, reg = tmp))
# results = unwrap(reduceResultsDataTable(
#   ids,
#   fun = function(res)
#     list(res = res),
#   reg = tmp
# ))
# tab = ljoin(pars, results)
# tab[, list(mres = mean(res)), by = c("n", "algorithm")]
