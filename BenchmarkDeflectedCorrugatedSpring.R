library(mlrMBO)
library(batchtools)
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(rgenoud)
library(parallelMap)
library(profvis)



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
    #parallelStartMulticore(detectCores() - 1)  # Use all but one core
    result = mbo(instance, design = design, control = ctrl)
    #parallelStop()
    run.time = result$final.opt.state$time.used
    best.y = result$y
    return(list(best.y, run.time))
  }

reg = makeExperimentRegistry(file.dir = NA, make.default = FALSE)

# tmp = makeExperimentRegistry("benchmarkInfill")


# add first problem
# fun = function(job, data, n, mean, sd, ...) rnorm(n, mean = mean, sd = sd)
# addProblem("rnorm", fun = fun, reg = tmp)


addProblem(name = "2DDeflectedCorrugatedSpring",
           data = makeDeflectedCorrugatedSpringFunction(2),
           reg = reg)

addProblem(name = "4DDeflectedCorrugatedSpring",
           data = makeDeflectedCorrugatedSpringFunction(4),
           reg = reg)

addProblem(name = "6DDeflectedCorrugatedSpring",
           data = makeDeflectedCorrugatedSpringFunction(6),
           reg = reg)


# # Define dimensions and functions
# dimensions <- c(2, 4, 6)
# functions <- list(
#   DeflectedCorrugatedSpring = makeDeflectedCorrugatedSpringFunction,
#   Schwefel = makeSchwefelFunction,
#   Rosenbrock = makeRosenbrockFunction,
#   Ackley = makeAckleyFunction
# )
# 
# 
# 
# # Loop through dimensions and functions to add problems
# for (dim in dimensions) {
#   for (func_name in names(functions)) {
#     addProblem(
#       reg = reg,
#       name = paste0(func_name, dim, "D"),
#       data = functions[[func_name]](dim)
#     )
#   }
# }





# fun = function(instance, ...) sd(instance)
# addAlgorithm("deviation", fun = fun, reg = tmp)






# addAlgorithm(
#   reg = tmp,
#   name = "maximin",
#   fun = function(job, data, instance, infill_crit, ...) {
#     # Define the design based on the parameter set of the instance
#     design = generateDesign(5 * getNumberOfParameters(instance),
#                             getParamSet(instance),
#                             fun = lhs::maximinLHS)
#
#     # Create the surrogate model using kriging with Matern 3/2 covariance
#     surrogate = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2")
#
#     # Set up MBO control with termination after 200 iterations
#     ctrl = makeMBOControl()
#     ctrl = setMBOControlTermination(ctrl, iters = 200)
#
#     # Set the specified infill criterion
#     ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritMaximin())
#
#     # Run the optimization process
#     mbo(instance, design = design, control = ctrl)
#     #run.time = result$final.opt.state$time.used
#     #best.y = result$y
#     #return(list(best.y, run.time))
#   }
# )

addAlgorithm(
  reg = reg,
  name = "maximax",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritMaximax(), ...)
)

addAlgorithm(
  reg = reg,
  name = "maximin",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritMaximin(), ...)
)

addAlgorithm(
  reg = reg,
  name = "hurwicz",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritHurwicz(hurwicz.lambda =  0.7), ...)
)

addAlgorithm(
  reg = reg,
  name = "bayes",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritBayes(), ...)
)

addAlgorithm(
  reg = reg,
  name = "hodgesLehmann",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritHodgesLehmann(hodges.factor = 0.7), ...)
)

addAlgorithm(
  reg = reg,
  name = "ei",
  fun = function(job, data, instance, ...)
    mbo_optimization(job, data, instance, infill_crit = makeMBOInfillCritEI(), ...)
)


# prob.designs = algo.designs = list()
# prob.designs$rnorm = CJ(n = 100, mean =-1:1, sd = 1:5)
# prob.designs$rexp = data.table(n = 100, lambda = 1:5)
# algo.designs$average = data.table(method = c("mean", "median"))
# algo.designs$deviation = data.table()

ids = addExperiments(reg = reg, repls = 40)

# getJobPars(reg = tmp)
summarizeExperiments(reg = reg, by = c("problem", "algorithm"))

# profvis({
#   submitJobs(reg = reg)
# })

submitJobs(reg = reg)
waitForJobs(reg = reg)

reg$algorithms



#########################

resOG = ijoin(unwrap(getJobPars(ids, reg = reg)),
              reduceResultsDataTable(
                reg = reg,
                fun = function(x)
                  list(res = x)
              ))
res <- resOG |> as.data.frame()


# Step 1: Extract values from the nested list
res_values <- lapply(res$result, function(x) {
  best_y <- x$res[[1]]  # Extract res[[1]]
  running_time <-
    as.numeric(x$res[[2]])  # Extract res[[2]] and convert to numeric
  return(c(best_y, running_time))
})

# Step 2: Convert the list to a data frame
res_df <- do.call(rbind, res_values)  # Combine all rows
res_df
colnames(res_df) <- c("best.y", "running.time")  # Name the columns
res <- cbind(res, res_df)
res$result = NULL




plot <- ggplot(res, aes(x = algorithm, y = best.y)) +
  geom_boxplot(aes(fill = algorithm)) +
  facet_wrap( ~ problem, ncol = 2, scales = "free") +
  labs(title = "Best objective value (on y axis) found by respective algorithms on respective test
function",
       x = "algorithms",
       y = "best y") + theme_bw()

ggsave(path = "Plots",
       filename = "Benchmark1.png")

avg_running_time <- res %>%
  group_by(problem, algorithm) %>%
  summarise(avg_time = mean(running.time))

# Print the result
avg_running_time



# Add a rank column based on the best.y values
res <- res %>%
  group_by(problem) %>%
  mutate(rank = dense_rank(best.y))



avg_rank <- res %>%
  group_by(algorithm) %>%
  summarise(avg_rank = mean(rank))

avg_rank
avg_running_time
res

# res <- res %>%
#   left_join(avg_running_time, by = c("problem", "algorithm")) |>
#   left_join(avg_rank, by = c("problem", "algorithm"))
#
#
# # Print the final result
# res