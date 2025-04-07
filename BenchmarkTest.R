library(mlrMBO)         # Load the mlrMBO package for model-based optimisation
library(batchtools)     # Load batchtools for parallel execution and experiment management
library(ggplot2)        # Load ggplot2 for visualisation
library(data.table)     # Load data.table for efficient data handling
library(tidyr)          # Load tidyr for data tidying
library(dplyr)          # Load dplyr for data manipulation
library(rgenoud)        # Load rgenoud for genetic optimisation

# ---- Set global seed for reproducibility ----
set.seed(1234)  # Critical for all stochastic steps

# Function to perform model-based optimisation (MBO)
mbo_optimisation <-
  function(job, data, instance, infill_crit, ...) {
    set.seed(job$job.id)
    # Generate the initial design using Latin Hypercube Sampling (LHS)
    design = generateDesign(
      5 * getNumberOfParameters(instance),
      getParamSet(instance),
      fun = lhs::maximinLHS,
    )
    
    # Define the surrogate model using Kriging with a Matern 3/2 covariance function
    surrogate = makeLearner(
      "regr.km",
      predict.type = "se",
      covtype = "matern3_2"
    )
    
    # Create MBO control settings
    ctrl = makeMBOControl()
    ctrl = setMBOControlTermination(ctrl, iters = 100)  # Terminate after 100 iterations
    ctrl = setMBOControlInfill(ctrl, crit = infill_crit)  # Set the infill criterion
    # Set a seed for the MBO optimization process (unique per job)
    #ctrl = setMBOControlSeed(ctrl, seed = job$job.id)
    
    # Run the MBO optimisation process
    result = mbo(instance, design = design, control = ctrl)
    run.time = result$final.opt.state$time.used  # Capture runtime
    best.y = result$y  # Capture the best function value found
    
    return(list(best.y, run.time))  # Return best value and runtime
  }

# Create an experiment registry
reg = makeExperimentRegistry(file.dir = NA,
                             seed = 1234,
                             # Seed for batchtools
                             make.default = FALSE)

# Add benchmark optimisation problems to the registry
addProblem(name = "2DAckley",
           data = makeAckleyFunction(2),
           reg = reg)
addProblem(name = "2DDeflected",
           data = makeDeflectedCorrugatedSpringFunction(2),
           reg = reg)
addProblem(name = "2DRosenbrock",
           data = makeRosenbrockFunction(2),
           reg = reg)
addProblem(name = "2DSchwefel",
           data = makeSchwefelFunction(2),
           reg = reg)

addProblem(name = "4DAckley",
           data = makeAckleyFunction(4),
           reg = reg)
addProblem(name = "4DDeflected",
           data = makeDeflectedCorrugatedSpringFunction(4),
           reg = reg)
addProblem(name = "4DRosenbrock",
           data = makeRosenbrockFunction(4),
           reg = reg)
addProblem(name = "4DSchwefel",
           data = makeSchwefelFunction(4),
           reg = reg)

addProblem(name = "6DAckley",
           data = makeAckleyFunction(6),
           reg = reg)
addProblem(name = "6DDeflected",
           data = makeDeflectedCorrugatedSpringFunction(6),
           reg = reg)
addProblem(name = "6DRosenbrock",
           data = makeRosenbrockFunction(6),
           reg = reg)
addProblem(name = "6DSchwefel",
           data = makeSchwefelFunction(6),
           reg = reg)

# Add different optimisation algorithms to the registry
addAlgorithm(
  reg = reg,
  name = "ei",
  fun = function(job, data, instance, ...)
    mbo_optimisation(job, data, instance, infill_crit = makeMBOInfillCritEI(), ...)
)

addAlgorithm(
  reg = reg,
  name = "hodgesLehmann",
  fun = function(job, data, instance, ...)
    mbo_optimisation(job, data, instance, infill_crit = makeMBOInfillCritHodgesLehmann(), ...)
)

addAlgorithm(
  reg = reg,
  name = "hurwicz",
  fun = function(job, data, instance, ...)
    mbo_optimisation(job, data, instance, infill_crit = makeMBOInfillCritHurwicz(), ...)
)

addAlgorithm(
  reg = reg,
  name = "maximax",
  fun = function(job, data, instance, ...)
    mbo_optimisation(job, data, instance, infill_crit = makeMBOInfillCritMaximax(), ...)
)

addAlgorithm(
  reg = reg,
  name = "maximin",
  fun = function(job, data, instance, ...)
    mbo_optimisation(job, data, instance, infill_crit = makeMBOInfillCritMaximin(), ...)
)

# Create experiment instances and execute them
ids = addExperiments(reg = reg, repls = 10)
summarizeExperiments(reg = reg, by = c("problem", "algorithm"))
submitJobs(reg = reg)
waitForJobs(reg = reg)


## Working with the results
# Retrieve experiment results
resOG = ijoin(unwrap(getJobPars(ids, reg = reg)),
              reduceResultsDataTable(
                reg = reg,
                fun = function(x)
                  list(res = x)
              ))
res <- resOG |> as.data.frame()

# Process experiment results
res_values <- lapply(res$result, function(x) {
  best_y <- x$res[[1]]  # Extract best function value found
  running_time <- as.numeric(x$res[[2]])  # Extract runtime
  return(c(best_y, running_time))
})

# Convert the results into a data frame
res_df <- do.call(rbind, res_values)
colnames(res_df) <- c("best.y", "running.time")  # Rename columns
res <- cbind(res, res_df)
res$result = NULL  # Remove unnecessary column

# Create boxplot of best function values found by algorithms
plot <- ggplot(res, aes(x = algorithm, y = best.y)) +
  geom_boxplot(aes(fill = algorithm)) +
  facet_wrap(~ problem, ncol = 2, scales = "free") +
  labs(title = "Best objective value found by respective algorithms",
       x = "Algorithms",
       y = "Best y") + theme_bw()

# Ensure the Plots directory exists
if (!dir.exists("Plots")) {
  dir.create("Plots")
}

ggsave(path = "Plots", filename = "BenchmarkOutcome.png")  # Save the plot

# Compute average running time per problem and algorithm
avg_running_time <- res %>%
  group_by(problem, algorithm) %>%
  summarise(avg_time = mean(running.time))

# Compute ranking based on best objective values
res <- res %>%
  group_by(problem) %>%
  mutate(rank = dense_rank(best.y))

# Compute average ranking per algorithm
res_summary <- res %>%
  group_by(algorithm) %>%
  summarise(avg_rank = mean(rank),
            avg_time = mean(running.time)) |>
  arrange(avg_rank, avg_time)

# Display results
avg_running_time |> arrange(avg_time) # Sort in ascending order
res_summary
plot
res |> arrange(rank, running.time) |> print(n = 600) 