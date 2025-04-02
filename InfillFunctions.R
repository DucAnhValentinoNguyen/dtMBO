# This file contains all the implemented infill functions
# Note: Without loss of generality, objective functions are minimised throughout
# this work.
# For the benchmark tets we will only use the first six infill functions.

library(mlrMBO)  # Load the mlrMBO library for Bayesian Optimisation

# Infill function 1
# --- Expected Improvement (EI) Criterion ---
# The EI criterion balances exploration and exploitation in Bayesian Optimisation.
makeMBOInfillCritEI = function(se.threshold = 1e-6) {
  # Ensure se.threshold is a positive number
  assertNumber(se.threshold, lower = 1e-20)
  force(se.threshold)  # Ensure argument is evaluated
  
  # Create and return the EI criterion function
  makeMBOInfillCrit(
    fun = function(points,
                   models,
                   control,
                   par.set,
                   designs,
                   iter,
                   progress,
                   attributes = FALSE) {
      model = models[[1L]]  # Extract first surrogate model
      design = designs[[1]]  # Extract design points
      maximize.mult = if (control$minimize)
        1
      else-1  # Adjust for minimisation or maximisation
      
      # Extract objective values from design
      y = maximize.mult * design[, control$y.name]
      assertNumeric(y, any.missing = FALSE)
      
      # Predict new values
      p = predict(model, newdata = points)$data
      p.mu = maximize.mult * p$response  # Mean prediction
      p.se = p$se  # Standard error
      
      # Compute EI components
      y.min = min(y)  # Best observed objective value
      d = y.min - p.mu  # Difference between best and predicted value
      xcr = d / p.se  # Standardised improvement
      xcr.prob = pnorm(xcr)  # Probability of improvement
      xcr.dens = dnorm(xcr)  # Density of normal distribution
      ei = d * xcr.prob + p.se * xcr.dens  # EI formula
      
      # Apply thresholding to avoid numerical issues
      res = ifelse(p.se < se.threshold, 0,-ei)
      return(res)
    },
    name = "Expected improvement",
    id = "ei",
    components = c("se", "mean"),
    params = list(se.threshold = se.threshold),
    opt.direction = "maximize",
    requires.se = TRUE
  )
}

# Infill function 2
# --- Hodges-Lehmann Criterion ---
# A robust criterion that combines the idea of Bayesian criterion
# (uncertainty type I) and Maximin criterion (uncertainty type II).
makeMBOInfillCritHodgesLehmann = function(hodges.factor = 0.5,
                                          cb.lambda = 1) {
  assertNumber(hodges.factor,
               lower = 0,
               upper = 1,
               null.ok = TRUE)
  assertNumber(cb.lambda,
               lower = 0,
               upper = 1,
               null.ok = TRUE)
  
  makeMBOInfillCrit(
    fun = function(points,
                   models,
                   control,
                   par.set,
                   designs,
                   iter,
                   progress,
                   attributes = FALSE) {
      model = models[[1L]]
      p = predict(model, newdata = points)$data
      
      # Compute the modified response
      if (control$minimize) {
        res = hodges.factor*p$response + (1 - hodges.factor) *(p$response + cb.lambda * p$se)
      } else {
        res = - hodges.factor*p$response - (1 - hodges.factor) *(p$response - cb.lambda * p$se)
      }
      return(res)
    },
    name = "Hodges Lehmann with confidence bound",
    id = "hodgeslehmann",
    components = c("se", "mean", "lambda"),
    params = list(hodges.factor = hodges.factor, cb.lambda = cb.lambda),
    requires.se = TRUE,
    opt.direction = "objective"
  )
}

# Infill function 3
# --- Hurwicz Criterion ---
# Inspired by a decision criterion based on optimism-pessimism balance.
makeMBOInfillCritHurwicz = function(hurwicz.lambda = 0.5,
                                    cb.lambda = 1) {
  assertNumber(
    hurwicz.lambda,
    lower = 0,
    upper = 1,
    null.ok = TRUE
  )
  assertNumber(cb.lambda,
               lower = 0,
               upper = 1,
               null.ok = TRUE)
  
  makeMBOInfillCrit(
    fun = function(points,
                   models,
                   control,
                   par.set,
                   designs,
                   iter,
                   progress,
                   attributes = FALSE) {
      model = models[[1L]]
      p = predict(model, newdata = points)$data
      
      # Compute Hurwicz adjustment
      if (control$minimize) {
        res = p$response + (1 - 2 * hurwicz.lambda) * cb.lambda * p$se
      } else {
        res = -p$response - (2 * hurwicz.lambda - 1) * cb.lambda * p$se
      }
      return(res)
    },
    name = "Hurwicz confidence bound",
    id = "hurwicz",
    components = c("se", "mean", "lambda"),
    params = list(hurwicz.lambda = hurwicz.lambda, cb.lambda = cb.lambda),
    requires.se = TRUE,
    opt.direction = "objective"
  )
}

# Infill function 4
# --- Maximax Criterion ---
# An infill function inspired by a highly optimistic criterion that seeks the maximum possible gain.
makeMBOInfillCritMaximax = function(cb.lambda = 1) {
  assertNumber(cb.lambda, lower = 0, null.ok = TRUE)
  makeMBOInfillCrit(
    fun = function(points,
                   models,
                   control,
                   par.set,
                   designs,
                   iter,
                   progress,
                   attributes = FALSE) {
      model = models[[1L]]
      p = predict(model, newdata = points)$data
      res = ifelse(control$minimize, 1,-1) * p$response - cb.lambda * p$se
      return(res)
    },
    name = "Maximax with confidence bound",
    id = "maximax",
    components = c("se", "mean", "lambda"),
    params = list(cb.lambda = cb.lambda),
    requires.se = TRUE,
    opt.direction = "objective"
  )
}

# Infill function 5
# --- Maximin Criterion ---
# inspired by the maximin criterion that focuses on the worst-case scenario and
# tries to optimise the outcome nevertheless
makeMBOInfillCritMaximin = function(cb.lambda = 1) {
  assertNumber(cb.lambda, lower = 0, null.ok = TRUE)
  makeMBOInfillCrit(
    fun = function(points,
                   models,
                   control,
                   par.set,
                   designs,
                   iter,
                   progress,
                   attributes = FALSE) {
      model = models[[1L]]
      p = predict(model, newdata = points)$data
      
      # Compute Maximin criterion
      if (control$minimize) {
        res = p$response + cb.lambda * p$se
      } else {
        res = -p$response + cb.lambda * p$se
      }
      return(res)
    },
    name = "Maximin with confidence bound",
    id = "maximin",
    components = c("se", "mean", "lambda"),
    params = list(cb.lambda = cb.lambda),
    requires.se = TRUE,
    opt.direction = "objective"
  )
}

# Infill function 6
# --- Bayes Criterion ---
# Inspired by Bayesian criterion, this infill function is similar to the pure
# mean infill function, which ignores uncertainty.
makeMBOInfillCritBayes = function() {
  makeMBOInfillCrit(
    fun = function(points,
                   models,
                   control,
                   par.set,
                   designs,
                   iter,
                   progress,
                   attributes = FALSE) {
      ifelse(control$minimize, 1,-1) * predict(models[[1L]], newdata = points)$data$response
    },
    name = "Bayes Criterion",
    id = "bayes",
    opt.direction = "objective"
  )
}