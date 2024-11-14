# 1. Maximax (optimistic)
# 2. Maximin (pessimistic)
# 3. Minimax regret

# 4. Equally likely (Laplace)
# 5. Criterion of realism (Hurwicz)
# 6. bayes decision theory

########################### CB based #########################
library(mlrMBO)

####################################### 1 maximin
# max of LCB
#' @export
#' @rdname infillcrits
makeMBOInfillCritMaximinCB = function() {
  # assertNumber(cb.lambda, lower = 0, null.ok = TRUE)
  # force(cb.lambda)
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
      res = p$response + p$se
      # if (attributes) {
      #   res = setAttribute(
      #     res,
      #     "crit.components",
      #     data.frame(
      #       se = p$se,
      #       mean = p$response,
      #       lambda = cb.lambda
      #     )
      #   )
      # }
      # return(res)
    },
    name = "Min upper confidence bound",
    id = "minucb",
    components = c("se", "mean"),
    # params = list(cb.lambda = cb.lambda),
    opt.direction = "objective",
    requires.se = TRUE
  )
}

##################### 2 minimax
# min of UCB
#' @export
#' @rdname infillcrits
makeMBOInfillCritMinimaxCB = function() {
  # assertNumber(cb.lambda, lower = 0, null.ok = TRUE)
  # force(cb.lambda)
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
      res =  p$response
      # if (attributes) {
      #   res = setAttribute(
      #     res,
      #     "crit.components",
      #     data.frame(
      #       se = p$se,
      #       mean = p$response,
      #       lambda = cb.lambda
      #     )
      #   )
      # }
      # return(res)
    },
    name = "Minimax Confidence Bound",
    id = "minimaxCB",
    components = c("mean"),
    # params = list(cb.lambda = cb.lambda),
    opt.direction = "objective",
    requires.se = TRUE
  )
}


##################################### maximax = UCB

##################################### Hurwicz

#' @export
#' @rdname infillcrits
makeMBOInfillCritHurwicz = function(hurwicz.alpha = 0.5) {
  assertNumber(hurwicz.alpha, lower = 0, upper = 1)
  force(hurwicz.alpha)
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
      # maximize.mult = if (control$minimize)
      #   - 1
      # else
      #   1
      p = predict(model, newdata = points)$data
      # Hurwicz criterion: a weighted combination of the mean and the pessimism factor (standard error)
      # hurwicz.alpha weights optimism (mean) vs pessimism (uncertainty)
      res = hurwicz.alpha * (p$response + p$se) + (1 - hurwicz.alpha) * (p$response - p$se)
      # if (attributes) {
      #   res = setAttribute(
      #     res,
      #     "crit.components",
      #     data.frame(
      #       se = p$se,
      #       mean = p$response,
      #       hurwicz.alpha = hurwicz.alpha
      #     )
      #   )
      # }
      # return(res)
    },
    name = "Hurwicz criterion",
    id = "hurwicz",
    components = c("se", "mean", "hurwicz.alpha"),
    params = list(hurwicz.alpha = hurwicz.alpha),
    opt.direction = "objective",
    requires.se = TRUE
  )
}


####################################

#' @export
#' @rdname infillcrits
makeMBOInfillCritBayesDecision = function(loss.fn = "squared") {
  assertChoice(loss.fn, choices = c("squared"))
  force(loss.fn)
  
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
      # maximize.mult = if (control$minimize)
      #   1
      # else
      #   - 1
      p = predict(model, newdata = points)$data
      # Bayes Decision Theory: minimize expected risk (Bayes risk)
      # For squared loss, the risk is mean^2 + variance
      if (loss.fn == "squared") {
        # Here we assume the loss is a squared error: mean^2 + variance (p$se^2)
        res = maximize.mult * (p$response ^ 2 + p$se ^ 2)
      }
      
      # if (attributes) {
      #   res = setAttribute(res,
      #                      "crit.components",
      #                      data.frame(
      #                        se = p$se,
      #                        mean = p$response,
      #                        bayes.risk = res
      #                      ))
      # }
      # return(res)
    },
    name = "Bayes Decision Criterion",
    id = "bayes_decision",
    components = c("se", "mean", "bayes.risk"),
    params = list(loss.fn = loss.fn),
    opt.direction = "objective",
    requires.se = TRUE
  )
}


#' @export
#' @rdname infillcrits
makeMBOInfillCritHodgesLehmann = function(hodges.factor = 0.5) {
  assertNumber(hodges.factor, lower = 0, upper = 1)
  force(hodges.factor)
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
      # maximize.mult = if (control$minimize)
      #   1
      # else-1
      p = predict(model, newdata = points)$data
      
      # Hodges-Lehmann estimator: A robust estimate of central tendency
      # Approximate it as a weighted combination of mean and uncertainty (standard error)
      # 
      # TODO: how to deal with finding minimum? take only minus by the part of maximin or 
      # take minus of the whole func?
      #
      res = hodges.factor * p$response + (1 - hodges.factor) * (p$response + p$se)
      # if (attributes) {
      #   res = setAttribute(
      #     res,
      #     "crit.components",
      #     data.frame(
      #       se = p$se,
      #       mean = p$response,
      #       hodges.factor = hodges.factor
      #     )
      #   )
      # }
      # return(res)
    },
    name = "Hodges-Lehmann Criterion",
    id = "hodges_lehmann",
    components = c("se", "mean", "hodges.factor"),
    params = list(hodges.factor = hodges.factor),
    opt.direction = "objective",
    requires.se = TRUE
  )
}
