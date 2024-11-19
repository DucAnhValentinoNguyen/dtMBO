library(mlrMBO)


#' @export
#' @rdname infillcrits
#' maximin
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
      if (control$minimize) {
        res = p$response + cb.lambda * p$se
      } else {
        res = -p$response + cb.lambda * p$se
      }
      return(res)
    },
    name = "Maximin confidence bound",
    id = "maximin",
    components = c("se", "mean", "lambda"),
    params = list(cb.lambda = cb.lambda),
    requires.se = TRUE,
    opt.direction = "objective"
  )
}


#' @export
#' @rdname infillcrits
#' maximax
#'
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
      maximize.mult = if (control$minimize)
        1
      else
        - 1
      p = predict(model, newdata = points)$data
      res = maximize.mult * p$response - cb.lambda * p$se
      return(res)
    },
    name = "Maximax",
    id = "maximax",
    components = c("se", "mean", "lambda"),
    params = list(cb.lambda = cb.lambda),
    opt.direction = "objective",
    requires.se = TRUE
  )
}

#' @export
#' @rdname infillcrits
#' Hurwicz
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
      if (control$minimize) {
        res = hurwicz.lambda * (p$response + cb.lambda * p$se)
      } else {
        res = -p$response + cb.lambda * p$se
      }
      return(res)
    },
    name = "Maximin confidence bound",
    id = "maximin",
    components = c("se", "mean", "lambda"),
    params = list(hurwicz.lambda = hurwicz.lambda, cb.lambda = cb.lambda),
    requires.se = TRUE,
    opt.direction = "objective"
  )
}


#' @export
#' @rdname infillcrits
#' Bayes
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
      ifelse(control$minimize, 1, -1) * predict(models[[1L]], newdata = points)$data$response
    },
    name = "Bayes Criterion",
    id = "bayes",
    opt.direction = "objective"
  )
}

#' @export
#' @rdname infillcrits
#' Hodges Lehmann Criterion
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
      if (control$minimize) {
        res = hodges.factor * p$response + (1 - hodges.factor) * (p$response + cb.lambda * p$se)
      } else {
        res = -hodges.factor * p$response - (1 - hodges.factor) * (p$response - cb.lambda * p$se)
      }
      return(res)
    },
    name = "Maximin confidence bound",
    id = "maximin",
    components = c("se", "mean", "lambda"),
    params = list(hodges.factor = hodges.factor, cb.lambda = cb.lambda),
    requires.se = TRUE,
    opt.direction = "objective"
  )
}

#' @export
#' @rdname infillcrits
makeMBOInfillCritEI = function(se.threshold = 1e-6) {
  assertNumber(se.threshold, lower = 1e-20)
  force(se.threshold)
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
      design = designs[[1]]
      maximize.mult = if (control$minimize)
        1
      else-1
      assertString(control$y.name)
      y = maximize.mult * design[, control$y.name]
      assertNumeric(y, any.missing = FALSE)
      p = predict(model, newdata = points)$data
      p.mu = maximize.mult * p$response
      p.se = p$se
      y.min = min(y)
      d = y.min - p.mu
      xcr = d / p.se
      xcr.prob = pnorm(xcr)
      xcr.dens = dnorm(xcr)
      ei = d * xcr.prob + p.se * xcr.dens
      res = ifelse(p.se < se.threshold, 0,-ei)
      # if (attributes) {
      #   res = setAttribute(res, "crit.components", data.frame(se = p$se, mean = p$response))
      # }
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