RFGLS_predict_spatial <- function(RFGLS_out, coords.0, Xtest, h = 1, verbose = FALSE){

  if(missing(coords.0)){stop("error: coords.0 must be specified\n")}
  if(!any(is.data.frame(coords.0), is.matrix(coords.0))){stop("error: coords.0 must be a data.frame or matrix\n")}
  if(!ncol(coords.0) == 2){stop("error: coords.0 must have two columns\n")}

  coords <- RFGLS_out$coords
  func_pred <- RFGLS_predict(RFGLS_out, Xtest, h = h, verbose = verbose)
  func_pred_input <- RFGLS_predict(RFGLS_out, RFGLS_out$X, h = h, verbose = verbose)
  rfgls_residual <- RFGLS_out$y - func_pred_input$predicted
  est_theta <- BRISC_estimation(coords, x = matrix(1,nrow(coords),1), y = rfgls_residual, verbose = verbose)
  corr_pred <- BRISC_prediction(est_theta, coords.0, X.0 = matrix(1,nrow(coords.0),1), verbose = verbose)

  prediction <- corr_pred$prediction +  func_pred$predicted
  RFGLS_prediction_spatial <- list()
  RFGLS_prediction_spatial$prediction <- prediction
  return(RFGLS_prediction_spatial)
}
