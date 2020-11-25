RFGLS_predict_spatial <- function(RFGLS_out, coords, coords.0, Xtest, h = 1, verbose = FALSE){
  func_pred <- RFGLS_predict(RFGLS_out, Xtest, h = h, verbose = FALSE)
  func_pred_input <- RFGLS_predict(RFGLS_out, RFGLS_out$X, h = h, verbose = FALSE)
  rfgls_residual <- RFGLS_out$y - func_pred_input$predicted
  est_theta <- BRISC_estimation(coords, x = matrix(1,nrow(coords),1), y = rfgls_residual, verbose = FALSE)
  corr_pred <- BRISC_prediction(est_theta, coords.0, X.0 = matrix(1,nrow(coords.0),1), verbose = FALSE)

  prediction <- corr_pred$prediction +  func_pred$predicted
  RFGLS_prediction_spatial <- list()
  RFGLS_prediction_spatial$prediction <- prediction
  return(RFGLS_prediction_spatial)
}
