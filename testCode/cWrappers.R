
#' @useDynLib MRTransmissionModel calc_phi_and_update_tran_

calc_phi_and_update_tran <- function(waifw,
                                     state,
                                     s_indsR,
                                     i_indsR,
                                     exponentR,
                                     denomR,
                                     tran_matrix){

  .Call(calc_phi_and_update_tran_,
        waifw,
        state,
        s_indsR,
        i_indsR,
        exponentR,
        denomR,
        tran_matrix)

}


#' @useDynLib MRTransmissionModel do_ID_transition_SIR_stochastic_moves_cl_

do_ID_transition_SIR_stochastic_moves_cl <- function(curstateR,
                                                     transmatrixR){

  .Call(do_ID_transition_SIR_stochastic_moves_cl_,
        curstateR,transmatrixR)

}




#' @useDynLib MRTransmissionModel calc_phi_and_update_tran
calc_phi_and_update_tran_wrapper <- function(waifw,
                                             state,
                                             s_indsR,
                                             i_indsR,
                                             exponentR,
                                             denomR,
                                             tran_matrix){

  result <- .Call(calc_phi_and_update_tran,
                  waifw,
                  state,
                  s_indsR,
                  i_indsR,
                  exponentR,
                  denomR,
                  tran_matrix)
  return(result)

}




