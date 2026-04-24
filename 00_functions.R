#This is convenient when you want to keep just the first pioneer of each tradition
keep_first_pioneer <- function(traditions) {
  lapply(traditions, function(tradition) {
    list(
      pioneers = tradition$pioneers[1],
      epigones = tradition$epigones
    )
  })
}

#This function prepares the data for Stan models
prepare_stan_data <- function(D, namlong, traditions) {
  dims <- ncol(D)
  mus  <- apply(D, 2, mean)
  sds  <- apply(D, 2, sd2)
  
  epi_vec <- mud_vec <- mup_vec <- sdd_vec <- sdr_vec <- sdp_vec <- NULL
  book_vec <- dim_vec <- rID_vec <- eID_vec <- NULL
  eID_cnt <- 0L
  
  sds_rm<-rep(NA, dims)
  sd_rm_avg <- NULL
  
  for (t in seq_along(traditions)) {
    tr    <- traditions[[t]]
    D_p   <- D[namlong %in% tr$pioneers, , drop=FALSE]
    mu_p  <- apply(D_p, 2, mean)
    sd_p  <- apply(D_p, 2, sd2)
    sd_rm <- sqrt(sum(sd_p^2) / length(sd_p))
    
    for (ep in tr$epigones) {
      D_e <- D[namlong == ep, , drop=FALSE]
      if (nrow(D_e)==0) next
      eID_cnt <- eID_cnt + 1L
      
      # 1) exact same flattening + names
      epi_block <- unlist(D_e)
      
      # 2) moments, all repeated in the same order as epi_block
      n_e       <- nrow(D_e)
      mud_block <- rep(mus, each = n_e)
      mup_block <- rep(mu_p, each = n_e)
      sdd_block <- rep(sds, each = n_e)
      sdr_block <- rep(rep(sd_rm,dims), each = n_e)
      sdp_block <- rep(sd_p, each = n_e)
      
      # 3) correct dim / id indices
      dim_block <- rep(seq_len(dims), each = n_e)
      book_block <- rep(seq_len(n_e), times = dims)
      rID_block <- rep(t, length(epi_block))
      eID_block <- rep(eID_cnt, length(epi_block))
      
      # 4) append
      epi_vec <- c(epi_vec, epi_block)
      mud_vec <- c(mud_vec, mud_block)
      mup_vec <- c(mup_vec, mup_block)
      sdd_vec <- c(sdd_vec, sdd_block)
      sdr_vec <- c(sdr_vec, sdr_block)
      sdp_vec <- c(sdp_vec, sdp_block)
      dim_vec <- c(dim_vec, dim_block)
      book_vec <- c(book_vec, book_block)
      rID_vec <- c(rID_vec, rID_block)
      eID_vec <- c(eID_vec, eID_block)
    }
    sds_rm <- rbind(sds_rm, sd_p)
    sd_rm_avg <- c(sd_rm_avg, sd_rm)
  }
  
  list(
  dA=list(
    epi = epi_vec,
    mud = mud_vec,
    mup = mup_vec,
    sdd = sdd_vec,
    sdr = sdr_vec,
    sdp = sdp_vec,
    dim = dim_vec,
    book = book_vec, #book within given epigone, this might not be used in models, but it helps to think about the dataset
    rID = rID_vec,
    eID = eID_vec,
    N   = length(epi_vec),
    Nr  = length(traditions),
    Ne  = max(eID_vec),
    Nd  = dims
  ),
  sds=sds,
  sds_rm=sds_rm[-1,],
  sd_rm_avg=sd_rm_avg)
}
