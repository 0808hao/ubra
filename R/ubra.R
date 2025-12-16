
## ==========================================================

UBRA_2D_V2A2_VERSION <- "UBRA 2D v2-A.2 (prod-grade): BIC+stability, SE-weighted, soft-MR (2025-12-12)"

## =========================
## 1) Internal helpers (non-exported)
## =========================

.safe_sd <- function(x){
  x <- x[is.finite(x)]
  if(length(x) < 2) return(0)
  stats::sd(x)
}

.quant <- function(x, p){
  x <- x[is.finite(x)]
  if(length(x) == 0) return(0)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, type = 7))
}

## --- CRAN-friendly seeding: preserve and restore RNG state
.with_seed <- function(seed, expr){
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if(has_seed){
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  force(expr)
}

## robust extraction of coef matrix from flexmix fit
.get_coef_matrix <- function(fm){
  cm <- try(flexmix::parameters(fm), silent = TRUE)
  if(!inherits(cm, "try-error") && is.matrix(cm) && ncol(cm) == fm@k) return(cm)
  
  cf <- try(stats::coef(fm), silent = TRUE)
  if(inherits(cf, "try-error") || !is.list(cf)) return(NULL)
  
  K <- fm@k
  mat <- matrix(NA_real_, nrow = 2, ncol = K)
  rownames(mat) <- c("(Intercept)", "beta_exposure")
  colnames(mat) <- paste0("Comp", seq_len(K))
  
  for(k in 1:min(K, length(cf))){
    v <- cf[[k]]
    if(is.null(v)) next
    if("(Intercept)" %in% names(v)) mat["(Intercept)", k] <- as.numeric(v["(Intercept)"])
    if("beta_exposure" %in% names(v)) mat["beta_exposure", k] <- as.numeric(v["beta_exposure"])
    if(is.na(mat["beta_exposure", k])){
      nm <- names(v)
      hit <- nm[grepl("beta_exposure", nm, fixed = TRUE)]
      if(length(hit) >= 1) mat["beta_exposure", k] <- as.numeric(v[hit[1]])
    }
  }
  mat
}

.coef_fallback_by_cluster <- function(beta_exposure0, beta_outcome0, cl_id, K){
  alpha <- rep(NA_real_, K)
  theta <- rep(NA_real_, K)
  
  for(k in 1:K){
    idx <- which(cl_id == k)
    if(length(idx) < 5){
      alpha[k] <- 0; theta[k] <- 0
      next
    }
    dd <- data.frame(y = beta_outcome0[idx], x = beta_exposure0[idx])
    fit <- try(stats::lm(y ~ x, data = dd), silent = TRUE)
    if(inherits(fit, "try-error")){
      alpha[k] <- 0; theta[k] <- 0
    } else {
      cc <- stats::coef(fit)
      alpha[k] <- ifelse("(Intercept)" %in% names(cc), as.numeric(cc["(Intercept)"]), 0)
      theta[k] <- ifelse("x" %in% names(cc), as.numeric(cc["x"]), 0)
    }
  }
  list(alpha = alpha, theta = theta)
}

## ARI (no extra package)
.ari <- function(x, y){
  x <- as.integer(as.factor(x))
  y <- as.integer(as.factor(y))
  n <- length(x)
  tab <- table(x, y)
  comb2 <- function(z) ifelse(z >= 2, z*(z-1)/2, 0)
  sum_ij <- sum(comb2(tab))
  sum_i  <- sum(comb2(rowSums(tab)))
  sum_j  <- sum(comb2(colSums(tab)))
  total  <- comb2(n)
  expected <- (sum_i * sum_j) / total
  maxval   <- 0.5 * (sum_i + sum_j)
  if(maxval - expected == 0) return(0)
  (sum_ij - expected) / (maxval - expected)
}

## Soft-IVW (posterior-weighted IVW)
.soft_ivw <- function(bx, by, se_by, post_w){
  ok <- is.finite(bx) & is.finite(by) & is.finite(se_by) & se_by > 0 & is.finite(post_w) & post_w > 0
  bx <- bx[ok]; by <- by[ok]; se <- se_by[ok]; w <- post_w[ok]
  if(length(bx) < 5) return(NULL)
  
  W <- w / (se^2)
  num <- sum(W * bx * by)
  den <- sum(W * bx^2)
  if(den <= 0) return(NULL)
  
  b <- num / den
  se_b <- sqrt(1 / den)
  z <- b / se_b
  p <- 2 * stats::pnorm(-abs(z))
  
  data.frame(
    method = "Soft-IVW",
    b = b, se = se_b, pval = p,
    lci = b - 1.96 * se_b,
    uci = b + 1.96 * se_b,
    nsnp = length(bx),
    stringsAsFactors = FALSE
  )
}

.fit_once_flexmix <- function(df, K, seed, use_weights = TRUE){
  fm <- .with_seed(seed, {
    try(
      flexmix::flexmix(
        beta_outcome ~ beta_exposure,
        data    = df,
        k       = K,
        model   = flexmix::FLXMRglm(family = "gaussian"),
        weights = if(use_weights && ".w" %in% names(df)) df$.w else NULL
      ),
      silent = TRUE
    )
  })
  if(inherits(fm, "try-error")) return(NULL)
  fm
}

.bind_rows_dt <- function(lst){
  if(length(lst) == 0) return(data.frame())
  if(!requireNamespace("data.table", quietly = TRUE)){
    ## base fallback (slower)
    return(do.call(rbind, lst))
  }
  as.data.frame(data.table::rbindlist(lst, fill = TRUE, use.names = TRUE))
}

.left_join_by_key <- function(x, y, key = "GRF2_label"){
  if(is.null(y) || nrow(y) == 0) return(x)
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  if(!(key %in% names(x)) || !(key %in% names(y))) return(x)
  merge(x, y, by = key, all.x = TRUE, sort = FALSE)
}

## =========================
## 2) Main: ubra_fit_2d_v2A2()
## =========================

#' Fit UBRA 2D v2-A.2 model (flexmix + stability; mclust fallback).
#'
#' @param beta_exposure Numeric vector of SNP-exposure effects.
#' @param beta_outcome Numeric vector of SNP-outcome effects.
#' @param se_outcome Optional numeric vector of outcome SEs (for SE weighting).
#' @param snp_id Optional SNP IDs; if NULL, auto-generated.
#' @param trait Trait label used in GRF labels.
#' @param G_range Candidate number of clusters (K). If NULL, adaptively chosen.
#' @param tail_p VaR quantile for |beta_outcome|.
#' @param tail_var_thr Tail-variance-share threshold for extreme-risk flag.
#' @param abs_beta_q_thr Quantile threshold on |mean_beta_out| for extreme-risk flag.
#' @param min_abs_gx Filter SNPs with |beta_exposure| < min_abs_gx.
#' @param min_cluster_n Mark clusters with n < min_cluster_n as degenerate.
#' @param use_se_weight If TRUE and se_outcome valid, use 1/se^2 weights.
#' @param n_start Number of random starts per K for BIC minimization.
#' @param n_resample Number of subsampling runs for ARI stability.
#' @param resample_frac Fraction of SNPs in each subsample.
#' @param stability_metric Currently only "ARI".
#' @param stability_weight Weight of stability term in selection score.
#' @param seed Base seed (RNG state is preserved and restored).
#' @param allow_fallback_mclust If TRUE, fallback to 2D mclust when flexmix fails.
#' @param verbose If TRUE, prints progress messages.
#' @return An object of class 'ubra_fit_2d_v2A2'.
#' @export
ubra_fit_2d_v2A2 <- function(
    beta_exposure,
    beta_outcome,
    se_outcome        = NULL,
    snp_id            = NULL,
    trait             = "Trait",
    G_range           = NULL,
    
    ## tail-risk
    tail_p            = 0.99,
    tail_var_thr      = 0.5,
    abs_beta_q_thr    = 0.75,
    
    ## filters
    min_abs_gx        = 0,
    min_cluster_n     = 10,
    
    ## weighting
    use_se_weight     = TRUE,
    
    ## stability
    n_start           = 5,
    n_resample        = 20,
    resample_frac     = 0.8,
    stability_metric  = c("ARI"),
    stability_weight  = 0.3,
    seed              = 2025,
    
    ## fallback
    allow_fallback_mclust = TRUE,
    verbose           = TRUE
){
  stability_metric <- match.arg(stability_metric)
  
  if(!requireNamespace("flexmix", quietly = TRUE)){
    stop("Package 'flexmix' is required. Please install it: install.packages('flexmix').")
  }
  if(!requireNamespace("mclust", quietly = TRUE)){
    stop("Package 'mclust' is required. Please install it: install.packages('mclust').")
  }
  
  bx <- as.numeric(beta_exposure)
  by <- as.numeric(beta_outcome)
  n  <- length(bx)
  if(length(by) != n) stop("beta_exposure and beta_outcome must have the same length.")
  
  if(is.null(snp_id)) snp_id <- paste0("SNP_", seq_len(n))
  if(length(snp_id) != n) stop("snp_id must have the same length as beta vectors.")
  snp_id <- as.character(snp_id)
  
  if(is.null(se_outcome)){
    se_outcome <- rep(NA_real_, n)
  } else {
    se_outcome <- as.numeric(se_outcome)
    if(length(se_outcome) != n) stop("se_outcome must have the same length as beta vectors.")
  }
  
  ## weights
  w <- rep(1, n)
  if(isTRUE(use_se_weight) && all(is.finite(se_outcome)) && all(se_outcome > 0)){
    w <- 1/(se_outcome^2)
  } else if(isTRUE(use_se_weight) && isTRUE(verbose)){
    message("[v2A.2] se_outcome not fully available; SE-weighting disabled (w=1).")
  }
  
  keep <- is.finite(bx) & is.finite(by) & is.finite(w) & w > 0
  if(!is.null(min_abs_gx) && min_abs_gx > 0) keep <- keep & (abs(bx) >= min_abs_gx)
  
  bx0 <- bx[keep]; by0 <- by[keep]; snp0 <- snp_id[keep]; se0 <- se_outcome[keep]; w0 <- w[keep]
  n0 <- length(bx0)
  if(n0 < 30) stop("Too few valid SNPs (<30) after filtering.")
  
  ## adaptive G_range
  if(is.null(G_range)){
    G_max <- if(n0 < 200) 3 else if(n0 < 500) 4 else if(n0 < 2000) 5 else 6
    G_range <- 1:G_max
  } else {
    G_range <- sort(unique(as.integer(G_range)))
    G_range <- G_range[G_range >= 1]
    if(length(G_range) == 0) stop("G_range must contain values >= 1.")
  }
  
  df <- data.frame(
    SNP          = snp0,
    beta_exposure= bx0,
    beta_outcome = by0,
    se_outcome   = se0,
    .w           = w0,
    stringsAsFactors = FALSE
  )
  
  ## -------------------------
  ## 2.1 select K by (BIC + stability)
  ## -------------------------
  if(isTRUE(verbose)) message("[v2A.2] Selecting K by BIC + stability...")
  
  cand <- data.frame(K = G_range, best_BIC = NA_real_, best_seed = NA_integer_, stringsAsFactors = FALSE)
  
  for(i in seq_along(G_range)){
    K <- G_range[i]
    bic_best  <- Inf
    seed_best <- NA_integer_
    
    for(s in 1:n_start){
      sd <- seed + 1000*K + s
      fm <- .fit_once_flexmix(df, K, sd, use_weights = TRUE)
      if(is.null(fm)) next
      b <- try(stats::BIC(fm), silent = TRUE)
      if(!inherits(b, "try-error") && is.finite(b) && b < bic_best){
        bic_best  <- b
        seed_best <- sd
      }
    }
    if(is.finite(bic_best)){
      cand$best_BIC[i]  <- bic_best
      cand$best_seed[i] <- seed_best
    }
  }
  
  ## ==========================================================
  ## FIX #1: fallback branch MUST return snp_table_2d/grf_table_2d
  ## ==========================================================
  if(all(!is.finite(cand$best_BIC))){
    if(!allow_fallback_mclust) stop("All flexmix fits failed for all K, and fallback is disabled.")
    if(isTRUE(verbose)) message("[v2A.2] flexmix failed for all K; fallback to mclust in 2D space.")
    
    mat <- cbind(scale(df$beta_exposure), scale(df$beta_outcome))
    mc <- mclust::Mclust(mat, G = G_range, verbose = FALSE)
    if(is.null(mc) || is.null(mc$classification)) stop("mclust fallback failed: empty classification.")
    
    cls  <- as.integer(mc$classification)
    Ksel <- as.integer(mc$G)
    if(!is.finite(Ksel) || Ksel < 1) stop("mclust fallback failed: invalid selected G.")
    
    post <- matrix(0, nrow = n0, ncol = Ksel)
    okc <- is.finite(cls) & cls >= 1 & cls <= Ksel
    post[cbind(which(okc), cls[okc])] <- 1
    colnames(post) <- paste0("post_G", seq_len(Ksel))
    
    GRF2_label <- paste0(trait, "_2DGRF", seq_len(Ksel))
    
    snp_table_2d <- data.frame(
      trait         = trait,
      SNP           = df$SNP,
      beta_exposure = df$beta_exposure,
      beta_outcome  = df$beta_outcome,
      se_outcome    = df$se_outcome,
      weights       = df$.w,
      cluster_2d    = cls,
      GRF2_label    = GRF2_label[pmax(1, pmin(Ksel, cls))],
      stringsAsFactors = FALSE
    )
    snp_table_2d <- cbind(snp_table_2d, post)
    
    var_beta_total <- stats::var(df$beta_outcome)
    grf_list <- vector("list", Ksel)
    
    for(k in 1:Ksel){
      idx <- which(cls == k)
      n_k <- length(idx)
      gx  <- df$beta_exposure[idx]
      byk <- df$beta_outcome[idx]
      
      mu_gx <- if(n_k > 0) mean(gx) else 0
      mu_by <- if(n_k > 0) mean(byk) else 0
      sd_gx <- .safe_sd(gx)
      sd_by <- .safe_sd(byk)
      
      cov_k0 <- if(n_k > 1) stats::cov(cbind(gx, byk))[1,2] else 0
      cor_k0 <- if(n_k > 2) suppressWarnings(stats::cor(gx, byk)) else NA_real_
      
      var_by_k   <- if(n_k > 1) stats::var(byk) else 0
      var_share_k <- if(is.finite(var_beta_total) && var_beta_total > 0) var_by_k/var_beta_total else 0
      
      abs_by <- abs(byk)
      VaR_k  <- .quant(abs_by, tail_p)
      tail_idx <- which(abs_by > VaR_k)
      if(length(tail_idx) > 0){
        CVaR_k <- mean(abs_by[tail_idx], na.rm = TRUE)
        tail_var_k  <- sum((byk[tail_idx] - mu_by)^2)
        total_var_k <- sum((byk - mu_by)^2)
        tail_var_share_k <- if(total_var_k > 0) tail_var_k/total_var_k else 0
      } else {
        CVaR_k <- 0
        tail_var_share_k <- 0
      }
      
      deg_flag <- (n_k < min_cluster_n)
      
      grf_list[[k]] <- data.frame(
        trait              = trait,
        GRF2_label         = GRF2_label[k],
        cluster            = k,
        n_snp              = n_k,
        weight_prior       = n_k / n0,
        alpha_intercept    = NA_real_,
        theta_causal       = NA_real_,
        sigma_resid        = NA_real_,
        degenerate_cluster = deg_flag,
        mean_beta_exp      = mu_gx,
        mean_beta_out      = mu_by,
        sd_beta_exp        = sd_gx,
        sd_beta_out        = sd_by,
        cov_exp_out        = cov_k0,
        cor_exp_out        = cor_k0,
        var_share_out      = var_share_k,
        var_share_out_pct  = var_share_k * 100,
        VaR_p              = VaR_k,
        CVaR_p             = CVaR_k,
        tail_var_share     = tail_var_share_k,
        stringsAsFactors   = FALSE
      )
    }
    
    grf_table_2d <- .bind_rows_dt(grf_list)
    
    abs_beta_mean <- abs(grf_table_2d$mean_beta_out)
    abs_beta_mean[grf_table_2d$degenerate_cluster] <- NA_real_
    thr_abs <- stats::quantile(abs_beta_mean, probs = abs_beta_q_thr, na.rm = TRUE)
    
    grf_table_2d$is_extreme <- ifelse(
      (!grf_table_2d$degenerate_cluster) &
        (grf_table_2d$tail_var_share >= tail_var_thr) &
        (abs(grf_table_2d$mean_beta_out) >= thr_abs),
      "Extreme-risk GRF", "Non-extreme GRF"
    )
    
    score_table <- data.frame(K = G_range, best_BIC = NA_real_, ARI = NA_real_, score = NA_real_)
    
    out <- list(
      trait         = trait,
      version       = UBRA_2D_V2A2_VERSION,
      df_used       = df,
      flexmix_fit   = NULL,
      K_selected    = Ksel,
      seed_selected = NA_integer_,
      BIC_table     = cand,
      stability     = data.frame(K = G_range, ARI = NA_real_),
      score_table   = score_table,
      G             = Ksel,
      cluster_2d    = cls,
      posterior     = post,
      grf_table_2d  = grf_table_2d,
      snp_table_2d  = snp_table_2d,
      tail_p        = tail_p,
      tail_var_thr  = tail_var_thr,
      abs_beta_q_thr= abs_beta_q_thr,
      min_abs_gx    = min_abs_gx,
      min_cluster_n = min_cluster_n,
      use_se_weight = use_se_weight,
      method        = "fallback_mclust_2d",
      call          = match.call()
    )
    class(out) <- "ubra_fit_2d_v2A2"
    return(out)
  }
  ## ===================== end FIX #1 =======================
  
  ## stability (ARI)
  stab <- data.frame(K = G_range, ARI = NA_real_, stringsAsFactors = FALSE)
  
  for(i in seq_along(G_range)){
    K <- G_range[i]
    if(!is.finite(cand$best_BIC[i])) next
    
    ref_seed <- cand$best_seed[i]
    fm_ref <- .fit_once_flexmix(df, K, ref_seed, use_weights = TRUE)
    if(is.null(fm_ref)) next
    
    ref_cls <- max.col(flexmix::posterior(fm_ref))
    
    ari_vec <- c()
    ## preserve RNG state during resampling loop
    ari_vec <- .with_seed(seed + 999 + K, {
      for(r in 1:n_resample){
        idx <- sample(seq_len(n0), size = max(30, floor(resample_frac * n0)), replace = FALSE)
        df_r <- df[idx, , drop = FALSE]
        fm_r <- .fit_once_flexmix(df_r, K, seed + 5000*K + r, use_weights = TRUE)
        if(is.null(fm_r)) next
        cls_r <- max.col(flexmix::posterior(fm_r))
        ari_vec <<- c(ari_vec, .ari(ref_cls[idx], cls_r))
      }
      ari_vec
    })
    if(length(ari_vec) > 0) stab$ARI[i] <- mean(ari_vec, na.rm = TRUE)
  }
  
  bic_ok <- cand$best_BIC
  bic_scaled <- (bic_ok - min(bic_ok, na.rm = TRUE)) / (stats::sd(bic_ok, na.rm = TRUE) + 1e-9)
  ari_ok <- stab$ARI
  ari_ok[!is.finite(ari_ok)] <- min(ari_ok, na.rm = TRUE)
  score <- -bic_scaled + stability_weight * ari_ok
  
  K_sel <- cand$K[which.max(score)]
  seed_sel <- cand$best_seed[cand$K == K_sel][1]
  
  if(isTRUE(verbose)){
    message("[v2A.2] Selected K = ", K_sel,
            " (best_BIC=", round(cand$best_BIC[cand$K == K_sel], 2),
            ", ARI=", round(stab$ARI[stab$K == K_sel], 3), ")")
  }
  
  fm_best <- .fit_once_flexmix(df, K_sel, seed_sel, use_weights = TRUE)
  if(is.null(fm_best)) stop("Internal error: flexmix fit failed at selected K.")
  
  post <- flexmix::posterior(fm_best)
  cl_id <- max.col(post)
  K <- fm_best@k
  
  ## component params
  coefs <- .get_coef_matrix(fm_best)
  alpha_k <- rep(NA_real_, K)
  theta_k <- rep(NA_real_, K)
  
  if(!is.null(coefs) && is.matrix(coefs)){
    rn <- rownames(coefs)
    if(!is.null(rn)){
      r_int <- which(rn %in% c("(Intercept)", "Intercept", "intercept"))
      r_bx  <- which(rn %in% c("beta_exposure", "beta.exposure", "x", "beta_exposure0"))
      if(length(r_int) >= 1) alpha_k <- as.numeric(coefs[r_int[1], ])
      if(length(r_bx)  >= 1) theta_k <- as.numeric(coefs[r_bx[1], ])
    }
  }
  
  if(any(!is.finite(alpha_k)) || any(!is.finite(theta_k))){
    fb <- .coef_fallback_by_cluster(df$beta_exposure, df$beta_outcome, cl_id, K)
    alpha_k[!is.finite(alpha_k)] <- fb$alpha[!is.finite(alpha_k)]
    theta_k[!is.finite(theta_k)] <- fb$theta[!is.finite(theta_k)]
  }
  
  w_k <- as.numeric(flexmix::prior(fm_best))
  sigma_k <- rep(NA_real_, K)
  for(k in 1:K){
    sk <- try(stats::sigma(fm_best@components[[k]]), silent = TRUE)
    if(!inherits(sk, "try-error")) sigma_k[k] <- as.numeric(sk)
  }
  
  GRF2_label <- paste0(trait, "_2DGRF", seq_len(K))
  var_beta_total <- stats::var(df$beta_outcome)
  
  deg_flag <- rep(FALSE, K)
  grf_list <- vector("list", K)
  
  for(k in 1:K){
    idx <- which(cl_id == k)
    n_k <- length(idx)
    if(n_k < min_cluster_n) deg_flag[k] <- TRUE
    
    gx  <- df$beta_exposure[idx]
    byk <- df$beta_outcome[idx]
    
    mu_gx <- if(n_k > 0) mean(gx) else 0
    mu_by <- if(n_k > 0) mean(byk) else 0
    sd_gx <- .safe_sd(gx)
    sd_by <- .safe_sd(byk)
    
    cov_k0 <- if(n_k > 1) stats::cov(cbind(gx, byk))[1,2] else 0
    cor_k0 <- if(n_k > 2) suppressWarnings(stats::cor(gx, byk)) else NA_real_
    
    var_by_k    <- if(n_k > 1) stats::var(byk) else 0
    var_share_k <- if(is.finite(var_beta_total) && var_beta_total > 0) var_by_k/var_beta_total else 0
    
    abs_by <- abs(byk)
    VaR_k  <- .quant(abs_by, tail_p)
    tail_idx <- which(abs_by > VaR_k)
    if(length(tail_idx) > 0){
      CVaR_k <- mean(abs_by[tail_idx], na.rm = TRUE)
      tail_var_k  <- sum((byk[tail_idx] - mu_by)^2)
      total_var_k <- sum((byk - mu_by)^2)
      tail_var_share_k <- if(total_var_k > 0) tail_var_k/total_var_k else 0
    } else {
      CVaR_k <- 0
      tail_var_share_k <- 0
    }
    
    slope_k <- if(deg_flag[k]) 0 else theta_k[k]
    alpha0  <- if(deg_flag[k]) 0 else alpha_k[k]
    sig0    <- if(deg_flag[k]) NA_real_ else sigma_k[k]
    
    grf_list[[k]] <- data.frame(
      trait              = trait,
      GRF2_label         = GRF2_label[k],
      cluster            = k,
      n_snp              = n_k,
      weight_prior       = w_k[k],
      alpha_intercept    = alpha0,
      theta_causal       = slope_k,
      sigma_resid        = sig0,
      degenerate_cluster = deg_flag[k],
      mean_beta_exp      = mu_gx,
      mean_beta_out      = mu_by,
      sd_beta_exp        = sd_gx,
      sd_beta_out        = sd_by,
      cov_exp_out        = cov_k0,
      cor_exp_out        = cor_k0,
      var_share_out      = var_share_k,
      var_share_out_pct  = var_share_k * 100,
      VaR_p              = VaR_k,
      CVaR_p             = CVaR_k,
      tail_var_share     = tail_var_share_k,
      stringsAsFactors   = FALSE
    )
  }
  
  grf_table_2d <- .bind_rows_dt(grf_list)
  
  abs_beta_mean <- abs(grf_table_2d$mean_beta_out)
  abs_beta_mean[grf_table_2d$degenerate_cluster] <- NA_real_
  thr_abs <- stats::quantile(abs_beta_mean, probs = abs_beta_q_thr, na.rm = TRUE)
  
  grf_table_2d$is_extreme <- ifelse(
    (!grf_table_2d$degenerate_cluster) &
      (grf_table_2d$tail_var_share >= tail_var_thr) &
      (abs(grf_table_2d$mean_beta_out) >= thr_abs),
    "Extreme-risk GRF", "Non-extreme GRF"
  )
  
  colnames(post) <- paste0("post_G", seq_len(K))
  
  snp_table_2d <- data.frame(
    trait         = trait,
    SNP           = df$SNP,
    beta_exposure = df$beta_exposure,
    beta_outcome  = df$beta_outcome,
    se_outcome    = df$se_outcome,
    weights       = df$.w,
    cluster_2d    = cl_id,
    GRF2_label    = GRF2_label[cl_id],
    stringsAsFactors = FALSE
  )
  snp_table_2d <- cbind(snp_table_2d, post)
  
  out <- list(
    trait         = trait,
    version       = UBRA_2D_V2A2_VERSION,
    df_used       = df,
    flexmix_fit   = fm_best,
    K_selected    = K_sel,
    seed_selected = seed_sel,
    BIC_table     = cand,
    stability     = stab,
    score_table   = data.frame(K = G_range, best_BIC = cand$best_BIC, ARI = stab$ARI, score = score),
    G             = K,
    cluster_2d    = cl_id,
    posterior     = post,
    grf_table_2d  = grf_table_2d,
    snp_table_2d  = snp_table_2d,
    tail_p        = tail_p,
    tail_var_thr  = tail_var_thr,
    abs_beta_q_thr= abs_beta_q_thr,
    min_abs_gx    = min_abs_gx,
    min_cluster_n = min_cluster_n,
    use_se_weight = use_se_weight,
    method        = "flexmix_2d",
    call          = match.call()
  )
  class(out) <- "ubra_fit_2d_v2A2"
  out
}

## =========================
## 3) S3 methods
## =========================

#' @export
print.ubra_fit_2d_v2A2 <- function(x, ...){
  cat("UBRA-2D v2-A.2 fit\n")
  cat("  trait      :", x$trait, "\n")
  cat("  version    :", x$version, "\n")
  cat("  method     :", ifelse(is.null(x$method), "NA", x$method), "\n")
  cat("  n_used     :", nrow(x$df_used), "\n")
  cat("  K_selected :", x$K_selected, " (G=", x$G, ")\n", sep = "")
  if(!is.null(x$stability) && "ARI" %in% names(x$stability)){
    ari_sel <- x$stability$ARI[x$stability$K == x$K_selected][1]
    cat("  ARI_mean   :", ifelse(is.finite(ari_sel), round(ari_sel, 3), NA), "\n")
  }
  invisible(x)
}

#' @export
summary.ubra_fit_2d_v2A2 <- function(object, ...){
  gt <- object$grf_table_2d
  if(is.null(gt) || nrow(gt) == 0){
    return(list(trait = object$trait, K = object$G, grf_table = NULL))
  }
  out <- gt[order(gt$is_extreme, -abs(gt$mean_beta_out), -gt$tail_var_share), ]
  rownames(out) <- NULL
  list(
    trait      = object$trait,
    version    = object$version,
    method     = object$method,
    K_selected = object$K_selected,
    G          = object$G,
    grf_table  = out
  )
}

## =========================
## 4) Plots (Suggests: ggplot2, ggrepel)
## =========================

#' @export
plot_ubra2d_scatter_v2A2 <- function(fit){
  stopifnot(inherits(fit, "ubra_fit_2d_v2A2"))
  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Package 'ggplot2' is required for plotting. Install it or set Suggests in your package.")
  }
  df <- fit$snp_table_2d
  ggplot2::ggplot(df, ggplot2::aes(x = beta_exposure, y = beta_outcome, color = GRF2_label)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.4, color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.4, color = "grey60") +
    ggplot2::geom_point(alpha = 0.75, size = 1.6) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    ggplot2::labs(
      x = "beta.exposure", y = "beta.outcome",
      title = paste0("UBRA-2D v2-A.2: ", fit$trait, " (BIC + stability)")
    )
}

#' @export
plot_ubra2d_tailrisk_v2A2 <- function(fit){
  stopifnot(inherits(fit, "ubra_fit_2d_v2A2"))
  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Package 'ggplot2' is required for plotting. Install it or set Suggests in your package.")
  }
  if(!requireNamespace("ggrepel", quietly = TRUE)){
    stop("Package 'ggrepel' is required for this plot. Install it or set Suggests in your package.")
  }
  gt <- fit$grf_table_2d
  ggplot2::ggplot(gt, ggplot2::aes(x = VaR_p, y = CVaR_p, label = GRF2_label, color = is_extreme)) +
    ggplot2::geom_point(size = 3.6, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 30, box.padding = 0.3, point.padding = 0.2) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = paste0("VaR(|beta.outcome|, p=", fit$tail_p, ")"),
      y = paste0("CVaR(|beta.outcome| | > VaR_", fit$tail_p, ")"),
      title = paste0("Tail risk by GRF (", fit$trait, ")")
    )
}

#' @export
plot_ubra2d_theta_alpha_v2A2 <- function(fit){
  stopifnot(inherits(fit, "ubra_fit_2d_v2A2"))
  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Package 'ggplot2' is required for plotting. Install it or set Suggests in your package.")
  }
  if(!requireNamespace("ggrepel", quietly = TRUE)){
    stop("Package 'ggrepel' is required for this plot. Install it or set Suggests in your package.")
  }
  gt <- fit$grf_table_2d
  ggplot2::ggplot(gt, ggplot2::aes(x = theta_causal, y = alpha_intercept, label = GRF2_label, size = n_snp, color = is_extreme)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey60") +
    ggplot2::geom_point(alpha = 0.85) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 30) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "theta (component slope): pathway causal effect",
      y = "alpha (intercept): mean offset / directional pleiotropy",
      title = paste0("theta vs alpha (", fit$trait, ")")
    )
}

#' @export
plot_ubra2d_K_selection_v2A2 <- function(fit){
  stopifnot(inherits(fit, "ubra_fit_2d_v2A2"))
  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Package 'ggplot2' is required for plotting. Install it or set Suggests in your package.")
  }
  sc <- fit$score_table
  sc$K <- as.factor(sc$K)
  ggplot2::ggplot(sc, ggplot2::aes(x = K, y = score, group = 1)) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_line() +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::labs(
      x = "Candidate K",
      y = "Selection score (-BIC_scaled + w*ARI)",
      title = paste0("K selection (selected K=", fit$K_selected, ")")
    )
}

## =========================
## 5) Runner (optional; uses data.table + ggplot2)
## =========================

#' Run UBRA-2D v2-A.2 from a harmonized data.frame and save tables/figures.
#' @param dat_harm data.frame with columns SNP, beta.exposure, beta.outcome, optionally se.outcome
#' @param trait Trait label
#' @param out_dir Output directory
#' @export
run_ubra2d_v2A2_from_dat_harm <- function(
    dat_harm, trait, out_dir,
    G_range = NULL,
    use_se_weight = TRUE,
    n_start = 5, n_resample = 20, resample_frac = 0.8,
    stability_weight = 0.3,
    min_cluster_n = 10,
    verbose = TRUE
){
  stopifnot(all(c("SNP", "beta.exposure", "beta.outcome") %in% names(dat_harm)))
  
  if(!requireNamespace("data.table", quietly = TRUE)){
    stop("Package 'data.table' is required for file writing in the runner. Install it or remove runner usage.")
  }
  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Package 'ggplot2' is required for figures in the runner. Install it or remove runner usage.")
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  se_out <- if("se.outcome" %in% names(dat_harm)) dat_harm$se.outcome else NULL
  
  fit <- ubra_fit_2d_v2A2(
    beta_exposure = dat_harm$beta.exposure,
    beta_outcome  = dat_harm$beta.outcome,
    se_outcome    = se_out,
    snp_id        = dat_harm$SNP,
    trait         = trait,
    G_range       = G_range,
    use_se_weight = use_se_weight,
    n_start       = n_start,
    n_resample    = n_resample,
    resample_frac = resample_frac,
    stability_weight = stability_weight,
    min_cluster_n = min_cluster_n,
    seed          = 2025,
    verbose       = verbose
  )
  
  ## tables
  data.table::fwrite(as.data.frame(fit$grf_table_2d), file.path(out_dir, paste0(trait, "_v2A2_GRF_table.tsv")), sep = "\t")
  data.table::fwrite(as.data.frame(fit$snp_table_2d), file.path(out_dir, paste0(trait, "_v2A2_SNP_table.tsv")), sep = "\t")
  data.table::fwrite(as.data.frame(fit$score_table),  file.path(out_dir, paste0(trait, "_v2A2_K_selection_table.tsv")), sep = "\t")
  
  ## figures (only if packages are present)
  p1 <- plot_ubra2d_scatter_v2A2(fit)
  ggplot2::ggsave(file.path(out_dir, paste0(trait, "_v2A2_scatter.png")), p1, width = 7.6, height = 5.6, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, paste0(trait, "_v2A2_scatter.pdf")), p1, width = 7.6, height = 5.6, bg = "white")
  
  p2 <- plot_ubra2d_tailrisk_v2A2(fit)
  ggplot2::ggsave(file.path(out_dir, paste0(trait, "_v2A2_tailrisk.png")), p2, width = 7.6, height = 5.6, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, paste0(trait, "_v2A2_tailrisk.pdf")), p2, width = 7.6, height = 5.6, bg = "white")
  
  p3 <- plot_ubra2d_theta_alpha_v2A2(fit)
  ggplot2::ggsave(file.path(out_dir, paste0(trait, "_v2A2_theta_alpha.png")), p3, width = 7.6, height = 5.6, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, paste0(trait, "_v2A2_theta_alpha.pdf")), p3, width = 7.6, height = 5.6, bg = "white")
  
  p4 <- plot_ubra2d_K_selection_v2A2(fit)
  ggplot2::ggsave(file.path(out_dir, paste0(trait, "_v2A2_K_selection.png")), p4, width = 7.2, height = 4.8, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, paste0(trait, "_v2A2_K_selection.pdf")), p4, width = 7.2, height = 4.8, bg = "white")
  
  ## soft MR per cluster
  snp_tbl <- fit$snp_table_2d
  K <- fit$G
  soft_mr_list <- list()
  
  for(k in 1:K){
    lab <- paste0(trait, "_2DGRF", k)
    wpost <- snp_tbl[[paste0("post_G", k)]]
    se_by <- snp_tbl$se_outcome
    if(!all(is.finite(se_by)) || any(se_by <= 0)) se_by <- rep(1, nrow(snp_tbl))
    
    res <- .soft_ivw(
      bx = snp_tbl$beta_exposure,
      by = snp_tbl$beta_outcome,
      se_by = se_by,
      post_w = wpost
    )
    if(is.null(res)) next
    res$GRF2_label <- lab
    res$K <- k
    soft_mr_list[[lab]] <- res
  }
  
  soft_mr_tbl <- if(length(soft_mr_list) > 0) .bind_rows_dt(soft_mr_list) else NULL
  if(!is.null(soft_mr_tbl)){
    data.table::fwrite(soft_mr_tbl, file.path(out_dir, paste0(trait, "_v2A2_SoftMR_by_GRF.tsv")), sep = "\t")
  }
  
  ## combined summary
  sum1 <- as.data.frame(fit$grf_table_2d)
  sum1$K_selected <- fit$K_selected
  ari_sel <- fit$stability$ARI[fit$stability$K == fit$K_selected][1]
  sum1$ARI_mean <- ari_sel
  sum1$version <- fit$version
  
  if(!is.null(soft_mr_tbl)){
    soft_keep <- soft_mr_tbl[, c("GRF2_label","b","se","pval","lci","uci","nsnp")]
    sum1 <- .left_join_by_key(sum1, soft_keep, key = "GRF2_label")
  }
  
  summary_file <- file.path(out_dir, paste0(trait, "_v2A2_ALL_SUMMARY.tsv"))
  data.table::fwrite(sum1, summary_file, sep = "\t")
  if(isTRUE(verbose)) message("[v2A.2] DONE. Summary saved: ", summary_file)
  
  fit
}

## ==========================================================
## 6) Optional "4D alias" (experimental; do NOT export initially)
## ==========================================================

UBRA_4D_VERSION <- UBRA_2D_V2A2_VERSION

ubra_fit_4d <- function(
    beta_exposure,
    beta_outcome,
    se_outcome        = NULL,
    snp_id            = NULL,
    trait             = "Trait",
    G_range           = NULL,
    
    tail_p            = 0.99,
    tail_var_thr      = 0.5,
    abs_beta_q_thr    = 0.75,
    
    min_abs_gx        = 0,
    min_cluster_n     = 10,
    
    use_se_weight     = TRUE,
    
    n_start           = 5,
    n_resample        = 20,
    resample_frac     = 0.8,
    stability_metric  = c("ARI"),
    stability_weight  = 0.3,
    seed              = 2025,
    
    allow_fallback_mclust = TRUE,
    verbose           = TRUE
){
  fit <- ubra_fit_2d_v2A2(
    beta_exposure   = beta_exposure,
    beta_outcome    = beta_outcome,
    se_outcome      = se_outcome,
    snp_id          = snp_id,
    trait           = trait,
    G_range         = G_range,
    tail_p          = tail_p,
    tail_var_thr    = tail_var_thr,
    abs_beta_q_thr  = abs_beta_q_thr,
    min_abs_gx      = min_abs_gx,
    min_cluster_n   = min_cluster_n,
    use_se_weight   = use_se_weight,
    n_start         = n_start,
    n_resample      = n_resample,
    resample_frac   = resample_frac,
    stability_metric= stability_metric,
    stability_weight= stability_weight,
    seed            = seed,
    allow_fallback_mclust = allow_fallback_mclust,
    verbose         = verbose
  )
  
  fit$method    <- "UBRA_4D_v2A2_alias"
  fit$version4d <- UBRA_4D_VERSION
  
  ## FIX #2: never call within() on NULL
  fit$grf_table <- fit$grf_table_2d
  
  if(is.null(fit$snp_table_2d)){
    df <- fit$df_used
    if(is.null(df)) stop("UBRA-4D alias: both snp_table_2d and df_used are NULL.")
    K <- if(!is.null(fit$G)) fit$G else 1
    GRF2_label <- paste0(trait, "_2DGRF", seq_len(K))
    cls <- if(!is.null(fit$cluster_2d)) fit$cluster_2d else rep(1L, nrow(df))
    fit$snp_table_2d <- data.frame(
      trait         = trait,
      SNP           = df$SNP,
      beta_exposure = df$beta_exposure,
      beta_outcome  = df$beta_outcome,
      se_outcome    = df$se_outcome,
      weights       = df$.w,
      cluster_2d    = cls,
      GRF2_label    = GRF2_label[pmax(1, pmin(K, cls))],
      stringsAsFactors = FALSE
    )
  }
  
  fit$snp_table <- within(fit$snp_table_2d, {
    cluster_4d <- cluster_2d
    GRF4_label <- GRF2_label
  })
  
  class(fit) <- unique(c("ubra_fit_4d_v2A2", class(fit)))
  fit
}

run_ubra4d_v2A2_from_dat_harm <- function(
    dat_harm, trait, out_dir,
    G_range = NULL,
    use_se_weight = TRUE,
    n_start = 5, n_resample = 20, resample_frac = 0.8,
    stability_weight = 0.3,
    min_cluster_n = 10,
    verbose = TRUE
){
  fit <- run_ubra2d_v2A2_from_dat_harm(
    dat_harm         = dat_harm,
    trait            = trait,
    out_dir          = out_dir,
    G_range          = G_range,
    use_se_weight    = use_se_weight,
    n_start          = n_start,
    n_resample       = n_resample,
    resample_frac    = resample_frac,
    stability_weight = stability_weight,
    min_cluster_n    = min_cluster_n,
    verbose          = verbose
  )
  
  fit$method    <- "UBRA_4D_v2A2_alias"
  fit$version4d <- UBRA_4D_VERSION
  fit$grf_table <- fit$grf_table_2d
  
  if(is.null(fit$snp_table_2d)){
    df <- fit$df_used
    if(is.null(df)) stop("UBRA-4D alias: both snp_table_2d and df_used are NULL.")
    K <- if(!is.null(fit$G)) fit$G else 1
    GRF2_label <- paste0(trait, "_2DGRF", seq_len(K))
    cls <- if(!is.null(fit$cluster_2d)) fit$cluster_2d else rep(1L, nrow(df))
    fit$snp_table_2d <- data.frame(
      trait         = trait,
      SNP           = df$SNP,
      beta_exposure = df$beta_exposure,
      beta_outcome  = df$beta_outcome,
      se_outcome    = df$se_outcome,
      weights       = df$.w,
      cluster_2d    = cls,
      GRF2_label    = GRF2_label[pmax(1, pmin(K, cls))],
      stringsAsFactors = FALSE
    )
  }
  
  fit$snp_table <- within(fit$snp_table_2d, {
    cluster_4d <- cluster_2d
    GRF4_label <- GRF2_label
  })
  
  class(fit) <- unique(c("ubra_fit_4d_v2A2", class(fit)))
  fit
}
