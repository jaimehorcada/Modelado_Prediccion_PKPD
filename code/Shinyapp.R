# app.R
# =========================================================
# Shiny QSP PF-06671008 (Humanos + Ratón)
# =========================================================

# Si no se tienen intalado, intalar las siguientes librerias:
# install.packages("shiny")
# install.packages("deSolve")
# install.packages("ggplot2")
# install.packages("pracma")

library(shiny)
library(deSolve)
library(ggplot2)
library(pracma)

theme_set(theme_minimal())

# =======================
# Helpers PK/PD
# =======================

.parse_doses <- function(x) {
  d <- as.numeric(unlist(strsplit(x, ",")))
  d <- sort(d[!is.na(d) & d >= 0])
  unique(d)
}

nm_to_ngml <- function(x, MW = 105000) x * (MW/1000)

exposure_metrics <- function(time_days, conc, tau_days = 7) {
  stopifnot(length(time_days) == length(conc))
  tmax <- max(time_days)
  tstart <- max(0, tmax - tau_days)
  keep <- time_days >= tstart & time_days <= tmax
  t <- time_days[keep]
  c <- pmax(conc[keep], 1e-12)
  
  AUC_tau <- trapz(t, c)
  Cmax <- max(c)
  Tmax <- t[which.max(c)]
  Ctrough <- c[length(c)]
  
  idx <- which(t >= (tmax - 0.3 * tau_days))
  if (length(idx) >= 3 && all(c[idx] > 0)) {
    fit <- lm(log(c[idx]) ~ t[idx])
    lambda_z <- -as.numeric(coef(fit)[2])
    t_half <- if (is.finite(lambda_z) && lambda_z > 0) log(2) / lambda_z else NA_real_
  } else {
    t_half <- NA_real_
  }
  data.frame(AUC_tau, Cmax, Tmax, Ctrough, t_half)
}

cycle_metrics <- function(time_days, conc, tau_days = 7) {
  o <- order(time_days)
  t <- time_days[o]; c <- conc[o]
  key <- round(t, 9)
  keep <- !duplicated(key)
  t <- t[keep]; c <- c[keep]
  
  tmax <- max(t)
  ncyc <- floor(tmax / tau_days)
  if (ncyc < 1) return(data.frame())
  
  res <- lapply(seq_len(ncyc), function(k){
    t0 <- (k - 1) * tau_days
    t1 <- k * tau_days
    idx <- t >= t0 & t <= t1 + 1e-12
    tt <- t[idx]; cc <- pmax(c[idx], 1e-12)
    data.frame(cycle = k, start = t0, end = t1,
               AUC = trapz(tt, cc),
               Ctrough = cc[length(cc)])
  })
  do.call(rbind, res)
}

# =======================
# 1) QSP HUMANO
# =======================
simulate_human <- function(dose_ugkg = 0.1, days_total = 45,
                           CD3 = 4.5e5, mPcad = 28706, Tcellsp = 5000,
                           et_den = NULL,
                           Tumor_cells_t = 1e7) {
  peso <- 70; MW <- 105000
  dose_nmol <- ((peso * dose_ugkg * 1e-6) / MW) * 1e9
  V1 <- 40.2*peso; V2 <- 211*peso
  CL <- 9.61*peso; CLD <- 20.2*peso
  k_el <- CL/V1; k_12 <- CLD/V1; k_21 <- CLD/V2
  
  kon_CD3 <- 7.72; koff_CD3 <- 14.66
  kon_Pcad <- 9.7;  koff_Pcad <- 0.24
  P <- 334; Rcap <- 8; Rkrogh <- 75; epsilon <- 0.24
  TV <- 1
  Tcellsta_base <- 5e10
  Tcellsta <- if (is.null(et_den)) Tcellsta_base else (Tumor_cells_t / et_den)
  D <- 0.022; Rtumor <- 1; Kint <- 0.01728/24
  k_degcx <- 0.00915; kdeg <- 0.015; sPcad <- 1.1
  
  ksyn <- kdeg * sPcad
  totCD3t <- (Tcellsta * CD3) / (6.023e23) * 1e9
  totCD3p <- (Tcellsp * CD3) / (6.023e23) * 1e9
  totPCAD <- (Tumor_cells_t * mPcad) / (6.023e23) * 1e9
  
  tumor_disposition <- function(t, y) {
    (((2*P*Rcap)/(Rkrogh^2)) + ((6*D)/(Rtumor^2))) * (y[1] - (y[3]/epsilon))
  }
  
  model <- function(t, y, parms) {
    td <- tumor_disposition(t, y)
    with(as.list(c(y)), {
      dy1 <- -k_el*y[1] - k_12*y[1] + k_21*y[2]*(V2/V1) -
        kon_CD3*y[1]*(totCD3p - y[4]) + koff_CD3*y[4] -
        kon_Pcad*y[1]*y[8] + koff_Pcad*y[9] -
        td*(TV/V1)
      dy2 <-  k_12*y[1]*(V1/V2) - k_21*y[2]
      dy3 <-  td - kon_Pcad*y[3]*((totPCAD - y[5] - y[6])/epsilon) + koff_Pcad*y[6] -
        kon_CD3*y[3]*((totCD3t - y[7] - y[6])/epsilon) + koff_CD3*y[7]
      dy4 <-  kon_CD3*y[1]*(totCD3p - y[4]) - koff_CD3*y[4]
      dy5 <-  kon_Pcad*y[3]*((totPCAD - y[5] - y[6])/epsilon) - koff_Pcad*y[5] -
        kon_CD3*y[5]*((totCD3t - y[7] - y[6])/epsilon) + koff_CD3*y[6] - Kint*y[5]
      dy6 <-  kon_CD3*y[5]*((totCD3t - y[7] - y[6])/epsilon) - koff_CD3*y[6] +
        kon_Pcad*y[7]*((totPCAD - y[5] - y[6])/epsilon) - koff_Pcad*y[6]
      dy7 <-  kon_CD3*y[3]*((totCD3t - y[7] - y[6])/epsilon) - koff_CD3*y[7] -
        kon_Pcad*y[7]*((totPCAD - y[5] - y[6])/epsilon) + koff_Pcad*y[6]
      dy8 <-  ksyn - kdeg*y[8] - kon_Pcad*y[1]*y[8] + koff_Pcad*y[9]
      dy9 <-  kon_Pcad*y[1]*y[8] - koff_Pcad*y[9] - k_degcx*y[9]
      list(c(dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9))
    })
  }
  
  times <- seq(0, 24*days_total, by = 1)
  yinit <- c((dose_nmol/V1)*1000, 0, 0, 0, 0, 0, 0, sPcad, 0)
  yout <- yinit; tout <- 0
  
  for (start in seq(0, 24*days_total - 1, by = 24*7)) {
    end_t <- min(start + 24*7, 24*days_total)
    times_sub <- seq(start, end_t, by = 1)
    if (length(times_sub) < 2) break
    
    out <- ode(y = yinit, times = times_sub, func = model, parms = NULL,
               method = "lsoda", rtol = 1e-6, atol = 1e-8)
    
    yinit <- out[nrow(out), -1]
    
    if (end_t < 24*days_total) {
      y_post <- yinit
      y_post[1] <- y_post[1] + (dose_nmol/V1)*1000
      tout <- c(tout, out[-1, 1], end_t)
      yout <- rbind(yout, out[-1, -1], y_post)
      yinit <- y_post
    } else {
      tout <- c(tout, out[-1, 1])
      yout <- rbind(yout, out[-1, -1])
    }
  }
  
  data.frame(time = tout/24, trimer = yout[,6], C1 = yout[,1])
}

# =======================
# 2) QSP RATÓN (tu bloque base + W-only)
# =======================
simulate_mouse <- function(dose_ugkg = 0.01, total_days = 45,
                           kc50 = 3e-3,       # nM del driver (Ce) para 50% kill
                           kmax = 1/24,       # 1/h máximo kill
                           kill_cap = 0.35/24,# 1/h techo diario
                           ke0  = 6/24,       # 1/h comp. efecto
                           h    = 1.3) {      # Hill
  # ===== Constantes y conversión de dosis =====
  peso <- 0.25
  MW   <- 105000
  dose <- (peso * dose_ugkg * 1e-6) / MW * 1e9 # nmol
  
  # ===== PK base ratón =====
  V1 <- 49.6*peso; V2 <- 60.7*peso
  CL <- 0.145*peso; CLD <- 4.95*peso
  k_el <- CL/V1; k_12 <- CLD/V1; k_21 <- CLD/V2
  
  # T-cells
  k_elT <- 0.01; k_12T <- 0.001; k_21T <- 0.002
  
  # Afinidades/binding
  kon_CD3 <- 9.72*1.6; koff_CD3 <- 6.66*1.6
  kon_Pcad <- 19.7*4;  koff_Pcad <- 8.24*4
  
  # Penetración tumoral
  P <- 334*1.2; D <- 0.022*1.5; Rcap <- 8; Rkrogh <- 75; Rtumor <- 1; epsilon <- 0.24
  
  # Targets / células
  mPcad <- 28706; CD3 <- 6e5
  Tcellsp0 <- 2.5e8; Tcelltm0 <- 8.49e6
  Tumor_cells <- 1e8; TV <- 2
  
  # Internalización
  Kint <- 0.91728*2.5
  
  # ===== Crecimiento (Simeoni) =====
  kg0 <- 0.012/24; kg <- 9.3/24
  psi <- 10; Mmax <- 5.8e3; tau <- 2.25/24
  
  # Proliferación T y activación con retardo (~5 días)
  ksyn <- 0.015*1.1
  P_rate <- ((0.014)/(4 + dose_ugkg) + 1.5e-5) * dose_ugkg
  tlag_days <- 5
  act_fun <- function(t){ 1/(1 + exp(-(t/24 - tlag_days)/0.5)) }
  
  # ===== Helpers =====
  tot_CD3t <- function(Tcellst, CD3) (Tcellst * CD3)/(6.023e23)*1e9
  tot_PCADt <- function(Tumor_cells_t, mPcad) (Tumor_cells_t * mPcad)/(6.023e23)*1e9
  Tcells_t <- function(t, y, P_rate){
    if (t < tlag_days*24) y[9] else y[9] * exp(P_rate*((t/24) - tlag_days))
  }
  tumor_disposition_fun <- function(t, y, P, Rcap, Rkrogh, D, Rtumor, epsilon){
    ((2*P*Rcap)/(Rkrogh^2) + (6*D)/(Rtumor^2)) * (y[1] - (y[3] / epsilon))
  }
  w_fun <- function(y) y[10] + y[11] + y[12] + y[13]
  
  # ===== Parámetros a pasar a la ODE =====
  parms <- list(CD3=CD3, mPcad=mPcad, P=P, Rcap=Rcap, Rkrogh=Rkrogh,
                D=D, Rtumor=Rtumor, epsilon=epsilon, V1=V1, V2=V2, TV=TV,
                k_el=k_el, k_12=k_12, k_21=k_21, kon_CD3=kon_CD3, koff_CD3=koff_CD3,
                kon_Pcad=kon_Pcad, koff_Pcad=koff_Pcad, Kint=Kint,
                k_elT=k_elT, k_12T=k_12T, k_21T=k_21T, ksyn=ksyn,
                kg0=kg0, kg=kg, Mmax=Mmax, psi=psi, tau=tau,
                P_rate=P_rate, Tumor_cells=Tumor_cells,
                kmax=kmax, kc50=kc50, ke0=ke0, kill_cap=kill_cap, h=h)
  
  # ===== ODE =====
  odefun <- function(t, y, parms){
    with(as.list(parms), {
      Tcellst_dyn <- Tcells_t(t, y, P_rate)
      TotCD3t <- tot_CD3t(Tcellst_dyn, CD3)
      TotPCAD <- tot_PCADt(Tumor_cells, mPcad)
      td <- tumor_disposition_fun(t, y, P, Rcap, Rkrogh, D, Rtumor, epsilon)
      
      Trimer <- max(y[6], 0)
      Ce     <- y[14]
      
      # Kill tipo Emax con Hill sobre el driver (Ce)
      kkill <- ((kmax) * (Ce^h) / (kc50^h + Ce^h + .Machine$double.eps))^0.5
      kill  <- act_fun(t) * pmin(kkill, kill_cap)
      
      W <- w_fun(y); epsW <- 1e-9
      dydt <- numeric(14)
      
      # PK central/periférico/tumor
      dydt[1] <- -k_el*y[1] - k_12*y[1] + k_21*y[2]*(V2/V1) -
        kon_CD3*y[1]*(((y[8]*CD3)/(6.023e23)*1e9) - y[4]) + koff_CD3*y[4] -
        td*(TV/V1)
      dydt[2] <-  k_12*y[1]*(V1/V2) - k_21*y[2]
      dydt[3] <-  td - kon_Pcad*y[3]*((TotPCAD - y[5] - y[6])/epsilon) + koff_Pcad*y[5] -
        kon_CD3*y[3]*((TotCD3t - y[7] - y[6])/epsilon) + koff_CD3*y[7]
      
      # Dímeros / trímero
      dydt[4] <-  kon_CD3*y[1]*(((y[8]*CD3)/(6.023e23)*1e9) - y[4]) - koff_CD3*y[4]
      dydt[5] <-  kon_Pcad*y[3]*((TotPCAD - y[5] - y[6])/epsilon) - koff_Pcad*y[5] -
        kon_CD3*y[5]*((TotCD3t - y[7] - y[6])/epsilon) + koff_CD3*y[6] - Kint*y[5]
      dydt[6] <-  kon_CD3*y[5]*((TotCD3t - y[7] - y[6])/epsilon) - koff_CD3*y[6] +
        kon_Pcad*y[7]*((TotPCAD - y[5] - y[6])/epsilon) - koff_Pcad*y[6]
      dydt[7] <-  kon_CD3*y[3]*((TotCD3t - y[7] - y[6])/epsilon) - koff_CD3*y[7] -
        kon_Pcad*y[7]*((TotPCAD - y[5] - y[6])/epsilon) + koff_Pcad*y[6]
      
      # T-cells
      dydt[8] <- -k_elT*y[8] - k_12T*y[8] + k_21T*y[9]*(TV/V1)
      dydt[9] <-  k_12T*y[8]*(V1/TV) - k_21T*y[9]
      
      # Crecimiento (Simeoni) y kill restado en M1
      growth <- ((kg0 * (1 - W/Mmax)) * max(y[10], epsW)) /
        ((1 + ((kg0/kg) * W)^psi)^(1/psi))
      dydt[10] <-  growth - kill * max(y[10], epsW)  # M1: crecimiento – kill
      dydt[11] <-  kill * max(y[10], epsW) - y[11]/tau
      dydt[12] <- (y[11] - y[12]) / tau
      dydt[13] <- (y[12] - y[13]) / tau
      
      # Compartimento de efecto: filtra el TRÍMERO (driver)
      dydt[14] <- ke0 * (Trimer - Ce)
      
      list(dydt)
    })
  }
  
  # ===== Redosis semanal y simulación =====
  interval_days  <- 7
  interval_hours <- interval_days*24
  n_doses <- ceiling(total_days/interval_days)
  times <- seq(0, interval_hours, by=1)
  
  # Estado inicial
  W0 <- 180
  y <- c((dose/V1)*1e3, 0, 0, 0, 0, 0, 0, Tcellsp0, Tcelltm0, W0, 0, 0, 0, 0)
  
  results <- NULL
  for (i in 1:n_doses) {
    out <- deSolve::ode(y=y, times=times, func=odefun, parms=parms,
                        method="lsoda", rtol=1e-9, atol=1e-12)
    y <- out[nrow(out), -1]
    if ((i * interval_hours) < (total_days * 24)) {
      y[1] <- y[1] + (dose/V1)*1e3
    }
    out[, "time"] <- out[, "time"] + (i - 1)*interval_hours
    results <- rbind(results, out[out[, "time"] <= total_days*24, , drop=FALSE])
  }
  
  df <- as.data.frame(results)
  colnames(df) <- c("time","C1","C2","C3","DCD3p","DPcadt","Trimer","DCD3t",
                    "Tcellsp","Tcelltm","M1","M2","M3","M4","Ce")
  df$days <- df$time/24
  df$W <- pmax(df$M1 + df$M2 + df$M3 + df$M4, 1e-9)
  df
}

# =======================
# NUEVO: Ratón 
# =======================
simulate_mouse_ce_sat <- function(dose_ugkg,
                                  total_days,
                                  kc50, kmax, ke0, tau, kg0,
                                  peso_kg, MW, Ce_sat) {
  
  dose <- (peso_kg * dose_ugkg) * (1e3 / MW)   # nmol
  
  # PK base
  V1 <- 49.6 * peso_kg; V2 <- 60.7 * peso_kg
  CL <- 0.145 * peso_kg; CLD <- 4.95 * peso_kg
  k_el <- CL/V1; k_12 <- CLD/V1; k_21 <- CLD/V2
  
  # T cells PK
  k_elT <- 0.01; k_12T <- 0.001; k_21T <- 0.002
  
  # Binding/transporte
  kon_CD3 <- 9.72*1.6; koff_CD3 <- 6.66*1.6
  kon_Pcad <- 19.7*4;  koff_Pcad <- 8.24*4
  P <- 334*1.2; D <- 0.022*1.5; Rcap <- 8; Rkrogh <- 75; Rtumor <- 1; epsilon <- 0.24
  
  # Targets/células
  mPcad <- 28706; CD3 <- 6e5
  Tcellsp0 <- 2.5e8; Tcelltm0 <- 8.49e6
  Tumor_cells <- 1e8
  
  # varios
  Kint <- 0.91728*2.5
  TV <- 180; psi <- 10; Mmax <- 5.8e3
  
  # activación T
  P_rate <- ((0.014)/(4 + dose_ugkg) + 1.5e-5) * dose_ugkg
  tlag_days <- 5
  Tcells_t <- function(t, y, P_rate) if (t < tlag_days*24) y[9] else y[9]*exp(P_rate*((t/24) - tlag_days))
  
  tot_CD3t  <- function(Tcellst, CD3) (Tcellst * CD3)/(6.023e23)*1e9
  tot_PCADt <- function(Tumor_cells_t, mPcad) (Tumor_cells_t * mPcad)/(6.023e23)*1e9
  tumor_disp_fun <- function(t, y) ((2*P*Rcap)/(Rkrogh^2) + (6*D)/(Rtumor^2)) * (y[1] - (y[3]/epsilon))
  
  parms <- list(CD3=CD3, mPcad=mPcad, P=P, Rcap=Rcap, Rkrogh=Rkrogh,
                D=D, Rtumor=Rtumor, epsilon=epsilon, V1=V1, V2=V2, TV=TV,
                k_el=k_el, k_12=k_12, k_21=k_21, kon_CD3=kon_CD3, koff_CD3=koff_CD3,
                kon_Pcad=kon_Pcad, koff_Pcad=koff_Pcad, Kint=Kint,
                k_elT=k_elT, k_12T=k_12T, k_21T=k_21T,
                kg0=kg0, Mmax=Mmax, psi=psi, tau=tau,
                P_rate=P_rate, Tumor_cells=Tumor_cells,
                kmax=kmax, kc50=kc50, ke0=ke0, Ce_sat=Ce_sat)
  
  odefun <- function(t, y, parms){
    with(as.list(parms), {
      Tcellst_dyn <- Tcells_t(t, y, P_rate)
      TotCD3t <- tot_CD3t(Tcellst_dyn, CD3)
      TotPCAD <- tot_PCADt(Tumor_cells, mPcad)
      td <- tumor_disp_fun(t, y)
      
      Trimer <- max(y[6], 0)
      Ce     <- y[14]
      Ce_eff <- Ce / (1 + Ce/Ce_sat)      
      kkill  <- (kmax * Ce_eff) / (kc50 + Ce_eff)
      
      dydt <- numeric(14)
      # PK + uniones
      dydt[1] <- -k_el*y[1] - k_12*y[1] + k_21*y[2]*(V2/V1) -
        kon_CD3*y[1]*(((y[8]*CD3)/(6.023e23)*1e9) - y[4]) + koff_CD3*y[4] -
        td*(TV/V1)
      dydt[2] <-  k_12*y[1]*(V1/V2) - k_21*y[2]
      dydt[3] <-  td - kon_Pcad*y[3]*((TotPCAD - y[5] - y[6])/epsilon) + koff_Pcad*y[5] -
        kon_CD3*y[3]*((TotCD3t - y[7] - y[6])/epsilon) + koff_CD3*y[7]
      dydt[4] <-  kon_CD3*y[1]*(((y[8]*CD3)/(6.023e23)*1e9) - y[4]) - koff_CD3*y[4]
      dydt[5] <-  kon_Pcad*y[3]*((TotPCAD - y[5] - y[6])/epsilon) - koff_Pcad*y[5] -
        kon_CD3*y[5]*((TotCD3t - y[7] - y[6])/epsilon) + koff_CD3*y[6] - Kint*y[5]
      dydt[6] <-  kon_CD3*y[5]*((TotCD3t - y[7] - y[6])/epsilon) - koff_CD3*y[6] +
        kon_Pcad*y[7]*((TotPCAD - y[5] - y[6])/epsilon) - koff_Pcad*y[6]
      dydt[7] <-  kon_CD3*y[3]*((TotCD3t - y[7] - y[6])/epsilon) - koff_CD3*y[7] -
        kon_Pcad*y[7]*((TotPCAD - y[5] - y[6])/epsilon) + koff_Pcad*y[6]
      
      # T cells PK
      dydt[8] <- -k_elT*y[8] - k_12T*y[8] + k_21T*y[9]*(TV/V1)
      dydt[9] <-  k_12T*y[8]*(V1/TV) - k_21T*y[9]
      
      # Tumor
      dydt[10] <- kg0 * y[10] - kkill * y[10]
      dydt[11] <- kkill * y[10] - y[11]/tau
      dydt[12] <- (y[11] - y[12]) / tau
      dydt[13] <- (y[12] - y[13]) / tau
      
      # Compartimento de efecto (dientes de sierra)
      dydt[14] <- ke0 * (Trimer - Ce)
      
      list(dydt)
    })
  }
  
  # redosis semanal
  interval_hours <- 7*24
  n_doses <- ceiling(total_days/7)
  times <- seq(0, interval_hours, by = 1)
  
  y <- c((dose/V1)*1e3, 0, 0, 0, 0, 0, 0,
         2.5e8, 8.49e6,
         180, 0, 0, 0,
         0)
  
  results <- NULL
  for (i in 1:n_doses) {
    out <- deSolve::ode(y=y, times=times, func=odefun, parms=parms,
                        method="lsoda", rtol=1e-9, atol=1e-12)
    y <- out[nrow(out), -1]
    if ((i * interval_hours) < (total_days * 24)) {
      y[1] <- y[1] + (dose/V1)*1e3
    }
    out[, "time"] <- out[, "time"] + (i - 1)*interval_hours
    results <- rbind(results, out[out[, "time"] <= total_days*24, , drop=FALSE])
  }
  
  df <- as.data.frame(results)
  colnames(df) <- c("time","C1","C2","C3","DCD3p","DPcadt","Trimer","DCD3t",
                    "Tcellsp","Tcelltm","M1","M2","M3","M4","Ce")
  df$days <- df$time/24
  df$W <- pmax(df$M1 + df$M2 + df$M3 + df$M4, 1e-9)
  df
}

simulate_mouse_ce_sat_sched <- function(dose_ugkg,
                                        dose_times_days = c(0),   # vector con los días de administración
                                        total_days,
                                        kc50, kmax, ke0, tau, kg0,
                                        peso_kg, MW, Ce_sat) {
  
  dose <- (peso_kg * dose_ugkg) * (1e3 / MW)   
  
  V1 <- 49.6 * peso_kg; V2 <- 60.7 * peso_kg
  CL <- 0.145 * peso_kg; CLD <- 4.95 * peso_kg
  k_el <- CL/V1; k_12 <- CLD/V1; k_21 <- CLD/V2
  
  k_elT <- 0.01; k_12T <- 0.001; k_21T <- 0.002
  kon_CD3 <- 9.72*1.6; koff_CD3 <- 6.66*1.6
  kon_Pcad <- 19.7*4;  koff_Pcad <- 8.24*4
  P <- 334*1.2; D <- 0.022*1.5; Rcap <- 8; Rkrogh <- 75; Rtumor <- 1; epsilon <- 0.24
  
  mPcad <- 28706; CD3 <- 6e5
  Tcellsp0 <- 2.5e8; Tcelltm0 <- 8.49e6
  Tumor_cells <- 1e8
  
  Kint <- 0.91728*2.5
  TV <- 180; psi <- 10; Mmax <- 5.8e3
  
  P_rate <- ((0.014)/(4 + dose_ugkg) + 1.5e-5) * dose_ugkg
  tlag_days <- 5
  Tcells_t <- function(t, y, P_rate) if (t < tlag_days*24) y[9] else y[9]*exp(P_rate*((t/24) - tlag_days))
  
  tot_CD3t  <- function(Tcellst, CD3) (Tcellst * CD3)/(6.023e23)*1e9
  tot_PCADt <- function(Tumor_cells_t, mPcad) (Tumor_cells_t * mPcad)/(6.023e23)*1e9
  tumor_disp_fun <- function(t, y) ((2*P*Rcap)/(Rkrogh^2) + (6*D)/(Rtumor^2)) * (y[1] - (y[3]/epsilon))
  
  parms <- list(CD3=CD3, mPcad=mPcad, P=P, Rcap=Rcap, Rkrogh=Rkrogh,
                D=D, Rtumor=Rtumor, epsilon=epsilon, V1=V1, V2=V2, TV=TV,
                k_el=k_el, k_12=k_12, k_21=k_21, kon_CD3=kon_CD3, koff_CD3=koff_CD3,
                kon_Pcad=kon_Pcad, koff_Pcad=koff_Pcad, Kint=Kint,
                k_elT=k_elT, k_12T=k_12T, k_21T=k_21T,
                kg0=kg0, Mmax=Mmax, psi=psi, tau=tau,
                P_rate=P_rate, Tumor_cells=Tumor_cells,
                kmax=kmax, kc50=kc50, ke0=ke0, Ce_sat=Ce_sat)
  
  odefun <- function(t, y, parms){
    with(as.list(parms), {
      Tcellst_dyn <- Tcells_t(t, y, P_rate)
      TotCD3t <- tot_CD3t(Tcellst_dyn, CD3)
      TotPCAD <- tot_PCADt(Tumor_cells, mPcad)
      td <- tumor_disp_fun(t, y)
      
      Trimer <- max(y[6], 0)
      Ce     <- y[14]
      Ce_eff <- Ce / (1 + Ce/Ce_sat)
      kkill  <- (kmax * Ce_eff) / (kc50 + Ce_eff)
      
      dydt <- numeric(14)
      dydt[1] <- -k_el*y[1] - k_12*y[1] + k_21*y[2]*(V2/V1) -
        kon_CD3*y[1]*(((y[8]*CD3)/(6.023e23)*1e9) - y[4]) + koff_CD3*y[4] -
        td*(TV/V1)
      dydt[2] <-  k_12*y[1]*(V1/V2) - k_21*y[2]
      dydt[3] <-  td - kon_Pcad*y[3]*((TotPCAD - y[5] - y[6])/epsilon) + koff_Pcad*y[5] -
        kon_CD3*y[3]*((TotCD3t - y[7] - y[6])/epsilon) + koff_CD3*y[7]
      dydt[4] <-  kon_CD3*y[1]*(((y[8]*CD3)/(6.023e23)*1e9) - y[4]) - koff_CD3*y[4]
      dydt[5] <-  kon_Pcad*y[3]*((TotPCAD - y[5] - y[6])/epsilon) - koff_Pcad*y[5] -
        kon_CD3*y[5]*((TotCD3t - y[7] - y[6])/epsilon) + koff_CD3*y[6] - Kint*y[5]
      dydt[6] <-  kon_CD3*y[5]*((TotCD3t - y[7] - y[6])/epsilon) - koff_CD3*y[6] +
        kon_Pcad*y[7]*((TotPCAD - y[5] - y[6])/epsilon) - koff_Pcad*y[6]
      dydt[7] <-  kon_CD3*y[3]*((TotCD3t - y[7] - y[6])/epsilon) - koff_CD3*y[7] -
        kon_Pcad*y[7]*((TotPCAD - y[5] - y[6])/epsilon) + koff_Pcad*y[6]
      
      dydt[10] <- kg0 * y[10] - kkill * y[10]
      dydt[11] <- kkill * y[10] - y[11]/tau
      dydt[12] <- (y[11] - y[12]) / tau
      dydt[13] <- (y[12] - y[13]) / tau
      
      dydt[14] <- ke0 * (Trimer - Ce)
      list(dydt)
    })
  }
  

  dose_times_days <- sort(unique(dose_times_days))
  dose_times_days <- dose_times_days[dose_times_days < total_days]
  
  y <- c(0,0,0,0,0,0,0, 2.5e8, 8.49e6, 180, 0,0,0, 0)
  
  cuts <- sort(unique(c(0, dose_times_days, total_days)))
  results <- NULL
  
  for (i in seq_len(length(cuts)-1)) {
    t0 <- cuts[i]; t1 <- cuts[i+1]
    times <- seq(t0*24, t1*24, by = 1)
    out <- deSolve::ode(y=y, times=times, func=odefun, parms=parms,
                        method="lsoda", rtol=1e-9, atol=1e-12)
    results <- rbind(results, out[-1, , drop=FALSE])
    y <- out[nrow(out), -1]
    if (t1 %in% dose_times_days) y[1] <- y[1] + (dose/V1)*1e3
  }
  
  df <- as.data.frame(results)
  colnames(df) <- c("time","C1","C2","C3","DCD3p","DPcadt","Trimer","DCD3t",
                    "Tcellsp","Tcelltm","M1","M2","M3","M4","Ce")
  df$days <- df$time/24
  df$W <- pmax(df$M1 + df$M2 + df$M3 + df$M4, 1e-9)
  df
}


# =======================
# Comparativa 3 dosis 
# =======================
compare_three_doses_ce <- function(
    doses_ugkg = c(0.01, 0.1, 1),
    total_days = 45,
    kc50 = 3e-5,    
    kmax = 24,       
    ke0  = 10,    
    tau  = 2.25/24,
    kg0  = 0.001,
    peso_kg = 0.025,
    MW = 105000,
    Ce_sat = 1.5e-8 # bajar => más solape entre 0.1 y 1
){
  dfs <- lapply(doses_ugkg, function(d){
    df <- simulate_mouse_ce_sat(d, total_days, kc50, kmax, ke0, tau, kg0,
                                peso_kg, MW, Ce_sat)
    df$Dosis <- paste0(formatC(d, format="fg", digits=3), " µg/kg")
    df
  })
  df_all <- do.call(rbind, dfs)
  
  levs <- paste0(formatC(sort(unique(doses_ugkg)), format="fg", digits=3), " µg/kg")
  df_all$Dosis <- factor(df_all$Dosis, levels = levs)
  
  df_all
}

# =======================
# UI
# =======================
ui <- fluidPage(
  titlePanel("Modelado y predicción PK/PD para PF-06671008"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("species", "Especie", choices = c("Humano","Ratón"), inline = TRUE),
      hr(),
      # ---- Controles HUMANOS ----
      conditionalPanel(
        condition = "input.tabs == 'Fármaco central (C1)'",
        textInput("pk_doses", "Dosis PK (µg/kg), separadas por coma (mín. 3):",
                  value = "0.01,0.1,1")
      ),
      conditionalPanel(
        condition = "input.tabs == 'Trímero vs mPcad' && input.species == 'Humano'",
        textInput("mpcad_vals", "mPcad (rec/cél), separados por coma",
                  value = "1000,5000,15000,28706")
      ),
      conditionalPanel(
        condition = "input.tabs == 'Trímero vs CD3' && input.species == 'Humano'",
        textInput("cd3_vals", "CD3 (rec/cél), separados por coma",
                  value = "4e4,7e4,1e5,1.2e5")
      ),
      conditionalPanel(
        condition = "input.tabs == 'Trímero vs E:T' && input.species == 'Humano'",
        textInput("et_user", "E:T (formato 1:x, separados por coma)",
                  value = "1:15, 1:150, 1:1500"),
        numericInput("dose_et", "Dosis para E:T (µg/kg)", value = 0.1, min = 1e-4, step = 0.01)
      ),
      conditionalPanel(
        condition = "input.tabs == 'Trímero vs Dosis' && input.species == 'Humano'",
        textInput("doses_h_vec", "Dosis (µg/kg), separadas por coma",
                  value = "0.01,0.1,1")
      ),
      # ---- CONTROLES RATÓN ----
      conditionalPanel(
        condition = "input.tabs == 'Ratón – C1/C2/C3 (7 días tras redosis)'",
        numericInput("dose_readmin_m", "Dosis para C1/C2/C3 (µg/kg)", value = 0.1, min = 1e-4, step = 0.01)
      ),
      conditionalPanel(
        condition = "input.tabs == 'Ratón – Análisis PK – Métricas'",
        numericInput("dose_metrics_m", "Dosis (µg/kg)", value = 0.1, min = 1e-4, step = 0.01),
        numericInput("days_metrics_m", "Días totales", value = 45, min = 7, step = 7)
      ),
      conditionalPanel(
        condition = "input.tabs == 'Ratón – PK – Normalizado por dosis'",
        textInput("doses_dn_m", "Dosis a comparar (µg/kg), separadas por coma",
                  value = "0.01,0.03,0.1,0.3,1")
      ),
      conditionalPanel(
        condition = "input.tabs == 'Ratón – Volumen tumoral (dosis)'",
        textInput("doses_w_compare", "Dosis (µg/kg), separadas por coma",
                  value = "0,0.01,0.1,1,3"),
        numericInput("single_day_b", "Dosis única – día", value = 14, min = 1, step = 1)
      )
    ),  
    
    mainPanel(
      uiOutput("tabs_ui")
    )
  )     
)      


# =======================
# SERVIDOR
# =======================
server <- function(input, output, session){
  
  # Tabs dinámicos por especie
  output$tabs_ui <- renderUI({
    if (input$species == "Humano") {
      tabsetPanel(id = "tabs",
                  tabPanel("Fármaco central (C1)", plotOutput("plot_pk")),
                  tabPanel("Trímero vs mPcad", plotOutput("plot_mPcad")),
                  tabPanel("Trímero vs CD3", plotOutput("plot_CD3")),
                  tabPanel("Trímero vs E:T", plotOutput("plot_trimer_et")),
                  tabPanel("Trímero vs Dosis", plotOutput("plot_trimer_dose")),
      )
    } else {
      tabsetPanel(id = "tabs",
                  tabPanel("Ratón – Volumen tumoral (dosis)",
                           plotOutput("plot_mouse_vol_compare_ce", height = "340px")
                  ),
                  tabPanel("Ratón – C1/C2/C3 (7 días tras redosis)",
                           fluidRow(
                             column(6, plotOutput("plot_mouse_c1_readmin", height = "260px")),
                             column(6, plotOutput("plot_mouse_c2_readmin", height = "260px"))
                           ),
                           fluidRow(
                             column(12, plotOutput("plot_mouse_c3_readmin", height = "280px"))
                           )
                  ),
                  tabPanel("Ratón – Análisis PK – Métricas",
                           h4("C1 (central) – último ciclo"),
                           tableOutput("tbl_metrics_m_c1"),
                           br(), h4("C3 (tumor) – último ciclo"),
                           tableOutput("tbl_metrics_m_c3")
                           ),
                  tabPanel("Ratón – PK – Normalizado por dosis",
                           fluidRow(
                             column(6, plotOutput("plot_dn_cmax_m", height="300px")),
                             column(6, plotOutput("plot_dn_auc_m",  height="300px"))
                           )
                  )
      )
    }
  })
  
  # --------- HUMANO: PK central (C1 en ng/mL)
  output$plot_pk <- renderPlot({
    req(input$species == "Humano")
    doses <- .parse_doses(input$pk_doses)
    validate( need(length(doses) >= 3, "Introduce al menos tres dosis válidas (separadas por coma).") )
    
    dfs <- lapply(doses, function(d){
      df <- simulate_human(dose_ugkg = d, days_total = 45,
                           CD3 = 4.5e5, mPcad = 28706, Tcellsp = 5000)
      df$C1_ngml <- nm_to_ngml(df$C1)
      df$Dosis <- paste0(d, " µg/kg")
      df
    })
    df_all <- do.call(rbind, dfs)
    
    ggplot(df_all, aes(x = time, y = C1_ngml, color = Dosis)) +
      geom_line(linewidth = 1.2, show.legend = TRUE) +
      scale_y_log10() +
      labs(title = "Humano: Fármaco en compartimento central (C1)",
           x = "Días", y = "C1 (ng/mL)", color = "Dosis (µg/kg)")
  })
  
  # --------- HUMANO: Trímero vs mPcad
  output$plot_mPcad <- renderPlot({
    req(input$species == "Humano")
    m_vals <- as.numeric(unlist(strsplit(input$mpcad_vals, ",")))
    m_vals <- m_vals[!is.na(m_vals) & m_vals > 0]
    dfs <- lapply(m_vals, function(m) {
      df <- simulate_human(dose_ugkg = 0.1, days_total = 45,
                           CD3 = 4.5e5, mPcad = m, Tcellsp = 5000)
      df$mPcad <- factor(paste0(m, " rec/cél"), levels = paste0(m_vals, " rec/cél"))
      df
    })
    validate( need(length(dfs) > 0, "Introduce valores válidos de mPcad") )
    ggplot(do.call(rbind, dfs), aes(x = time, y = trimer, color = mPcad)) +
      geom_line(linewidth = 1.2, show.legend = TRUE) + scale_y_log10() +
      labs(title = "Humano: Trímero vs mPcad", x = "Días", y = "Trímero (nM)", color = "mPcad")
  })
  
  # --------- HUMANO: Trímero vs CD3
  output$plot_CD3 <- renderPlot({
    req(input$species == "Humano")
    c_vals <- as.numeric(unlist(strsplit(input$cd3_vals, ",")))
    c_vals <- c_vals[!is.na(c_vals) & c_vals > 0]
    dfs <- lapply(c_vals, function(c) {
      df <- simulate_human(dose_ugkg = 0.1, days_total = 45,
                           CD3 = c, mPcad = 28706, Tcellsp = 5000)
      df$CD3 <- paste0("CD3=", formatC(c, format = "e", digits = 1))
      df
    })
    validate( need(length(dfs) > 0, "Introduce valores válidos de CD3") )
    ggplot(do.call(rbind, dfs), aes(x = time, y = trimer, color = CD3)) +
      geom_line(linewidth = 1.2, show.legend = TRUE) + scale_y_log10() +
      labs(title = "Humano: Trímero vs CD3", x = "Días", y = "Trímero (nM)", color = "CD3")
  })
  
  # --------- HUMANO: Trímero vs E:T
  output$plot_trimer_et <- renderPlot({
    req(input$species == "Humano")
    dose_use <- if (!is.null(input$dose_et)) as.numeric(input$dose_et) else 0.1
    validate( need(is.finite(dose_use) && dose_use > 0, "Dosis inválida") )
    
    raw_tokens <- unlist(strsplit(ifelse(is.null(input$et_user), "", input$et_user), "[,;\\s]+"))
    raw_tokens <- raw_tokens[nchar(raw_tokens) > 0]
    
    parse_den <- function(s){
      s <- trimws(s)
      s2 <- gsub("[^0-9:\\.]", "", s)
      if (grepl(":", s2)) {
        den <- sub(".*:", "", s2)
      } else if (grepl("^1\\.[0-9]{3,}$", s2)) {
        den <- gsub("\\.", "", s2); den <- sub("^1", "", den)
      } else {
        den <- s2
      }
      as.numeric(den)
    }
    
    dens <- vapply(raw_tokens, parse_den, numeric(1))
    dens <- unique(dens[is.finite(dens) & dens > 0])
    
    validate( need(length(dens) >= 1, "Especifica al menos un E:T válido (por ejemplo, 1:15).") )
    
    labels <- paste0("1:", formatC(dens, format = "fg", digits = 6))
    names(dens) <- labels
    
    dfs <- lapply(seq_along(dens), function(i){
      den <- dens[i]; lbl <- names(dens)[i]
      df <- simulate_human(dose_ugkg = dose_use, days_total = 45, Tumor_cells_t = 8e8,
                           CD3 = 4.5e5, mPcad = 28706, Tcellsp = 5000,
                           et_den = den)
      df$ET <- lbl; df
    })
    df_all <- do.call(rbind, dfs)
    df_all$trimer <- pmax(df_all$trimer, 1e-12)
    
    ggplot(df_all, aes(x = time, y = trimer, color = ET)) +
      geom_line(linewidth = 1.2, show.legend = TRUE) +
      scale_y_log10() +
      labs(title = "Humano: Trímero vs E:T",
           x = "Días", y = "Trímero (nM)", color = "E:T")
  })
  
  # --------- HUMANO: Trímero vs Dosis
  output$plot_trimer_dose <- renderPlot({
    req(input$species == "Humano")
    dosis_values <- .parse_doses(input$doses_h_vec)
    dfs <- lapply(dosis_values, function(d) {
      df <- simulate_human(dose_ugkg = d, days_total = 45,
                           CD3 = 4.5e5, mPcad = 28706, Tcellsp = 5000)
      df$Dosis <- paste0(formatC(d, format = "fg", digits = 3), " \u00B5g/kg"); df
    })
    validate( need(length(dfs) > 0, "Introduce dosis válidas") )
    df_all <- do.call(rbind, dfs)
    df_all$Dosis <- factor(df_all$Dosis, levels = unique(df_all$Dosis))
    ggplot(df_all, aes(x = time, y = trimer, color = Dosis)) +
      geom_line(linewidth = 1.2, show.legend = TRUE) + scale_y_log10() +
      labs(title = "Humano: Trímero vs Dosis", x = "Días", y = "Trímero (nM)",
           color = "Dosis (µg/kg)") +
      theme(legend.position = "right")
  })
  
  # --------- RATÓN: Volumen tumoral – comparativa dosis
  
  output$plot_mouse_vol_compare_ce <- renderPlot({
    req(input$species == "Ratón")
    
    doses <- .parse_doses(if (is.null(input$doses_w_compare)) "0.01,0.1,1" else input$doses_w_compare)
    validate(need(length(doses) >= 2, "Introduce al menos dos dosis válidas."))
    
    # Curvas con redosis semanal (las “normales”)
    df_main <- compare_three_doses_ce(doses_ugkg = doses)
    
    
    dayB <- if (is.null(input$single_day_b)) 14 else input$single_day_b
    dfB <- simulate_mouse_ce_sat_sched(
      dose_ugkg = 0.1,
      dose_times_days = c(dayB),
      total_days = 45,
      kc50 = 3e-5, kmax = 24, ke0 = 10, tau = 2.25/24,
      kg0  = 0.001, peso_kg = 0.025, MW = 105000, Ce_sat = 1.5e-8
    )
    dfB$Dosis <- paste0("0.1 µg/kg (única día ", dayB, ")")
    
    
    df_all <- rbind(
      df_main[, c("days","W","Dosis")],
      dfB[,    c("days","W","Dosis")]
    )
    
    num_levels   <- paste0(formatC(sort(unique(doses)), format = "fg", digits = 3), " µg/kg")
    extra_levels <- unique(dfB$Dosis)
    df_all$Dosis <- factor(df_all$Dosis, levels = c(num_levels, extra_levels))
    
    ggplot(df_all, aes(x = days, y = W, color = Dosis)) +
      geom_line(linewidth = 1.2) +
      labs(title = "Ratón – Volumen tumoral (Ce con saturación)",
           x = "Tiempo (días)", y = "W (mm³)", color = "Dosis") +
      theme_minimal()
  })
  
  
  # --------- RATÓN: C1/C2/C3  (ventana 7–14 días tras 2ª dosis)
  get_mouse_window <- function(df) {
    win <- subset(df, days >= 7 & days <= 13.9)
    win$days_since_redose <- win$days - 7
    win
  }
  
  output$plot_mouse_c1_readmin <- renderPlot({
    req(input$species == "Ratón")
    dose_val <- if (is.null(input$dose_readmin_m)) 0.1 else input$dose_readmin_m
    df <- simulate_mouse(dose_ugkg = dose_val, total_days = 45)
    win <- get_mouse_window(df)
    win$Dosis <- paste0(formatC(dose_val, format = "fg", digits = 3), " µg/kg")
    ggplot(win, aes(x = days_since_redose, y = C1, color = Dosis)) +
      geom_line(linewidth = 1.2, show.legend = TRUE) +
      labs(title = "Compartimento central (C1)", x = "Días desde redosis", y = "C1 (nM)",
           color = "Dosis (µg/kg)")
  })
  
  output$plot_mouse_c2_readmin <- renderPlot({
    req(input$species == "Ratón")
    dose_val <- if (is.null(input$dose_readmin_m)) 0.1 else input$dose_readmin_m
    df <- simulate_mouse(dose_ugkg = dose_val, total_days = 45)
    win <- get_mouse_window(df)
    win$Dosis <- paste0(formatC(dose_val, format = "fg", digits = 3), " µg/kg")
    ggplot(win, aes(x = days_since_redose, y = C2, color = Dosis)) +
      geom_line(linewidth = 1.2, show.legend = TRUE) +
      labs(title = "Compartimento periférico (C2)", x = "Días desde redosis", y = "C2 (nM)",
           color = "Dosis (µg/kg)")
  })
  
  output$plot_mouse_c3_readmin <- renderPlot({
    req(input$species == "Ratón")
    dose_val <- if (is.null(input$dose_readmin_m)) 0.1 else input$dose_readmin_m
    df <- simulate_mouse(dose_ugkg = dose_val, total_days = 45)
    win <- get_mouse_window(df)
    win$Dosis <- paste0(formatC(dose_val, format = "fg", digits = 3), " µg/kg")
    ggplot(win, aes(x = days_since_redose, y = C3, color = Dosis)) +
      geom_line(linewidth = 1.2, show.legend = TRUE) +
      labs(title = "Microambiente tumoral (C3)", x = "Días desde redosis", y = "C3 (nM)",
           color = "Dosis (µg/kg)")
  })
  
  # ----------------------- RATÓN: métricas y acumulación
  output$tbl_metrics_m_c1 <- renderTable({
    req(input$species == "Ratón")
    df <- simulate_mouse(dose_ugkg = input$dose_metrics_m, total_days = input$days_metrics_m)
    exposure_metrics(df$time/24, df$C1, tau_days = 7)
  }, digits = 6)
  
  output$tbl_metrics_m_c3 <- renderTable({
    req(input$species == "Ratón")
    df <- simulate_mouse(dose_ugkg = input$dose_metrics_m, total_days = input$days_metrics_m)
    exposure_metrics(df$time/24, df$C3, tau_days = 7)
  }, digits = 6)
  
  # --------- RATÓN: PK normalizado por dosis
  output$plot_dn_cmax_m <- renderPlot({
    req(input$species == "Ratón")
    doses <- .parse_doses(input$doses_dn_m)
    validate(need(length(doses) >= 2, "Introduce al menos dos dosis."))
    
    dd <- do.call(rbind, lapply(doses, function(d){
      df <- simulate_mouse(dose_ugkg = d, total_days = 45)
      data.frame(dose = d, Cmax = max(df$C1))
    }))
    dd$DN_Cmax <- dd$Cmax / dd$dose
    
    ggplot(dd, aes(x = dose, y = DN_Cmax)) +
      geom_line(aes(group = 1), color = "black", linewidth = 0.8, linetype = "dotted") +
      geom_point(aes(color = factor(dose)), size = 3) +
      scale_x_log10() + scale_y_log10() +
      labs(title = "Cmax normalizado por dosis",
           x = expression("Dosis ("*mu*"g/kg)"),
           y = expression("Cmax/dosis (nM·(µg/kg)^"-1*")"),
           color = "Dosis (µg/kg)") +
      guides(color = guide_legend(override.aes = list(linewidth = 1.5)))
  })
  
  output$plot_dn_auc_m <- renderPlot({
    req(input$species == "Ratón")
    doses <- .parse_doses(input$doses_dn_m)
    validate(need(length(doses) >= 2, "Introduce al menos dos dosis."))
    
    dd <- do.call(rbind, lapply(doses, function(d){
      df <- simulate_mouse(dose_ugkg = d, total_days = 45)
      data.frame(dose = d, AUC = pracma::trapz(df$time/24, df$C1))
    }))
    dd$DN_AUC <- dd$AUC / dd$dose
    
    ggplot(dd, aes(x = dose, y = DN_AUC)) +
      geom_line(aes(group = 1), color = "black", linewidth = 0.8, linetype = "dotted") +
      geom_point(aes(color = factor(dose)), size = 3) +
      scale_x_log10() + scale_y_log10() +
      labs(title = "AUC normalizada por dosis",
           x = expression("Dosis ("*mu*"g/kg)"),
           y = expression("AUC/dosis (nM·día·(µg/kg)^"-1*")"),
           color = "Dosis (µg/kg)") +
      guides(color = guide_legend(override.aes = list(linewidth = 1.5)))
  })
}

# Lanzar app
shinyApp(ui = ui, server = server)
