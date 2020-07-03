
remove(list = ls())

library(vars)
library(pls)
library(missMDA)
library(seasonal)
library(CCA)
library(mixOmics)
library(tempdisagg)
library(forecast)
library(seasonal)

# difference
d <- 12L
# icc_vars <- "COMP"  # "ICC" / "COMP"
  
# path
path <- "C:/Users/jesus.lopezp/Desktop/ICC/ICC_2020/"
source(paste(path, "Functions.R", sep = ""))

include_sav <- TRUE 
dlog_icc <- FALSE

my_date_econ_start <- "2005/01"
ts_econ_start <- as.numeric(c(substr(my_date_econ_start,1,4),
                              substr(my_date_econ_start,6,7)))

# data
icc_nsa <- read.csv(paste(path, "inputs/ICC_NSA.csv", sep = ""), row.names = 1)
icc_sa <- read.csv(paste(path, "inputs/ICC_SA.csv", sep = ""), row.names = 1)
comp_nsa <- read.csv(paste(path, "inputs/COMP_NSA.csv", sep = ""), row.names = 1)
comp_sa <- read.csv(paste(path, "inputs/COMP_SA.csv", sep = ""), row.names = 1)
etco_icc <- read.csv(paste(path, "inputs/ETCO_ICC.csv", sep = ""), row.names = 1)
etco_comp <- read.csv(paste(path, "inputs/ETCO_COMP.csv", sep = ""), row.names = 1)

# leemos variables del consumo
econ_nsa <- read.csv(paste(path, "inputs/Econ3.csv", sep = ""), row.names = 1)
econ_cat <- read.csv(paste(path, "inputs/econ_cat.csv", sep = ""))

# leemos variables de ahorro
sav <- read.csv(paste(path, "inputs/sav.csv", sep = ""), row.names = 1)


# modelos de ajuste estacional icc INEGI 
modelos_icc <- read.csv(paste(path, "inputs/modelos-SA-icc.csv", sep = ""), 
                    stringsAsFactors = FALSE)
modelos_comp <- read.csv(paste(path, "inputs/modelos-SA-comp.csv", sep = ""), 
                        stringsAsFactors = FALSE)
modelos_cons <- read.csv(paste(path, "inputs/modelos-SA-consumo.csv", sep = ""), 
                         stringsAsFactors = FALSE)

# leemos datos de variables auxiliares
var_aux <- read.csv(paste(path, "inputs/var_aux.csv", sep = ""), 
                    row.names = 1)

# juntamos enco y etco
index_icc_nd <- which(is.na(icc_nsa[,1]))
index_comp_nd <- which(is.na(comp_nsa[,1]))
icc_nsa <- rbind(icc_nsa[-index_icc_nd,], etco_icc)
comp_nsa <- rbind(comp_nsa[-index_comp_nd,], etco_comp)

# obtenemos vectores de fechas Utiles
dates_icc <- rownames(icc_nsa)
dates_comp <- row.names(comp_nsa)
dates_econ <- row.names(econ_nsa)
index_start_icc <- which(dates_icc == dates_comp[1])
dates_var_aux <- row.names(var_aux)
#index_start_econ <- which(dates_econ == dates_comp[1])

# donde empiezan los datos
econ_start <- as.numeric(c(substr(dates_econ[1],1,4),
                           substr(dates_econ[1],6,7)))
icc_start <- as.numeric(c(substr(dates_icc[1],1,4),
                          substr(dates_icc[1],6,7)))
comp_start <- as.numeric(c(substr(dates_comp[1],1,4),
                           substr(dates_comp[1],6,7)))
var_aux_start <- as.numeric(c(substr(dates_var_aux[1],1,4),
                              substr(dates_var_aux[1],6,7)))
# convertimos a ts
icc_nsa <- ts(icc_nsa, frequency = 12, start = icc_start)
comp_nsa <- ts(comp_nsa, frequency = 12,start = comp_start)
econ_nsa <- ts(econ_nsa, frequency = 12,start = econ_start)

#
# aplicamos ajuste estacional
#

# utilizamos los modelos de INEGI as of 25/11/2019
sa_INEGI <- function(x, models){ # x <- comp_nsa; models <- modelos_comp; i <- 1
  x_sa <- matrix(NA, nrow(x), ncol(x))
  for(i in 1:ncol(x)){ 
    if(models[i,"seasonalma"] == ""){
      x_sa[,i] <- x[,i]
    }else{
    variables <- vector(mode = "character", length = 7)
    variables <- c(td = models[i, "td"], 
                   lpyear = models[i, "lpyear"],
                   easter = models[i, "easter"], 
                   ao = unlist(strsplit(models[i, "AO"], ' ')),
                   tc = unlist(strsplit(  models[i, "TC"], ' ')), 
                   ls = models[i, "LS"],
                   rp = models[i, "RP"]) 
    variables <- variables[variables != ""]
    variables <- variables[!is.na(variables)]
    variables <- paste(variables, collapse = ",")
    
    if(variables == ""){
      x11_adjustment <- seas(x[,i], 
                             arima.model = models[i, "arimamodel"],
                             seats.noadmiss = "no", x11 = "", 
                             outlier = NULL,
                             transform.function = models[i, "transformation"], 
                             regression.aictest = NULL,
                             x11.seasonalma = models[i,"seasonalma"])
      }else{
      x11_adjustment <- seas(x[,i], 
                          regression.variables = variables,
                          arima.model = models[i, "arimamodel"],
                          seats.noadmiss = "no", x11 = "", 
                          outlier = NULL,
                          transform.function = models[i, "transformation"], 
                          regression.aictest = NULL,
                          x11.seasonalma = models[i,"seasonalma"])
    }
    x_sa[,i] <- x11_adjustment$series$d11
    }
  }
  
  x_sa <- ts(x_sa, frequency = 12, start = start(x[,1]))
  colnames(x_sa) <- colnames(x)
  row.names(x_sa) <- row.names(x)
  return(x_sa)
}

icc_sa <- sa_INEGI(icc_nsa, modelos_icc)
colnames(icc_sa) <- colnames(icc_nsa)
row.names(icc_sa) <- dates_icc
comp_sa <- sa_INEGI(comp_nsa, modelos_comp)
colnames(comp_sa) <- colnames(comp_nsa)
row.names(comp_sa) <- dates_comp

# nos quedamos con las variables con a.e.
icc <- icc_sa
comp <- comp_sa

# aplicamos ajuste estacional a las variables Econ que lo requieran
econ_sa <- matrix(NA, nrow(econ_nsa), ncol(econ_nsa))
colnames(econ_sa) <- colnames(econ_nsa)
for(i in 1:ncol(econ_nsa)){
  if(econ_cat[i,"AE"] == 0 ){
    s_start <- min(is.na(econ_nsa[1:(nrow(econ_nsa)-10),i]))
    if(s_start == 0) {
      s_start <- as.numeric(c(substr(dates_econ[1], 1, 4), 
                   substr(dates_econ[1], 6, 7)))
    }else{
      s_start <- as.numeric(c(substr(s_start, 1, 4), 
                              substr(s_start, 6, 7)))
    }
    mod <- seas(ts(econ_nsa[,i], start = s_start,
                   frequency = 12))
    s_sa <- mod$series$s11
    econ_sa[1:length(s_sa),i] <- s_sa
  } else{
    econ_sa[,i] <- ts(econ_nsa[,i], frequency = 12)
  }
}
econ <- ts(econ_sa, start = start(econ_nsa),
           frequency = 12)
row.names(econ) <- dates_econ

# periodo relevante variables econ
econ <- econ[(which(dates_econ == dates_comp[1])):length(dates_econ),]
econ_start <- as.numeric(c(substr(row.names(econ)[1],1,4),
                           substr(row.names(econ)[1],6,7)))
econ <- ts(econ, start = econ_start, frequency = 12)

# convertimos ahorro a variable mensual
# quitamos ahorro ext
sav <- sav[,-which(colnames(sav) == "AHORRO_EXT")]
sav <- sav[(which(row.names(sav) == dates_comp[1]):nrow(sav)),]
sav <- ts(sav, frequency = 4, start = comp_start)
sav_m <- sapply(1:ncol(sav),
                function(x) td(sav[,x] ~ ts(econ[!is.na(econ[,1]),1],
                                            frequency = 12, start = econ_start),
                               to = "monthly", conversion = "average")$values)
colnames(sav_m) <- colnames(sav)
sav_m <- ts(sav_m, frequency = 12, start = comp_start)
sav_m2 <- rbind(sav_m, 
      matrix(NA, dim(econ)[1] - dim(sav_m)[1], ncol(sav_m)))
row.names(sav_m2) <- row.names(econ)

# agregamos ahorro a las variables econOmicas
if(include_sav) 
  econ <- data.frame(econ, sav_m2)
econ <- ts(econ, frequency = 12, start = econ_start)
row.names(econ) <-
  dates_econ[(which(dates_econ == dates_comp[1])):length(dates_econ)]
row.names(econ) <- row.names(sav_m2)

# excluimos icc y lo guardamos en una variable
icc_inegi <- data.frame(icc[,which(colnames(icc) == "ICC")])
icc_inegi <- ts(icc_inegi, frequency = 12, start = icc_start)
row.names(icc_inegi) <- dates_icc
icc <- icc[,which(colnames(icc) != "ICC")]
row.names(icc) <- row.names(icc_inegi)

# nos quedamos con el periodo relevante
icc <- icc[index_start_icc:nrow(icc),]

# construimos la matriz de componentes de icc
icc_mat <- data.frame(cbind(icc, comp))
icc_mat <- ts(icc_mat, frequency = 12, start = comp_start)
row.names(icc_mat) <- row.names(comp)
colnames(icc_mat) <- paste("P",1:15, sep = "")

# aplicamos diferencias logarItmicas !! solo a variables Econ !!!
if(!is.null(d)){
  econ_dl <- apply(econ, 2, function(x) diff(log(x), d))
  if(dlog_icc){
    icc_mat <- apply(icc_mat, 2, function(x) diff(log(x), d))
  }else{
    icc_mat <- icc_mat[-(1:d),]
  }
}
row.names(econ_dl) <- 
  dates_econ[(which(dates_econ == dates_comp[1])+12):length(dates_econ)]
econ_dl <- 
  ts(econ_dl, start = as.numeric(c(substr(row.names(econ_dl)[1],1,4),
               substr(row.names(econ_dl)[1],6,7))),frequency = 12)
row.names(econ_dl) <- 
  dates_econ[(which(dates_econ == dates_comp[1])+12):length(dates_econ)]


sed <- 1-det(cor(econ_dl, use = "complete.obs"))^(1/(ncol(econ_dl)-1))
# including sav   0.8227016 #2020/07/02
# witout sav   0.5026158

# imput NA's with pca
icc_mat[is.na(icc_mat)] <- imputePCA(icc_mat)$fitted[is.na(icc_mat)]
#cor(icc_mat)

dim(econ_dl) #; head(econ_dl); tail(econ_dl)
dim(icc_mat) #; head(icc_mat) ;tail(icc_mat)

#figure1
# mostramos  variables dep
par(mfrow = c(6,5), mai = c(0.3,.3,.3,.3))
for(i in 1:ncol(econ_dl)){
  ts.plot(ts(scale(econ_dl[,i]), frequency = 12, 
             start = start(econ_dl)),
             col = c(4), 
             main = colnames(econ_dl)[i], xlab = "", ylab ="")
  #abline(h = 0)
  segments(2002.5, 0, 2020.5, 0, col = 1, lty = 1)
}

#figure5 
# mostramos  variables indep
par(mfrow = c(4,4), mai = c(0.3,.3,.3,.1))
ts.plot(icc_inegi,# frequency = 12, start = ), 
        main = "ICC", ylab = "", xlab = "", col = 2, ylim = c(0,65) ) # ylim = c(-.5,.5)
#segments(2002.5, 0, 2020.5, 0, col = 2, lty = 1)
for(i in 1:ncol(icc_mat)){
  ts.plot(ts(icc_mat[,i], frequency = 12, 
             start = c(as.numeric(substr(row.names(icc_mat)[1], 1,4)),
                       as.numeric(substr(row.names(econ_dl)[1], 6,7)))),
          col = c(4), main = colnames(icc_mat)[i], xlab = "", 
          ylab ="", ylim = c(0,65) ) # ylim = c(-.5,.5)
  #abline(h = 0)
  lines(icc_inegi, col = 2)
 #segments(2002.5, 0, 2020.5, 0, col = 1, lty = 1)
}

var_aux$short <- 
  rowMeans(cbind(var_aux[,"bonos_0_3"], var_aux[,"bonos_3_5"]), 
           na.rm = TRUE)
var_aux$long <- 
  rowMeans(cbind(var_aux[,"Bonos_10_20"], var_aux[,"bonos_20_30"]),
           na.rm = TRUE)
var_aux$spread <- var_aux[, "long"] - var_aux[,"short"]

# miramos las tasas de interEs
par(mfrow = c(3,1))
ts.plot(ts(var_aux[,c("short","bonos_0_3", "bonos_3_5")],
           frequency = 12, start = var_aux_start), 
        col = c(4,1,2), main = "Short Interest Rates")
ts.plot(ts(var_aux[,c("long","Bonos_10_20", "bonos_20_30")],
           frequency = 12, start = var_aux_start), 
        col = c(4,1,2), main = "Long Interest Rates")
ts.plot(ts(var_aux[,c("spread","short", "long")],
           frequency = 12, start = var_aux_start), 
        col = c(4,1,2), main = "Spread: Long - Short interest rates")


# anAlisis del consumo 
par(mfrow = c(3,2), mai = c(.5,.5,.5,.5))
Acf(econ_dl[, "CONS"])
Acf(econ_dl[, "CONS"], type = "partial")
Acf(econ_dl[, "CONS"], plot = FALSE)
Acf(diff(econ_dl[, "CONS"]))
Acf(diff(econ_dl[, "CONS"]), type = "partial")
Acf(diff(econ_dl[, "CONS"]), plot = FALSE)
Acf(diff(diff(econ_dl[, "CONS"])))
Acf(diff(diff(econ_dl[, "CONS"])), type = "partial")
Acf(diff(diff(econ_dl[, "CONS"])), plot = FALSE)

# modelacion de outliers consumo
mx13 <- seas(ts(econ_dl[-is.na(econ_dl[,"CONS"]),"CONS"], 
                frequency = 12, start = ts_econ_start))

par(mfrow = c(1,1))
plot(mx13)

x13_vars <- mx13$model$regression$variables
 if(any(grep("easter",x13_vars))){
   x13_outliers <- x13_vars[-grep("easter",x13_vars)]
 } else{
   x13_outliers <- x13_vars
}

AO_outliers <- modelos_cons[modelos_cons[,"Indicador"]=="CONS","AO"]
AO_outliers <- strsplit(AO_outliers,",")[[1]]
TC_outliers <- modelos_cons[modelos_cons[,"Indicador"]=="CONS","TC"]
TC_outliers <- strsplit(TC_outliers,",")[[1]]
LS_outliers <- modelos_cons[modelos_cons[,"Indicador"]=="CONS","LS"]
LS_outliers <- strsplit(LS_outliers,",")[[1]]
inegi_outliers <- c(AO_outliers, TC_outliers, LS_outliers)

# outliers solo > periodo relevante
inegi_outliers <- 
  inegi_outliers[as.numeric(substr(inegi_outliers, 3,6)) > as.numeric(substr(my_date_econ_start,1,4))]

month.abb.num <- c("01","02","03","04","05","06",
                   "07","08","09","10", "11","12")
out_year <- substr(inegi_outliers, 3,6)
out_type <- substr(inegi_outliers, 1,2)
out_month <- substr(inegi_outliers, 8,10)
out_month_num <- month.abb.num[sapply(1:length(out_month),
                                      function(x) which(out_month[x] == month.abb))]
out_yr_month <- paste(out_year, out_month_num, sep = "/")

mat_outliers <- matrix(0,nrow = nrow(econ_dl), 
                       ncol = length(inegi_outliers))
row.names(mat_outliers) <- row.names(econ_dl)
colnames(mat_outliers) <- inegi_outliers

for(i in 1:length(inegi_outliers)) {
  if (out_type[i] == "LS") {
    mat_outliers[which(row.names(mat_outliers) == out_yr_month[i]):nrow(mat_outliers), 
                 inegi_outliers[i]] <- 1
  } else if (out_type[i] == "AO") {
    mat_outliers[which(row.names(mat_outliers) == out_yr_month[i]), 
                 inegi_outliers[i]] <- 1
  }
}

#
# analisis de predicciOn 
#


# 
dim(econ_dl)
econ_dl <- 
  econ_dl[(which(row.names(econ_dl) == my_date_econ_start):nrow(econ_dl)),]
dim(econ_dl)

# recortamos variables auxiliares a periodo relevante
var_aux_rel <- 
  var_aux[which(row.names(var_aux) %in% row.names(econ_dl)),]
dim(var_aux_rel)

dim(mat_outliers)
mat_outliers <- 
  mat_outliers[which(row.names(mat_outliers) %in% row.names(econ_dl)),]
dim(mat_outliers)

# nos quedamos con el periodo relevante de las series ICC
icc_mat <- 
  icc_mat[which(row.names(icc_mat) == row.names(econ_dl)[1]):nrow(icc_mat),]
dim(icc_mat) # 184 

icc_inegi_bis <- 
  data.frame(icc_inegi[which(row.names(icc_inegi) == row.names(econ_dl)[1]):nrow(icc_inegi),])
colnames(icc_inegi_bis) <- "ICC"
row.names(icc_inegi_bis) <- 
  dates_icc[which(dates_icc == row.names(econ_dl)[1]):nrow(icc_inegi)]

df_icc_lags <- 
  data.frame(econ_dl[, "CONS"],icc_inegi_bis, 
        y_lag(icc_inegi_bis, 1)[1:(length(y_lag(icc_inegi_bis,1))-1)],
        y_lag(icc_inegi_bis, 2)[1:(length(y_lag(icc_inegi_bis,2))-2)],
        y_lag(icc_inegi_bis, 3)[1:(length(y_lag(icc_inegi_bis,3))-3)],
        y_lag(icc_inegi_bis, 12)[1:(length(y_lag(icc_inegi_bis,12))-12)])
colnames(df_icc_lags) <- c("CONS", "ICC", "ICC_lag1", 
                           "ICC_lag2","ICC_lag3","ICC_lag12")


# funcion para generar modelos del consumo

# forecast horizon
# Iniciamos estimacion de modelos, pronosticos y calculo de error
H <- 24
mat_rmse <- matrix(NA, H, 4)
colnames(mat_rmse) <- c("base_w_ICC","base_w_o_ICC",
                        "regarima_w_ICC", "regarima_w_o_ICC")
summary_list <- list()
k <- 1

for(h in H:1 ){ # h <- 24
  
  # recortamos las variables de acuerdo a h
  # periodo de entrenamiento
  econ_training <- econ_dl[1:(nrow(econ_dl) - h),] # dim(econ_sum_sample)
  econ_training <- ts(econ_training, frequency = 12, start = ts_econ_start)
  row.names(econ_training) <- row.names(econ_dl)[1:(nrow(econ_dl) - h)]
  
  icc_training <- icc_inegi_bis[1:(nrow(icc_inegi_bis) - h),]
  icc_mat_training <- icc_mat[1:(nrow(icc_mat) - h),]
  df_icc_lags_training <- df_icc_lags[1:(nrow(df_icc_lags) - h),]
  var_aux_training <- var_aux_rel[1:(nrow(var_aux_rel) - h),]
  mat_outliers_training <- mat_outliers[1:(nrow(mat_outliers) - h),]

  # series fuera de muestra
  econ_test <- econ_dl[(nrow(econ_dl) - h + 1):nrow(econ_dl),]
  icc_test <- 
    icc_inegi_bis[(nrow(icc_inegi_bis) - h + 1):nrow(icc_inegi_bis),]  
  icc_mat_test <- icc_mat[(nrow(icc_mat) - h + 1):nrow(icc_mat),]
  df_icc_lags_test <- 
    df_icc_lags[(nrow(df_icc_lags) - h+1):nrow(df_icc_lags),]
  var_aux_test <- 
    var_aux_rel[(nrow(var_aux_rel) - h + 1):nrow(var_aux_rel),]
  mat_outliers_test <- 
    mat_outliers[(nrow(mat_outliers) - h + 1):nrow(mat_outliers),]

  # dim(econ_sub_sample); length(icc_sub_sample); dim(icc_mat_sub_sample)

  # Estimacion de Modelos
  ## 1
  xreg <- cbind(ICC = df_icc_lags_training[,"ICC"], 
                mat_outliers_training)
  # Base incluyendo ICC
  base_mod <- auto.arima(ts(econ_training[,"CONS"], frequency = 12,
               start = ts_econ_start), allowdrift = FALSE,
               allowmean = FALSE, xreg = xreg)
  # Base quitando  ICC
  base_mod_noicc <- 
    arima(ts(econ_training[,"CONS"], frequency = 12,start = ts_econ_start),
          order = c(base_mod$arma[1],base_mod$arma[6],base_mod$arma[2]), 
          seasonal = list(order = c(base_mod$arma[3],base_mod$arma[7],
                      base_mod$arma[4]), period = base_mod$arma[5]),
                     xreg = as.matrix(cbind(mat_outliers_training)))
  
  # RegARIMA incluyendo ICC
  regarima_icc <- 
    auto.arima(ts(econ_training[,"CONS"],frequency = 12,
              start = ts_econ_start),
              allowdrift = FALSE, allowmean = FALSE, 
              xreg = ts(cbind(mat_outliers_training,
                    df_icc_lags_training[,"ICC"],
                    var_aux_training[,c("TD_O","INF_ANU","spread")]),
                    frequency = 12, start = ts_econ_start))
  
  # RegARIMA quitando ICC
  regarima_noicc <- 
    arima(ts(econ_training[,"CONS"],frequency = 12,start = ts_econ_start),
          order = c(regarima_icc$arma[1],regarima_icc$arma[6],regarima_icc$arma[2]),
          seasonal = list(order = c(regarima_icc$arma[3],regarima_icc$arma[7],
                      regarima_icc$arma[4]),period = regarima_icc$arma[5]),
          xreg = ts(cbind(mat_outliers_training,
                          var_aux_training[,c("TD_O","INF_ANU","spread")]),
                                      frequency = 12,
                                      start = ts_econ_start) )
  

  # summary(base_mod)
  # summary(base_mod_noicc)
  # summary(naive_icc_lag1)
  # summary(naive_lag1_noicc)
  # summary(naive_icc_lag12)
  # summary(naive_lag12_noicc)
  
  rmse_in_sample <- matrix(NA, 2,2)
  colnames(rmse_in_sample) <- c("w_ICC", "w_o_ICC")
  rmse_in_sample[1,"w_ICC"] <- summary(base_mod)[2]
  rmse_in_sample[1,"w_o_ICC"] <- summary(base_mod_noicc)[2]
  rmse_in_sample[2,"w_ICC"] <- summary(regarima_icc)[2]
  rmse_in_sample[2,"w_o_ICC"] <- summary(regarima_noicc)[2]
  
  mape_in_sample <- matrix(NA, 2,2)
  colnames(mape_in_sample) <- c("w_ICC", "w_o_ICC")
  mape_in_sample[1,"w_ICC"] <- summary(base_mod)[5]
  mape_in_sample[1,"w_o_ICC"] <- summary(base_mod_noicc)[5]
  mape_in_sample[2,"w_ICC"] <- summary(regarima_icc)[5]
  mape_in_sample[2,"w_o_ICC"] <- summary(regarima_noicc)[5]
  
  mase_in_sample <- matrix(NA, 2,2)
  colnames(mase_in_sample) <- c("w_ICC", "w_o_ICC")
  mase_in_sample[1,"w_ICC"] <- summary(base_mod)[6]
  mase_in_sample[1,"w_o_ICC"] <- summary(base_mod_noicc)[6]
  mase_in_sample[2,"w_ICC"] <- summary(regarima_icc)[6]
  mase_in_sample[2,"w_o_ICC"] <- summary(regarima_noicc)[6]
  
  # new_xreg <- cbind(ICC = df_icc_lags_test[,"ICC"],mat_outliers_test)
  # fore_h <- forecast(base_mod, xreg = as.matrix(new_xreg))
  # mat_rmse[k,1] <- 
  #   sqrt(sum((fore_h$mean - econ_test[,"CONS"])^2,
  #            na.rm = TRUE)) / length(fore_h$mean)
  # 
  # new_xreg <- mat_outliers_test
  # predi_h <- predict(base_mod_noicc, newxreg = new_xreg)
  # fore_h <- forecast(base_mod_noicc, xreg = new_xreg)
  # mat_rmse[k,2] <- 
  #   sqrt(sum((fore_h$mean - econ_test[,"CONS"])^2,
  #            na.rm = TRUE)) / length(fore_h$mean)
  # 
  # new_xreg <- cbind(mat_outliers_test,
  #                 ICC = df_icc_lags_test[,"ICC"],
  #                 var_aux_test[,c("TD_O","INF_ANU","spread")])
  # fore_h <- forecast(regarima_icc, xreg = as.matrix(new_xreg))
  # mat_rmse[k,3] <- 
  #   sqrt(sum((fore_h$mean - econ_test[,"CONS"])^2,
  #            na.rm = TRUE)) / length(fore_h$mean)
  # 
  # new_xreg <- cbind(mat_outliers_test,
  #                   var_aux_test[,c("TD_O","INF_ANU","spread")])
  # fore_h <- forecast(regarima_noicc, xreg = as.matrix(new_xreg))
  # mat_rmse[k,4] <- 
  #   sqrt(sum((fore_h$mean - econ_test[,"CONS"])^2,
  #            na.rm = TRUE)) / length(fore_h$mean)
  # 
  summary_list[[k]] <- list(RMSE = rmse_in_sample, 
                            MAPE = mape_in_sample,
                            MASE = mase_in_sample)
  k <- k + 1
  
}

mat_rmse2 <- mat_mape2 <- mat_mase2  <- matrix(0,24,4)
colnames(mat_rmse2) <- colnames(mat_mape2) <- 
  colnames(mat_mase2) <- c("base_icc", "base_no_icc", 
                     "regarima_icc","regarima_icc_noicc")
for(j in 1:length(summary_list)){
  list_temp <- summary_list[[j]]
  mat_rmse2[j, c(1,2)] <- list_temp$RMSE[1,]
  mat_rmse2[j, c(3,4)] <- list_temp$RMSE[2,]
  mat_mape2[j, c(1,2)] <- list_temp$MAPE[1,]
  mat_mape2[j, c(3,4)] <- list_temp$MAPE[2,]
  mat_mase2[j, c(1,2)] <- list_temp$MASE[1,]
  mat_mase2[j, c(3,4)] <- list_temp$MASE[2,]
}

par(mfrow = c(3,1), mai = c(1,1,.5, 3), xpd = TRUE)
ts.plot(mat_rmse2, col = c(4,2,3,1), 
        xlab = "Forecast horizon", ylab = "RMSE")
legend(x = 25,y = .0112,
  legend = c("base_icc", "base_no_icc", 
                     "regarima_icc","regarima_icc_noicc"),
       bty = "n", lty = 1, col = c(4,2,3,1))
ts.plot(mat_mape2, col = c(4,2,3,1), 
        xlab = "Forecast horizon", ylab = "MAPE")
legend(x = 25,y = 57,
       legend = c("base_icc", "base_no_icc", 
                  "regarima_icc","regarima_icc_noicc"),
       bty = "n", lty = 1, col = c(4,2,3,1))
ts.plot(mat_mase2, col = c(4,2,3,1), 
        xlab = "Forecast horizon", ylab = "MASE")
legend(x = 25,y = .9,
       legend = c("base_icc", "base_no_icc", 
                  "regarima_icc","regarima_icc_noicc"),
       bty = "n", lty = 1, col = c(4,2,3,1))




par(mfrow = c(3,1))
barplot(mat_h[,"rmse"], 
        beside = TRUE, main = "RMSE out of sample")
barplot(mat_h[,"mape"], 
        beside = TRUE, main = "MAPE out of sample")
barplot(mat_h[,"mase"],
        beside = TRUE, main = "MASE out of sample")


par(mfrow = c(3,1))
barplot(rmse_naive, beside = TRUE, main = "RMSE",
        legend.text = c("ICC", "ICC-lag1", "ICC-lag12"), 
        args.legend = list(x = 6, y =0.06,horiz = TRUE, 
                           bty = "n"))
barplot(mape_naive, beside = TRUE, main = "MAPE",
        legend.text = c("ICC", "ICC-lag1", "ICC-lag12"), 
        args.legend = list(x = 6, y =0.06,horiz = TRUE, 
                           bty = "n"))
barplot(mase_naive, beside = TRUE, main = "MASE",
        legend.text = c("ICC", "ICC-lag1", "ICC-lag12"), 
        args.legend = list(x = 6, y =0.06,horiz = TRUE, 
                           bty = "n"))



AIC(naive_icc, naive_noicc)
AIC(naive_icc_lag1, naive_lag1_noicc)
AIC(naive_icc_lag12, naive_lag12_noicc)

BIC(naive_icc, naive_noicc)
BIC(naive_icc_lag1, naive_lag1_noicc)
BIC(naive_icc_lag12, naive_lag12_noicc)

par(mfrow = c(3,2))
Acf(resid(naive_icc))
Acf(resid(naive_noicc))
Acf(resid(naive_icc_lag1))
Acf(resid(naive_lag1_noicc))
Acf(resid(naive_icc_lag12))
Acf(resid(naive_lag12_noicc))

par(mfrow = c(3,2))
plot(density(resid(naive_icc)))
plot(density(resid(naive_noicc)))
plot(density(resid(naive_icc_lag1), na.rm = TRUE))
plot(density(resid(naive_lag1_noicc), na.rm = TRUE))
plot(density(resid(naive_icc_lag12), na.rm = TRUE))
plot(density(resid(naive_lag12_noicc), na.rm = TRUE))

Box.test(resid(naive_icc), type = "Ljung-Box")
Box.test(resid(naive_noicc), type = "Ljung-Box")
Box.test(resid(naive_icc_lag1), type = "Ljung-Box")
Box.test(resid(naive_lag1_noicc), type = "Ljung-Box")
Box.test(resid(naive_icc_lag12), type = "Ljung-Box")
Box.test(resid(naive_lag12_noicc), type = "Ljung-Box")

# pronosticos 1

m1 <- naive_icc
m2 <- naive_noicc
m3 <- naive_icc_lag1
m4 <- naive_lag1_noicc
m5 <- naive_icc_lag12
m6 <- naive_lag12_noicc


fcst1 <- predict(m1, h, 
                 newxreg = cbind(df_icc_lags_test[,"ICC"], 
                                 mat_outliers_test))
fcst2 <- predict(m2, h, 
                 newxreg = cbind(mat_outliers_test))
fcst3 <- predict(m3, h, 
                 newxreg = cbind(df_icc_lags_test[,"ICC_lag1"],
                                 mat_outliers_test))
fcst4 <- predict(m4, h, 
                 newxreg = cbind(mat_outliers_test))
fcst5 <- predict(m5, h, 
                 newxreg = cbind(df_icc_lags_test[,"ICC_lag12"],
                                 mat_outliers_test))
fcst6 <- predict(m6, h, 
                 newxreg = cbind(mat_outliers_test))


ts_start <- start(m1$x)

par(mfrow = c(3,2), mai = c(.5, .5, .5,.5))
# m1
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m1),rep(NA,h))), 
           frequency = 12,start = ts_start), 
        ylab = "", xlab = "", col = c(4,2),
        main = "CONS con ICC",
        ylim = c(min(c(fcst1$pred - 1.96*fcst1$se,econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst1$pred + 1.96*fcst1$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst1$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst1$pred - 1.96*fcst1$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst1$pred + 1.96*fcst1$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
# m2
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m2),rep(NA,h))), 
           frequency = 12,start = ts_start), 
        ylab = "", xlab = "", col = c(4, 2),
        main = "CONS sin ICC",
        ylim = c(min(c(fcst2$pred - 1.96*fcst2$se,econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst2$pred + 1.96*fcst2$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst2$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst2$pred - 1.96*fcst2$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst2$pred + 1.96*fcst2$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
# m3
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m3),rep(NA,h))), 
           frequency = 12, start = ts_start), 
        ylab = "", xlab = "", col = c(4,2), 
        main = "CONS con ICC lag1",
        ylim = c(min(c(fcst3$pred - 1.96*fcst3$se, econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst3$pred + 1.96*fcst3$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst3$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst3$pred - 1.96*fcst3$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst3$pred + 1.96*fcst3$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
# m4
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m4),rep(NA,h))), 
           frequency = 12, start = ts_start), 
        ylab = "", xlab = "", col = c(4,2), 
        main = "CONS sin ICC lag 1",
        ylim = c(min(c(fcst4$pred - 1.96*fcst4$se, econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst4$pred + 1.96*fcst4$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst4$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst4$pred - 1.96*fcst4$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst4$pred + 1.96*fcst4$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
# m5
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m5),rep(NA,h))), 
           frequency = 12, start = ts_start), 
        ylab = "", xlab = "", col = c(4,2), 
        main = "naive con ICC lag12",
        ylim = c(min(c(fcst5$pred - 1.96*fcst5$se, econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst5$pred + 1.96*fcst5$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst5$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst5$pred - 1.96*fcst5$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst5$pred + 1.96*fcst5$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
# m6
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m6),rep(NA,h))), 
           frequency = 12, start = ts_start), 
        ylab = "", xlab = "", col = c(4,2), 
        main = "naive lag12 sin ICC",
        ylim = c(min(c(fcst6$pred - 1.96*fcst6$se, econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst6$pred + 1.96*fcst6$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst6$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst6$pred - 1.96*fcst6$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-h), fcst6$pred + 1.96*fcst6$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)



## 2
# full autoarima con variables auxiliares e ICC
naive_plus_icc <- 
  auto.arima(ts(econ_training[,"CONS"],frequency = 12,
      start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
      as.numeric(substr(row.names(econ_training)[1], 6,7))) ),
    xreg = ts(cbind(mat_outliers_training,
       df_icc_lags_training[,"ICC"],
       var_aux_training[,c("TD_O","INF_ANU","spread")]),
      frequency = 12,
       start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
       as.numeric(substr(row.names(econ_training)[1], 6,7))) ) )

# arima variables auxiliares quitando el ICC
naive_plus_noicc <- arima(ts(econ_training[,"CONS"],frequency = 12,
       start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
       as.numeric(substr(row.names(econ_training)[1], 6,7))) ),
      order = c(naive_plus_icc$arma[1],naive_plus_icc$arma[6],
                naive_plus_icc$arma[2]),
        seasonal = list(order = c(naive_plus_icc$arma[3],
                naive_plus_icc$arma[7],naive_plus_icc$arma[4]),
        period = naive_plus_icc$arma[5]),
      xreg = ts(cbind(mat_outliers_training,
        var_aux_training[,c("TD_O","INF_ANU","spread")]),
      frequency = 12,
        start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
        as.numeric(substr(row.names(econ_training)[1], 6,7))) ) )

# autoarima con variables auxiliares e ICC lag1
naive_plus_icc_lag1 <- auto.arima(ts(econ_training[,"CONS"],frequency = 12,
        start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
        as.numeric(substr(row.names(econ_training)[1], 6,7))) ),
      xreg = ts(cbind(mat_outliers_training,
        df_icc_lags_training[,"ICC_lag1"],
        var_aux_training[,c("TD_O","INF_ANU","spread")]),
      frequency = 12,
        start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
        as.numeric(substr(row.names(econ_training)[1], 6,7))) ) )

# autoarima con variables auxiliares e ICC lag1 quitando el ICC
naive_plus_lag1_noicc <- arima(ts(econ_training[,"CONS"],frequency = 12,
        start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
        as.numeric(substr(row.names(econ_training)[1], 6,7))) ),
      order = c(naive_plus_icc_lag1$arma[1],
                naive_plus_icc_lag1$arma[6],naive_plus_icc_lag1$arma[2]),
        seasonal = list(order = c(naive_plus_icc_lag1$arma[3],
                naive_plus_icc_lag1$arma[7],naive_plus_icc_lag1$arma[4]),
                period = naive_plus_icc_lag1$arma[5]),
      xreg = ts(cbind(mat_outliers_training,
        var_aux_training[,c("TD_O","INF_ANU","spread")]),
      frequency = 12,
        start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
        as.numeric(substr(row.names(econ_training)[1], 6,7))) ) )

# arima con variables auxiliares e ICC lag 12
naive_plus_icc_lag12 <- auto.arima(ts(econ_training[,"CONS"],frequency = 12,
        start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
        as.numeric(substr(row.names(econ_training)[1], 6,7))) ),
      xreg = ts(cbind(mat_outliers_training,
        df_icc_lags_training[,"ICC_lag12"],
        var_aux_training[,c("TD_O","INF_ANU","spread")]),
      frequency = 12,
        start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
        as.numeric(substr(row.names(econ_training)[1], 6,7))) ) )

# arima con variables auxiliares e ICC lag 12 quitando ICC
naive_plus_lag12_noicc <- arima(ts(econ_training[,"CONS"],frequency = 12,
        start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
        as.numeric(substr(row.names(econ_training)[1], 6,7))) ),
      order = c(naive_plus_icc_lag12$arma[1],naive_plus_icc_lag12$arma[6],
                naive_plus_icc_lag12$arma[2]),
      seasonal = list(order = c(naive_plus_icc_lag12$arma[3],
                naive_plus_icc_lag12$arma[7],naive_plus_icc_lag12$arma[4]),
        period = naive_plus_icc_lag12$arma[5]),
      xreg = ts(cbind(mat_outliers_training,
        var_aux_training[,c("TD_O","INF_ANU","spread")]),
      frequency = 12,
        start = c(as.numeric(substr(row.names(econ_training)[1],1,4)),
        as.numeric(substr(row.names(econ_training)[1], 6,7))) ) )

summary(naive_plus_icc)
summary(naive_plus_noicc)
summary(naive_plus_icc_lag1)
summary(naive_plus_lag1_noicc)
summary(naive_plus_icc_lag12)
summary(naive_plus_lag12_noicc)

rmse_plus <- matrix(NA, 3,2)
colnames(rmse_plus) <- c("w_ICC", "w_o_ICC")
row.names(rmse_plus) <- c("ICC", "ICC-lag1", "ICC-lag12")
rmse_plus[1,"w_ICC"] <- summary(naive_plus_icc)[2]
rmse_plus[1,"w_o_ICC"] <- summary(naive_plus_noicc)[2]
rmse_plus[2,"w_ICC"] <- summary(naive_plus_icc_lag1)[2]
rmse_plus[2,"w_o_ICC"] <- summary(naive_plus_lag1_noicc)[2]
rmse_plus[3,"w_ICC"] <- summary(naive_plus_icc_lag12)[2]
rmse_plus[3,"w_o_ICC"] <- summary(naive_plus_lag12_noicc)[2]


mape_plus <- matrix(NA, 3,2)
colnames(mape_plus) <- c("w_ICC", "w_o_ICC")
row.names(mape_plus) <- c("ICC", "ICC-lag1", "ICC-lag12")
mape_plus[1,"w_ICC"] <- summary(naive_plus_icc)[5]
mape_plus[1,"w_o_ICC"] <- summary(naive_plus_noicc)[5]
mape_plus[2,"w_ICC"] <- summary(naive_plus_icc_lag1)[5]
mape_plus[2,"w_o_ICC"] <- summary(naive_plus_lag1_noicc)[5]
mape_plus[3,"w_ICC"] <- summary(naive_plus_icc_lag12)[5]
mape_plus[3,"w_o_ICC"] <- summary(naive_plus_lag12_noicc)[5]

mase_plus <- matrix(NA, 3,2)
colnames(mase_plus) <- c("w_ICC", "w_o_ICC")
row.names(mase_plus) <- c("ICC", "ICC-lag1", "ICC-lag12")
mase_plus[1,"w_ICC"] <- summary(naive_plus_icc)[6]
mase_plus[1,"w_o_ICC"] <- summary(naive_plus_noicc)[6]
mase_plus[2,"w_ICC"] <- summary(naive_plus_icc_lag1)[6]
mase_plus[2,"w_o_ICC"] <- summary(naive_plus_lag1_noicc)[6]
mase_plus[3,"w_ICC"] <- summary(naive_plus_icc_lag12)[6]
mase_plus[3,"w_o_ICC"] <- summary(naive_plus_lag12_noicc)[6]

par(mfrow = c(3,1))
barplot(rmse_plus, beside = TRUE, main = "RMSE", 
        legend.text = c("ICC", "ICC-lag1", "ICC-lag12"), 
        args.legend = list(x = 6, y =0.06,horiz = TRUE, 
                           bty = "n"))
barplot(mape_plus, beside = TRUE, main = "MAPE", 
        legend.text = c("ICC", "ICC-lag1", "ICC-lag12"), 
        args.legend = list(x = 6, y =0.06,horiz = TRUE, 
                           bty = "n"))
barplot(mase_plus, beside = TRUE, main = "MASE", 
        legend.text = c("ICC", "ICC-lag1", "ICC-lag12"), 
        args.legend = list(x = 6, y =0.06,horiz = TRUE, 
                           bty = "n"))


AIC(naive_plus_icc, naive_plus_noicc)
AIC(naive_plus_icc_lag1, naive_plus_lag1_noicc)
AIC(naive_plus_icc_lag12, naive_plus_lag12_noicc)

BIC(naive_plus_icc, naive_plus_noicc)
BIC(naive_plus_icc_lag1, naive_plus_lag1_noicc)
BIC(naive_plus_icc_lag12, naive_plus_lag12_noicc)


par(mfrow = c(3,2))
Acf(resid(naive_plus_icc))
Acf(resid(naive_plus_noicc))
Acf(resid(naive_plus_icc_lag1))
Acf(resid(naive_plus_lag1_noicc))
Acf(resid(naive_plus_icc_lag12))
Acf(resid(naive_plus_lag12_noicc))

par(mfrow = c(3,2))
plot(density(resid(naive_plus_icc), na.rm = TRUE))
plot(density(resid(naive_plus_noicc), na.rm = TRUE))
plot(density(resid(naive_plus_icc_lag1), na.rm = TRUE))
plot(density(resid(naive_plus_lag1_noicc), na.rm = TRUE))
plot(density(resid(naive_plus_icc_lag12), na.rm = TRUE))
plot(density(resid(naive_plus_lag12_noicc), na.rm = TRUE))

Box.test(resid(naive_plus_icc), type = "Ljung-Box")
Box.test(resid(naive_plus_noicc), type = "Ljung-Box")
Box.test(resid(naive_plus_icc_lag1), type = "Ljung-Box")
Box.test(resid(naive_plus_lag1_noicc), type = "Ljung-Box")
Box.test(resid(naive_plus_icc_lag12), type = "Ljung-Box")
Box.test(resid(naive_plus_lag12_noicc), type = "Ljung-Box")


# pronosticos 2

m1 <- naive_plus_icc
m2 <- naive_plus_noicc
m3 <- naive_plus_icc_lag1
m4 <- naive_plus_lag1_noicc
m5 <- naive_plus_icc_lag12
m6 <- naive_plus_lag12_noicc


fcst1 <- predict(m1, H, 
           newxreg = cbind(df_icc_lags_test[,"ICC"], 
                           mat_outliers_test))
fcst2 <- predict(m2, H, 
           newxreg = cbind(mat_outliers_test))
fcst3 <- predict(m3, H, 
          newxreg = cbind(df_icc_lags_test[,"ICC_lag1"],
                          mat_outliers_test))
fcst4 <- predict(m4, H, 
          newxreg = cbind(mat_outliers_test))
fcst5 <- predict(m5, H, 
          newxreg = cbind(df_icc_lags_test[,"ICC_lag12"],
                                 mat_outliers_test))
fcst6 <- predict(m6, H, 
          newxreg = cbind(mat_outliers_test))


ts_start <- start(m1$x)
  
par(mfrow = c(3,2), mai = c(.5, .5, .5,.5))
# m1
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m1),rep(NA,H))), 
           frequency = 12,start = ts_start), 
        ylab = "", xlab = "", col = c(4,2),
        main = "naive + TD + INF + SPREAD  + ICC",
        ylim = c(min(c(fcst1$pred - 1.96*fcst1$se,econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst1$pred + 1.96*fcst1$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst1$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst1$pred - 1.96*fcst1$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst1$pred + 1.96*fcst1$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
# m2
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m2),rep(NA,H))), 
           frequency = 12,start = ts_start), 
        ylab = "", xlab = "", col = c(4, 2),
        main = "naive + TD + INF + SPREAD ",
        ylim = c(min(c(fcst2$pred - 1.96*fcst2$se,econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst2$pred + 1.96*fcst2$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst2$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst2$pred - 1.96*fcst2$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst2$pred + 1.96*fcst2$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
# m3
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m3),rep(NA,H))), 
           frequency = 12, start = ts_start), 
        ylab = "", xlab = "", col = c(4,2), 
        main = "naive + TD + INF + SPREAD  + ICC lag 1",
        ylim = c(min(c(fcst3$pred - 1.96*fcst3$se, econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst3$pred + 1.96*fcst3$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst3$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst3$pred - 1.96*fcst3$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst3$pred + 1.96*fcst3$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
# m4
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m4),rep(NA,H))), 
           frequency = 12, start = ts_start), 
        ylab = "", xlab = "", col = c(4,2), 
        main = "naive + TD + INF + SPREAD , sin ICC lag1 ",
       ylim = c(min(c(fcst4$pred - 1.96*fcst4$se, econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst4$pred + 1.96*fcst4$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst4$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst4$pred - 1.96*fcst4$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst4$pred + 1.96*fcst4$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
# m5
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m5),rep(NA,H))), 
           frequency = 12, start = ts_start), 
        ylab = "", xlab = "", col = c(4,2),
        main = "naive + TD + INF + SPREAD  + ICC lag12",
        ylim = c(min(c(fcst5$pred - 1.96*fcst5$se, econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst5$pred + 1.96*fcst5$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst5$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst5$pred - 1.96*fcst5$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst5$pred + 1.96*fcst5$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
# m6
ts.plot(ts(cbind(econ_dl[,"CONS"],c(fitted(m6),rep(NA,H))), 
           frequency = 12, start = ts_start), 
        ylab = "", xlab = "", col = c(4,2), 
        main = "naive + TD + INF + SPREAD, sin ICC lag12",
        ylim = c(min(c(fcst6$pred - 1.96*fcst6$se, econ_dl[,"CONS"]), na.rm = TRUE),
                 max(fcst6$pred + 1.96*fcst6$se)))
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst6$pred),
         frequency = 12, start = ts_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst6$pred - 1.96*fcst6$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-H), fcst6$pred + 1.96*fcst6$se),
         frequency = 12, start = ts_start), col = "lightblue", lwd = 2)



#
## Modelos CC y PLS
#


# number of variables
N <- ncol(icc_mat)

# matrix of models 
n_comb <- 0
for(i in 1 : N)
  n_comb <- n_comb + ncol(combn(N, i))

mat_models <- matrix(FALSE, n_comb, N)

h <- 0
for(i in 1 : N){
  combn_i <- combn(N, i)

  for(j in 1 : ncol(combn_i)){
    J <- j + h
    mat_models[J, combn_i[,j]] <- TRUE
  }
  h <- h + ncol(combn_i)
}

# matrix of correlations
rho_mat_pls <- rho_mat_cc <- matrix(NA, nrow(mat_models), ncol(econ_dl))
colnames(rho_mat_pls) <- colnames(rho_mat_cc) <- 
  paste("rho_", colnames(econ_dl), sep = "")

#econ <- econ[index_start_econ:nrow(econ),]
matcor(econ_dl, icc_mat[1:nrow(econ_dl), ])



# estimate CC
for(i in 1 : nrow(mat_models)){ print(i)
  
  icc_mat_i <- 
    as.matrix(icc_mat[1:nrow(econ_dl),mat_models[i,],drop = FALSE])
  
  print(paste("n = ", nrow(econ_dl), "; p = ", ncol(econ_dl), 
              "; q = ", ncol(icc_mat_i), sep = "") )
  
  cc_i <- cc(scale(econ_dl), scale(icc_mat_i[,,drop=FALSE ]))
  ft_cc <- -cc_i$scores$yscores[,1]
  P_cc <- cc_i$scores$corr.Y.yscores[,1]
  
  ft_cc_starts <- min(as.numeric(names(ft_cc)))
  ft_cc_ends <- max(as.numeric(row.names(ft_cc)))
  
  if(ncol(icc_mat_i) == 1){
    rho_mat_cc[i,] <- cor(econ_dl, ft_cc, use = "complete.obs")
  }else{
    if(all(P_cc < 0)){
      rho_mat_cc[i,] <- cor(econ_dl, -ft_cc, use = "complete.obs")
    }
    if(all(P_cc > 0)){
      rho_mat_cc[i,] <- cor(econ_dl, ft_cc, use = "complete.obs")
    }
  }
}

# estimate PLS
for(i in 1 : nrow(mat_models)){ print(i)
  
  icc_mat_i <- 
    as.matrix(icc_mat[1:nrow(econ_dl),mat_models[i,],drop=FALSE])
  
  pls_regre <- plsr(scale(econ_dl) ~ scale(icc_mat_i[,,drop=FALSE]))
  ft_pls <- -scores(pls_regre)[, 1, drop=FALSE]
  P_pls <- pls_regre$loadings[, 1, drop=FALSE]
  
  ft_pls_starts <- 
    which(row.names(econ_dl) == row.names(ft_pls)[1])
  ft_pls_ends <- 
    which(row.names(econ_dl) == row.names(ft_pls)[length(ft_pls)])
  
  if(ncol(icc_mat_i) == 1){
    rho_mat_pls[i,] <- cor(econ_dl[ft_pls_starts:ft_pls_ends,],
                           ft_pls, use = "complete.obs")
  }else{
    if(all(P_pls < 0)){
      rho_mat_pls[i,] <- cor(econ_dl[ft_pls_starts:ft_pls_ends,],
                             -ft_pls, use = "complete.obs")
    }
    if(all(P_pls > 0)){
      rho_mat_pls[i,] <- cor(econ_dl[ft_pls_starts:ft_pls_ends,],
                             ft_pls, use = "complete.obs")
    }
  }
}

save.image(paste(path, "icc-savf-iccDL.RData"))


