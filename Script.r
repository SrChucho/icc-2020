
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
library(imputeTS)
library(lmtest)

#library(fable) check 


# path
path <- "C:/Users/jesus.lopezp/Desktop/ICC/ICC_2020/"
source(paste(path, "Functions.R", sep = ""))


# difference
d <- 12L
# Seleccionar variable a pronosticar
# CONS / BYS_NAL  / B_NAL  / SERV_NAL / CONS_M
Y <- "CONS" 
sel_var_aux <- c("TD_SA","INF_ANU","short", 
                 "ISBS_MAY", "ISBS_MEN", "IPM")
include_sav <- FALSE 
my_date_econ_start <- "2005/01"

get_date <- function(x){
  as.numeric(c(substr(x,1,4),
               substr(x,6,7)))
}

# data
icc_nsa <- read.csv(paste(path, "inputs/ICC_NSA.csv", sep = ""), row.names = 1)
icc_sa <- read.csv(paste(path, "inputs/ICC_SA.csv", sep = ""), row.names = 1)
comp_nsa <- read.csv(paste(path, "inputs/COMP_NSA.csv", sep = ""), row.names = 1)
comp_sa <- read.csv(paste(path, "inputs/COMP_SA.csv", sep = ""), row.names = 1)
etco_icc <- read.csv(paste(path, "inputs/ETCO_ICC.csv", sep = ""), row.names = 1)
etco_comp <- read.csv(paste(path, "inputs/ETCO_COMP.csv", sep = ""), row.names = 1)

# leemos variables del consumo
econ_nsa <- read.csv(paste(path, "inputs/Econ4.csv", sep = ""), row.names = 1)
econ_cat <- read.csv(paste(path, "inputs/econ_cat2.csv", sep = ""))

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
var_aux <- read.csv(paste(path, "inputs/var_aux2.csv", sep = ""), 
                    row.names = 1)

# leemos datos de google trends
bcd <- read.csv(paste(path, "inputs/gtrends/BCD3.csv", sep = ""), 
                row.names = 1)
non_durables <- read.csv(paste(path, "inputs/gtrends/non-durables3.csv", sep = ""), 
                         row.names = 1)
services <- read.csv(paste(path, "inputs/gtrends/services3.csv", sep = ""), 
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
econ_start <- get_date(dates_econ[1])
icc_start <- get_date(dates_icc[1])
comp_start <- get_date(dates_comp[1])
var_aux_start <- get_date(dates_var_aux[1])
  
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
      s_start <- get_date(dates_econ[1])
        
    }else{
      s_start <- get_date(s_start)
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


# ajuste estacional Var_Aux, 
# solo para TD_O, las demAS variables no requieren a.e.
td_o <- ts(var_aux[(max(which(is.na(var_aux[-nrow(var_aux),"TD_O"])))+1):(nrow(var_aux)-1),"TD_O"], 
           start = get_date(dates_var_aux[max(which(is.na(var_aux[-nrow(var_aux),"TD_O"])))]), 
           frequency = 12)

var_aux$TD_SA <- NA
var_aux[(max(which(is.na(var_aux[-nrow(var_aux),"TD_O"])))+1):(nrow(var_aux)-1),"TD_SA"] <- 
  series(seas(td_o),"s11")

# imputamos NA's con NOCB     
var_aux <- na.locf(var_aux, fromLast = TRUE, na.rm = FALSE)
row.names(var_aux) <- dates_var_aux

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
# imputamos NA's con NOCB
icc_mat <- na.locf(icc_mat, fromLast = TRUE, na.rm = FALSE)
icc_mat <- ts(icc_mat, frequency = 12, start = comp_start)
row.names(icc_mat) <- row.names(comp)
colnames(icc_mat) <- paste("P",1:15, sep = "")

# aplicamos diferencias logarItmicas 
if(!is.null(d)){
  econ_dl <- apply(econ, 2, function(x) diff(log(x), d))
}

# imputamos NA's con NOCB
econ_dl <- na.locf(econ_dl, fromLast = TRUE, na.rm = FALSE)
row.names(econ_dl) <- 
  dates_econ[(which(dates_econ == dates_comp[1])+12):length(dates_econ)]
econ_dl <- 
  ts(econ_dl, start = get_date(row.names(econ_dl)[1]),frequency = 12)
row.names(econ_dl) <- 
  dates_econ[(which(dates_econ == dates_comp[1])+12):length(dates_econ)]


#sed <- 1-det(cor(econ_dl, use = "complete.obs"))^(1/(ncol(econ_dl)-1))
# including sav   0.8227016 #2020/07/02
# including sav   0.8693612 #2020/07/14
# witout sav   0.5026158


#cor(icc_mat)

dim(econ_dl) #; head(econ_dl); tail(econ_dl)
dim(icc_mat) #; head(icc_mat) ;tail(icc_mat)

#figure1
# mostramos  variables dep
par(mfrow = c(2,3), mai = c(0.3,.3,.3,.3))
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
             start = get_date(row.names(icc_mat)[1])),
          col = c(4), main = colnames(icc_mat)[i], xlab = "", 
          ylab ="", ylim = c(0,65) ) # ylim = c(-.5,.5)
  #abline(h = 0)
  lines(icc_inegi, col = 2)
 #segments(2002.5, 0, 2020.5, 0, col = 1, lty = 1)
}

# imputamos var_aux
var_aux_names <- colnames(var_aux)
var_aux <- sapply(1:ncol(var_aux), function(x) na.interp(var_aux[,x]))
row.names(var_aux) <- dates_var_aux
colnames(var_aux) <- var_aux_names

var_aux <- data.frame(var_aux)
var_aux$short <- 
  rowMeans(cbind(var_aux[,"bonos_3"], var_aux[,"bonos_5"]), 
           na.rm = TRUE)
var_aux$long <- 
  rowMeans(cbind(var_aux[,"Bonos_10"], var_aux[,"bonos_20"],
                 var_aux[,"bonos_30"]),
           na.rm = TRUE)
var_aux$spread <- var_aux[, "long"] - var_aux[,"short"]

# miramos las tasas de interEs
par(mfrow = c(3,1))
ts.plot(ts(var_aux[,c("short","bonos_3", "bonos_5")],
           frequency = 12, start = var_aux_start), 
        col = c(3,1,1), main = "Short Interest Rates", lwd = c(2,1,1))
ts.plot(ts(var_aux[,c("long","Bonos_10", "bonos_20", "bonos_30")],
           frequency = 12, start = var_aux_start), 
        col = c(4,1,1,1), main = "Long Interest Rates", lwd = c(2,1,1,1))
ts.plot(ts(var_aux[,c("spread","short", "long")],
           frequency = 12, start = var_aux_start), 
        col = c(2,3,4), main = "Spread: Long - Short interest rates")

# miramos var aux
par(mfrow = c(3,3))
for(i in 1:length(sel_var_aux)){
ts.plot(ts(var_aux[,sel_var_aux[i]],
           frequency = 12, start = var_aux_start),
  main = sel_var_aux[i], ylab = "", xlab = "")
}


# anAlisis de Y 
par(mfrow = c(3,2), mai = c(.5,.5,.5,.5))
Acf(econ_dl[, Y])
Acf(econ_dl[, Y], type = "partial")
Acf(econ_dl[, Y], plot = FALSE)
Acf(diff(econ_dl[, Y]))
Acf(diff(econ_dl[, Y]), type = "partial")
Acf(diff(econ_dl[, Y]), plot = FALSE)
Acf(diff(diff(econ_dl[, Y])))
Acf(diff(diff(econ_dl[, Y])), type = "partial")
Acf(diff(diff(econ_dl[, Y])), plot = FALSE)


# figure 11
par(mfrow = c(1,1))
ts.plot(ts(econ_dl[,Y], 
   frequency = 12, start = c(2004,1)), 
   col = "red", lwd = 2, ylab = "")
   
# modelacion de outliers consumo
mx13 <- seas(ts(econ_dl[,Y], 
                frequency = 12, start = c(2004,1)))
plot(mx13)
cons_s11 <- series(mx13, "s11")

x13_vars <- mx13$model$regression$variables
 if(any(grep("easter",x13_vars))){
   x13_outliers <- x13_vars[-grep("easter",x13_vars)]
 } else{
   x13_outliers <- x13_vars
}

AO_outliers <- modelos_cons[modelos_cons[,"Indicador"]==Y,"AO"]
AO_outliers <- strsplit(AO_outliers,",")[[1]]
TC_outliers <- modelos_cons[modelos_cons[,"Indicador"]==Y,"TC"]
TC_outliers <- strsplit(TC_outliers,",")[[1]]
LS_outliers <- modelos_cons[modelos_cons[,"Indicador"]==Y,"LS"]
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

# Google Trends section
gt <- data.frame(matrix(NA, dim(bcd),4))
row.names(gt) <- row.names(bcd)
colnames(gt) <- c("consumo_comp","bcd_comp", 
                        "nd_comp", "serv_comp")

gt_start <- get_date(row.names(gt)[1])
gt_dates <- paste(substr(row.names(gt), 1, 4), "/",
                  substr(row.names(gt), 6, 7),sep = "") 

bcd <- bcd[,!is.na(apply(bcd,2,var))]
non_durables <- non_durables[,!is.na(apply(non_durables,2,var))]
services <- services[,!is.na(apply(services,2,var))]

gt[,"bcd_comp"] <- princomp(bcd)$scores[,"Comp.1"]
gt[,"nd_comp"] <- -princomp(non_durables)$scores[,"Comp.1"]
gt[,"serv_comp"] <- princomp(services)$scores[,"Comp.1"]
gt[,"consumo_comp"] <- rowMeans(gt[,-1])


# google trends a.e.
gt_seas <- data.frame(matrix(NA, nrow = nrow(gt),
                             ncol = ncol(gt)))
row.names(gt_seas) <- gt_dates
colnames(gt_seas) <- colnames(gt)

for(i in 2:ncol(gt)){
   x13_temp <- seas(ts(gt[,i], frequency = 12, start = gt_start))
   gt_seas[,i] <- series(x13_temp, "s11")
}
gt_seas[,"consumo_comp"] <- rowMeans(gt_seas[,-1])

# plot google trends
ts.plot(ts(gt_seas, start = gt_start, frequency = 12), 
        col = c(4, 2, 3, 1), lwd = c(2, 1,1,1) )
legend("topleft", col = c(4,2,3,1), lty = 1, 
       horiz = TRUE, legend = colnames(gt_seas), cex = 0.8,
       lwd = c(2,1,1,1), bty = "n")


#
# analisis de predicciOn 
#


# 
dim(econ_dl)
econ_dl <- 
  econ_dl[(which(row.names(econ_dl) == my_date_econ_start):nrow(econ_dl)),]
dim(econ_dl)
econ_dl_start <- get_date(row.names(econ_dl)[1])
econ_dl_row_names <- row.names(econ_dl)
econ_dl <- ts(econ_dl, frequency = 12, start = econ_dl_start)
row.names(econ_dl) <- econ_dl_row_names

# icc
icc_inegi_bis <- 
  icc_inegi[which(row.names(icc_inegi) %in% row.names(econ_dl)),]

#gt
gt_seas <- gt_seas[which(row.names(gt_seas) %in% row.names(econ_dl)),]
dim(gt_seas)


# data frame con lags a probar
cons_lags <- 
  data.frame(econ_dl[-(1:12), Y], 
             y_lag(econ_dl[, Y], 1)[13:(length(y_lag(econ_dl[, Y],1))-1)],
             y_lag(econ_dl[, Y], 2)[13:(length(y_lag(econ_dl[, Y],2))-2)],
             y_lag(econ_dl[, Y], 3)[13:(length(y_lag(econ_dl[, Y],3))-3)],
             y_lag(econ_dl[, Y], 12)[13:(length(y_lag(econ_dl[, Y],12))-12)])

colnames(cons_lags) <- c("Y", "Y_lag1", 
                        "Y_lag2","Y_lag3","Y_lag12")


icc_lags <- 
  data.frame(icc_inegi_bis[-c(1:12)], 
        y_lag(icc_inegi_bis, 1)[13:(length(y_lag(icc_inegi_bis,1))-1)],
        y_lag(icc_inegi_bis, 2)[13:(length(y_lag(icc_inegi_bis,2))-2)],
        y_lag(icc_inegi_bis, 3)[13:(length(y_lag(icc_inegi_bis,3))-3)],
        y_lag(icc_inegi_bis, 12)[13:(length(y_lag(icc_inegi_bis,12))-12)])
        
colnames(icc_lags) <- c("ICC", "ICC_lag1", 
                        "ICC_lag2","ICC_lag3","ICC_lag12")
gt_lags <- 
  data.frame(gt_seas[-(1:12),-1], 
        sapply(2:4,function(x) y_lag(gt_seas[,x],1))[13:(length(y_lag(gt_seas[,1],1))-1),],
        sapply(2:4,function(x) y_lag(gt_seas[,x],2))[13:(length(y_lag(gt_seas[,1],2))-2),],
        sapply(2:4,function(x) y_lag(gt_seas[,x],3))[13:(length(y_lag(gt_seas[,1],3))-3),],
        sapply(2:4,function(x) y_lag(gt_seas[,x],12))[13:(length(y_lag(gt_seas[,1],12))-12),])

colnames(gt_lags) <- c("bcd","nd","serv",
                       "bcd_lag1","nd_lag1","serv_lag1",
                       "bcd_lag2","nd_lag2","serv_lag2",
                       "bcd_lag3","nd_lag3","serv_lag3",
                       "bcd_lag12","nd_lag12","serv_lag12")
df_lags <- data.frame(cons_lags, icc_lags, gt_lags)
colnames(df_lags) <- c(colnames(cons_lags), colnames(icc_lags), 
                       colnames(gt_lags))
                           
ts_icc_lags_start <- get_date(row.names(df_lags)[1])

# nos quedamos con el periodo relevante de las series outliers
dim(mat_outliers)
mat_outliers <- 
  mat_outliers[which(row.names(mat_outliers) %in% row.names(df_lags)),]
dim(mat_outliers)

# nos quedamos con el periodo relevante de las series ICC
icc_mat <- 
  icc_mat[which(row.names(icc_mat) %in% row.names(df_lags)),]
dim(icc_mat) # 184 


# recortamos variables auxiliares a periodo relevante
var_aux_rel <- 
  var_aux[which(row.names(var_aux) %in% row.names(df_lags)),]
dim(var_aux_rel)

# parametros del modelo
Ty <- length(df_lags[!is.na(df_lags[,Y]),Y]) 


# funcion para generar modelos del consumo

# forecast horizon
# Iniciamos estimacion de modelos, pronosticos y calculo de error
Ht <- 24
mat_rmse <- matrix(NA, Ht, 6)
colnames(mat_rmse) <- c("base_w_ICC","base_w_o_ICC",
                        "regarima_w_ICC", "regarima_w_o_ICC",
                        "regarima_w_gt", "regarima_w_o_gt")
row.names(mat_rmse) <- row.names(econ_dl)[(nrow(econ_dl)-Ht + 1):nrow(econ_dl)]
summary_list <- list()

training_forecasts <- training_lower <- training_upper <- 
  matrix(NA, Ht, 6)
colnames(training_forecasts) <- colnames(training_lower) <- 
  colnames(training_upper) <- c("base_w_ICC","base_w_o_ICC",
                                "regarima_w_ICC", "regarima_w_o_ICC",
                                "regarima_w_gt", "regarima_w_o_gt")

k <- 1
h <- 1

coefs_base <- coefs_base_bis <- 
  coefs_reg_icc <- coefs_reg_icc_bis <- 
  coefs_reg_gt <- coefs_reg_gt_bis <- list()
  

icc_lags_list <- list("ICC_lag1",
                  c("ICC_lag1", "ICC_lag2"),
                  c("ICC_lag1", "ICC_lag2", "ICC_lag3"),
                  c("ICC_lag1", "ICC_lag2", "ICC_lag3", "ICC_lag12"))

y_lags_list <- list("Y_lag1",
                    c("Y_lag1", "Y_lag2"),
                    c("Y_lag1", "Y_lag2", "Y_lag3"),
                    c("Y_lag1", "Y_lag2", "Y_lag3", "Y_lag12"))


Sys.time()
for(h in 1:Ht){ print(h) 
  # h <- 1; print(h)
  
  # dato a pronosticar 
  #econ_dl[nrow(econ_dl),Y] <- NA
  
  # Estimacion de Modelos
  ## 1
  # Base incluyendo ICC
  check_outliers <- sapply(1:ncol(mat_outliers), 
                           function(x) sd(mat_outliers[1:(Ty - (Ht - h + 1)),x]) == 0 )
  
  curr_mat_outliers <- mat_outliers[,!check_outliers]
  
  i <- 1
  df_model <- data.frame(cbind(df_lags[1:(Ty - (Ht - h + 1)), "Y"] ,
    df_lags[1:(Ty - (Ht - h + 1)), as.character(y_lags_list[[i]])],
    df_lags[1:(Ty - (Ht - h + 1)), "ICC"],
    df_lags[1:(Ty - (Ht - h + 1)), as.character(icc_lags_list[[i]])],
    curr_mat_outliers[1:(Ty - (Ht - h + 1)), ]))
  
  colnames(df_model) <- c("Y", as.character(y_lags_list[[i]]), 
                          "ICC", as.character(icc_lags_list[[i]]),
                          colnames(curr_mat_outliers))
  
  big_mod <- lm(Y ~ ., data = df_model)
  
  mod_base_list <- list()
  for(i in 1:length(lags_list)){
    mod_base_list[[i]] <-
      auto.arima(ts(df_lags[1:(Ty - (Ht - h + 1)), Y],
                    frequency = 12, start = ts_icc_lags_start),
                 allowdrift = FALSE,d = 0, max.q = 0, max.Q = 0, max.P = 0, 
                 xreg = as.matrix(cbind(df_lags[1:(Ty - (Ht - h + 1)),as.character(lags_list[[i]])],
                                        curr_mat_outliers[1:(Ty - (Ht - h + 1)), ]))
      )
  }
  
  
  
  model_sel <- 
    mod_base_list[[which.min(sapply(mod_base_list, BIC))]]
  lag_sel <- lags_list[[which.min(sapply(mod_base_list, BIC))]]
  print(paste("base_mod best lag: ",lag_sel, 
              "-", model_sel, sep = ""))
  
  
  p <- model_sel$arma[1]
  d <- model_sel$arma[6]
  q <- model_sel$arma[2]
  P <- model_sel$arma[3]
  D <- model_sel$arma[7]
  Q <- model_sel$arma[4]
  
  coefs_base[[h]] <- coeftest(model_sel)  
  
  # forecast
  n <- length(lag_sel) + ncol(curr_mat_outliers)
  new_xreg <- matrix(NA, 1, n)
  new_xreg[, 1:(1+length(lag_sel))] <- 
    df_lags[(Ty - (Ht - h)):(Ty - (Ht - h)),as.character(lag_sel)]
  new_xreg[, -1] <- curr_mat_outliers[(Ty - (Ht - h)):(Ty - (Ht - h)),]
  colnames(new_xreg) <- c(as.character(lag_sel), colnames(curr_mat_outliers))
  
  fore_h <- forecast(model_sel, xreg = new_xreg)
  mat_rmse[k,"base_w_ICC"] <- 
    sqrt((fore_h$mean - 
            econ_dl[(Ty - (Ht - h )):(Ty - (Ht - h)),Y])^2)
  
  training_forecasts[h, "base_w_ICC"] <- fore_h$mean
  training_lower[h, "base_w_ICC"] <- fore_h$lower[,"95%"]
  training_upper[h, "base_w_ICC"] <- fore_h$upper[,"95%"]
}


Sys.time()
for(h in 1:Ht){ print(h) 
  # h <- 1; print(h)
  
  # dato a pronosticar 
  #econ_dl[nrow(econ_dl),Y] <- NA
  
  # Estimacion de Modelos
  ## 1
    # Base incluyendo ICC
    check_outliers <- sapply(1:ncol(mat_outliers), 
         function(x) sd(mat_outliers[1:(Ty - (Ht - h + 1)),x]) == 0 )
  
    curr_mat_outliers <- mat_outliers[,!check_outliers]
    
    lags_list <- list(c("ICC"),
         c("ICC", "ICC_lag1"),
         c("ICC", "ICC_lag1", "ICC_lag2"),
         c("ICC", "ICC_lag1", "ICC_lag2", "ICC_lag3"),
         c("ICC", "ICC_lag1", "ICC_lag2", "ICC_lag3", "ICC_lag12"))
    
    mod_base_list <- list()
    for(i in 1:length(lags_list)){
      mod_base_list[[i]] <-
        auto.arima(ts(df_lags[1:(Ty - (Ht - h + 1)), Y],
             frequency = 12, start = ts_icc_lags_start),
          allowdrift = FALSE,d = 0, max.q = 0, max.Q = 0, max.P = 0, 
          xreg = as.matrix(cbind(df_lags[1:(Ty - (Ht - h + 1)),as.character(lags_list[[i]])],
                       curr_mat_outliers[1:(Ty - (Ht - h + 1)), ]))
          )
    }
    
    
    
    model_sel <- 
      mod_base_list[[which.min(sapply(mod_base_list, BIC))]]
    lag_sel <- lags_list[[which.min(sapply(mod_base_list, BIC))]]
    print(paste("base_mod best lag: ",lag_sel, 
                "-", model_sel, sep = ""))
    
    
    p <- model_sel$arma[1]
    d <- model_sel$arma[6]
    q <- model_sel$arma[2]
    P <- model_sel$arma[3]
    D <- model_sel$arma[7]
    Q <- model_sel$arma[4]

    coefs_base[[h]] <- coeftest(model_sel)  

    # forecast
    n <- length(lag_sel) + ncol(curr_mat_outliers)
    new_xreg <- matrix(NA, 1, n)
    new_xreg[, 1:(1+length(lag_sel))] <- 
      df_lags[(Ty - (Ht - h)):(Ty - (Ht - h)),as.character(lag_sel)]
    new_xreg[, -1] <- curr_mat_outliers[(Ty - (Ht - h)):(Ty - (Ht - h)),]
    colnames(new_xreg) <- c(as.character(lag_sel), colnames(curr_mat_outliers))
    
    fore_h <- forecast(model_sel, xreg = new_xreg)
    mat_rmse[k,"base_w_ICC"] <- 
      sqrt((fore_h$mean - 
              econ_dl[(Ty - (Ht - h )):(Ty - (Ht - h)),Y])^2)
    
    training_forecasts[h, "base_w_ICC"] <- fore_h$mean
    training_lower[h, "base_w_ICC"] <- fore_h$lower[,"95%"]
    training_upper[h, "base_w_ICC"] <- fore_h$upper[,"95%"]
  
  # restimacion del modelo sin ICC
    base_mod_noicc_h <- 
      Arima(ts(df_lags[1:(Ty - (Ht - h + 1)),Y], 
            start = ts_icc_lags_start, frequency = 12), 
            order = c(p, d, q), seasonal= c(P, D, Q), 
            xreg = cbind(curr_mat_outliers[1:(Ty - (Ht - h + 1)),]))
  
    coefs_base_bis[[h]] <- coeftest(base_mod_noicc_h)
    
    
    
    # forecast
    n <- ncol(curr_mat_outliers)
    new_xreg <- matrix(NA, 1, n)
    new_xreg[,1:n] <-
      unlist(curr_mat_outliers[(Ty - (Ht - h)):(Ty - (Ht - h)),])
    colnames(new_xreg) <- colnames(curr_mat_outliers)
  
    fore_h <- forecast(base_mod_noicc_h, xreg = new_xreg)
    mat_rmse[k,"base_w_o_ICC"] <- 
    sqrt((fore_h$mean - 
            econ_dl[(Ty - (Ht - h )):(Ty - (Ht - h)),Y])^2)
  
    training_forecasts[h, "base_w_o_ICC"] <- fore_h$mean
    training_lower[h, "base_w_o_ICC"] <- fore_h$lower[,"95%"]
    training_upper[h, "base_w_o_ICC"] <- fore_h$upper[,"95%"]
    
  
    # 2
    
    # RegARIMA incluyendo ICC
    lags_list <- list(c("ICC"),
                      c("ICC", "ICC_lag1"),
                      c("ICC", "ICC_lag1", "ICC_lag2"),
                      c("ICC", "ICC_lag1", "ICC_lag2", "ICC_lag3"),
                      c("ICC", "ICC_lag1", "ICC_lag2", "ICC_lag3", "ICC_lag12"))
    
    regarima_icc_list <- list()
    for(i in 1:length(lags_list)){
      xreg <- cbind(df_lags[1:(Ty - (Ht - h + 1)),as.character(lags_list[[i]])], 
                    curr_mat_outliers[1:(Ty - (Ht - h + 1)),], 
                    var_aux_rel[1:(Ty - (Ht - h + 1)),sel_var_aux])
      regarima_icc_list[[i]] <- 
        auto.arima(ts(econ_dl[1:(Ty - (Ht - h + 1)),Y], 
                   frequency = 12, start = ts_icc_lags_start),
                   allowdrift = FALSE, d = 0,
                   xreg = as.matrix(xreg))
    }
    
    model_sel <- 
      regarima_icc_list[[which.min(sapply(regarima_icc_list, BIC))]]
    lag_sel <- lags_list[[which.min(sapply(regarima_icc_list, BIC))]]
    print(paste("regarima_icc best lag: ",lag_sel, 
                "-", model_sel, sep = ""))
    
    p <- model_sel$arma[1]
    d <- model_sel$arma[6]
    q <- model_sel$arma[2]
    P <- model_sel$arma[3]
    D <- model_sel$arma[7]
    Q <- model_sel$arma[4]

    coefs_reg_icc[[h]] <- coeftest(model_sel)
      
    
    # forecast
    n <- length(lag_sel) + ncol(curr_mat_outliers) + length(sel_var_aux)
    ind1 <- c(1, 
             1 + length(lag_sel) ,
             1 + length(lag_sel) + ncol(curr_mat_outliers))
               
    ind2 <- c(length(lag_sel),
              length(lag_sel) + ncol(curr_mat_outliers),
              length(lag_sel) + ncol(curr_mat_outliers) + 
                length(sel_var_aux))
    
    new_xreg <- matrix(NA, 1, n)
    new_xreg[,ind1[1]:ind2[1]] <- 
      unlist(df_lags[(Ty - (Ht - h)):(Ty - (Ht - h)),as.character(lag_sel)])
    new_xreg[,ind1[2]:ind2[2]] <- 
      curr_mat_outliers[(Ty - (Ht - h)):(Ty - (Ht - h)),]
    new_xreg[, ind1[3]:ind2[3]] <- 
      unlist(var_aux_rel[(Ty - (Ht - h)):(Ty - (Ht - h)), sel_var_aux])
    colnames(new_xreg) <- c(as.character(lag_sel), colnames(curr_mat_outliers), sel_var_aux)
    
    fore_h <- forecast(model_sel, xreg = new_xreg)
    mat_rmse[k,"regarima_w_ICC"] <- 
      sqrt((fore_h$mean - 
              econ_dl[(Ty - (Ht - h )):(Ty - (Ht - h)),Y])^2)
    
    training_forecasts[h, "regarima_w_ICC"] <- fore_h$mean
    training_lower[h, "regarima_w_ICC"] <- fore_h$lower[,"95%"]
    training_upper[h, "regarima_w_ICC"] <- fore_h$upper[,"95%"]
    
    
    # reestimamos RegARIMA quitando ICC
    xreg <- cbind(curr_mat_outliers[1:(Ty - (Ht - h + 1)),], 
                  var_aux_rel[1:(Ty - (Ht - h + 1)),sel_var_aux])
    regarima_noicc_h <- 
      Arima(ts(econ_dl[1:(Ty - (Ht - h + 1)),Y], 
                    start = ts_icc_lags_start, frequency = 12), 
                              order = c(p, d, q), seasonal= c(P, D, Q), 
                              xreg = as.matrix(xreg)  )
    
    coefs_reg_icc_bis[[h]] <- coeftest(regarima_noicc_h)
      
    # forecast
    n <- ncol(curr_mat_outliers) + length(sel_var_aux)
    ind1 <- c(1, 1 + ncol(curr_mat_outliers))
    
    ind2 <- c(ncol(curr_mat_outliers),
              ncol(curr_mat_outliers) + length(sel_var_aux))
    
    new_xreg <- matrix(NA, 1, n)
    new_xreg[,ind1[1]:ind2[1]] <- 
      unlist(curr_mat_outliers[(Ty - (Ht - h)):(Ty - (Ht - h)),])
    new_xreg[,ind1[2]:ind2[2]] <- 
      unlist(var_aux_rel[(Ty - (Ht - h)):(Ty - (Ht - h)), sel_var_aux])
    colnames(new_xreg) <- c(colnames(curr_mat_outliers), sel_var_aux)
    
    fore_h <- forecast(regarima_noicc_h, xreg = new_xreg)
    mat_rmse[k,"regarima_w_o_ICC"] <- 
      sqrt((fore_h$mean - 
              econ_dl[(Ty - (Ht - h )):(Ty - (Ht - h)),Y])^2)
    
    training_forecasts[h, "regarima_w_o_ICC"] <- fore_h$mean
    training_lower[h, "regarima_w_o_ICC"] <- fore_h$lower[,"95%"]
    training_upper[h, "regarima_w_o_ICC"] <- fore_h$upper[,"95%"]
    
    # 3
    
    # RegARIMA incluyendo google trends
    lags_list <- list(c("bcd", "nd", "serv"),
                      c("bcd", "nd", "serv","bcd_lag1", "nd_lag1", "serv_lag1"),
                      c("bcd", "nd", "serv",
                        "bcd_lag1", "nd_lag1", "serv_lag1",
                        "bcd_lag2", "nd_lag2", "serv_lag2"),
                      c("bcd", "nd", "serv",
                        "bcd_lag1", "nd_lag1", "serv_lag1",
                        "bcd_lag2", "nd_lag2", "serv_lag2",
                        "bcd_lag3", "nd_lag3", "serv_lag3"),
                      c("bcd", "nd", "serv",
                        "bcd_lag1", "nd_lag1", "serv_lag1",
                        "bcd_lag2", "nd_lag2", "serv_lag2",
                        "bcd_lag3", "nd_lag3", "serv_lag3",
                        "bcd_lag12", "nd_lag12", "serv_lag12"))
    
    regarima_gt_list <- list()
    for(i in 1:length(lags_list)){
      xreg <- cbind(df_lags[1:(Ty - (Ht - h + 1)),as.character(lags_list[[i]])], 
                    curr_mat_outliers[1:(Ty - (Ht - h + 1)),], 
                    var_aux_rel[1:(Ty - (Ht - h + 1)),sel_var_aux])
      regarima_gt_list[[i]] <- 
        auto.arima(ts(econ_dl[1:(Ty - (Ht - h + 1)),Y], 
                      frequency = 12, start = ts_icc_lags_start),
                   allowdrift = FALSE, d = 0,
                   xreg = as.matrix(xreg))
    }
    
    model_sel <- 
      regarima_gt_list[[which.min(sapply(regarima_gt_list, BIC))]]
    lag_sel <- lags_list[[which.min(sapply(regarima_gt_list, BIC))]]
    print(paste("regarima_gt best lag: ",paste(lag_sel, collapse = "-"), 
                "-", model_sel, sep = ""))
    
    p <- model_sel$arma[1]
    d <- model_sel$arma[6]
    q <- model_sel$arma[2]
    P <- model_sel$arma[3]
    D <- model_sel$arma[7]
    Q <- model_sel$arma[4]
    
    coefs_reg_gt[[h]] <- coeftest(model_sel)
    
    
    # forecast
    n <- length(lag_sel)  + ncol(curr_mat_outliers) + length(sel_var_aux)
    ind1 <- c(1, 
              1 + length(lag_sel) ,
              1 + length(lag_sel) + ncol(curr_mat_outliers))
    
    ind2 <- c(length(lag_sel),
              length(lag_sel) + ncol(curr_mat_outliers),
              length(lag_sel) + ncol(curr_mat_outliers) + 
                length(sel_var_aux))
    
    new_xreg <- matrix(NA, 1, n)
    new_xreg[,ind1[1]:ind2[1]] <- 
      unlist(df_lags[(Ty - (Ht - h)):(Ty - (Ht - h)),as.character(lag_sel)])
    new_xreg[,ind1[2]:ind2[2]] <- 
      curr_mat_outliers[(Ty - (Ht - h)):(Ty - (Ht - h)),]
    new_xreg[,ind1[3]:ind2[3]] <- 
      unlist(var_aux_rel[(Ty - (Ht - h)):(Ty - (Ht - h)), sel_var_aux])
    colnames(new_xreg) <- c(colnames(gt_seas)[-1], 
                            colnames(curr_mat_outliers), sel_var_aux)
    
    fore_h <- forecast(model_sel, xreg = new_xreg)
    mat_rmse[k,"regarima_w_gt"] <- 
      sqrt((fore_h$mean - 
              econ_dl[(Ty - (Ht - h )):(Ty - (Ht - h)),Y])^2)
    
    training_forecasts[h, "regarima_w_gt"] <- fore_h$mean
    training_lower[h, "regarima_w_gt"] <- fore_h$lower[,"95%"]
    training_upper[h, "regarima_w_gt"] <- fore_h$upper[,"95%"]
    
    
    # reestimamos RegARIMA quitando google trends
    xreg <- cbind(curr_mat_outliers[1:(Ty - (Ht - h + 1)),], 
                  var_aux_rel[1:(Ty - (Ht - h + 1)),sel_var_aux])
    regarima_nogt_h <- 
      Arima(ts(df_lags[1:(Ty - (Ht - h + 1)),Y], 
               start = ts_icc_lags_start, frequency = 12), 
            order = c(p, d, q), seasonal= c(P, d, Q), 
            xreg = as.matrix(xreg)  )
    
    coefs_reg_gt_bis[[h]] <- coeftest(regarima_nogt_h)
      
    # forecast
    n <- ncol(curr_mat_outliers) + length(sel_var_aux)
    ind1 <- c(1, 
              1 + ncol(curr_mat_outliers))
    
    ind2 <- c(ncol(curr_mat_outliers),
              ncol(curr_mat_outliers) + length(sel_var_aux))
    
    new_xreg <- matrix(NA, 1, n)
    new_xreg[,ind1[1]:ind2[1]] <- 
      unlist(curr_mat_outliers[(Ty - (Ht - h)):(Ty - (Ht - h)),])
    new_xreg[,ind1[2]:ind2[2]] <- 
      unlist(var_aux_rel[(Ty - (Ht - h)):(Ty - (Ht - h)), sel_var_aux])
    colnames(new_xreg) <- c(colnames(curr_mat_outliers), sel_var_aux)
    
    fore_h <- forecast(regarima_nogt_h, xreg = new_xreg)
    mat_rmse[k,"regarima_w_o_gt"] <- 
      sqrt((fore_h$mean - 
              econ_dl[(Ty - (Ht - h )):(Ty - (Ht - h)),Y])^2)
    
    training_forecasts[h, "regarima_w_o_gt"] <- fore_h$mean
    training_lower[h, "regarima_w_o_gt"] <- fore_h$lower[,"95%"]
    training_upper[h, "regarima_w_o_gt"] <- fore_h$upper[,"95%"]
    
  k <- k + 1
  
}
Sys.time()

rmse_models <- apply(mat_rmse, 2, 
       function (x) sum(x, na.rm = TRUE) / sum(!is.na(x)== TRUE))

# figure8
par(mfrow = c(1,1), mai = c(1.5,.8,.5,.5))
barplot(rmse_models, beside = TRUE, 
        col = gray.colors(2), las = 2, cex.names = 0.7, 
        main = "Rolling RMSE fuera de muestra (h = {1,...,24})")

# figure 9 
par(mfrow = c(1,1), mai = c(.5,.5,.8,1.5), xpd = TRUE)
ts.plot(mat_rmse, xlab = "Forecast horizon", 
        ylab = "RMSE", col = c(4,2,3,1),
        lty = c(1,2,1,2), lwd = c(1,2,1,2))
legend(x = 25, y = max(mat_rmse, na.rm = TRUE) , 
       legend = colnames(mat_rmse),
         bty = "n",
       col = c(4,2,3,1), cex = 0.7,
       lty = c(1,2,1,2), lwd = c(1,2,1,2))


# figure 10
par(mfrow = c(3,2), mai = c(.5, .5, .5,.5))
# m1
ts.plot(ts(econ_dl[,Y], frequency = 12,start = ts_icc_lags_start), 
        ylab = "", xlab = "", col = c(4,2),
        main = paste(Y," con ICC", sep = ""))
lines(ts(c(rep(NA, nrow(econ_dl)-Ht), training_forecasts[,"base_w_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_lower[,"base_w_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_upper[,"base_w_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
#m2
ts.plot(ts(econ_dl[,Y], frequency = 12,start = ts_icc_lags_start), 
        ylab = "", xlab = "", col = c(4,2),
        main = paste(Y, " sin ICC", sep = ""))
lines(ts(c(rep(NA, nrow(econ_dl)-Ht), training_forecasts[, "base_w_o_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_lower[,"base_w_o_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_upper[,"base_w_o_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
# m3
ts.plot(ts(econ_dl[,Y], frequency = 12,start = ts_icc_lags_start), 
        ylab = "", xlab = "", col = c(4,2),
        main = paste(Y," +RegARIMA con ICC", sep = ""))
lines(ts(c(rep(NA, nrow(econ_dl)-Ht), training_forecasts[,"regarima_w_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_lower[,"regarima_w_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_upper[,"regarima_w_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
#m4
ts.plot(ts(econ_dl[,Y], frequency = 12,start = ts_icc_lags_start), 
        ylab = "", xlab = "", col = c(4,2),
        main = paste(Y, " +RegARIMA sin ICC", sep = ""))
lines(ts(c(rep(NA, nrow(econ_dl)-Ht), training_forecasts[,"regarima_w_o_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_lower[,"regarima_w_o_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_upper[,"regarima_w_o_ICC"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
# m5
ts.plot(ts(econ_dl[,Y], frequency = 12,start = ts_icc_lags_start), 
        ylab = "", xlab = "", col = c(4,2),
        main = paste(Y, "+RegARIMA con Google Trends", sep = ""))
lines(ts(c(rep(NA, nrow(econ_dl)-Ht), training_forecasts[,"regarima_w_gt"]),
         frequency = 12, start = ts_icc_lags_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_lower[,"regarima_w_gt"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_upper[,"regarima_w_gt"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
#m6
ts.plot(ts(econ_dl[,Y], frequency = 12,start = ts_icc_lags_start), 
        ylab = "", xlab = "", col = c(4,2),
        main = paste(Y, "+RegARIMA sin Google Trends"), sep = "")
lines(ts(c(rep(NA, nrow(econ_dl)-Ht), training_forecasts[, "regarima_w_o_gt"]),
         frequency = 12, start = ts_icc_lags_start), col = 2, lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_lower[,"regarima_w_o_gt"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)
lines(ts(c(rep(NA, nrow(econ_dl)-Ht),  training_upper[,"regarima_w_o_gt"]),
         frequency = 12, start = ts_icc_lags_start), col = "lightblue", lwd = 2)


# figure 11
par(mfrow = c(1,1), mai = c(1.5,.8,.5,.5))
barplot(c(obs_apr_20 = econ_dl[nrow(econ_dl)-1,Y],
          training_forecasts[nrow(training_forecasts),]), 
        beside = TRUE, col = gray.colors(2), las = 2, 
        cex.names = 0.7, main = "Pronosticos del consumo")


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

#save.image(paste(path, "icc-savf-iccDL.RData"))


