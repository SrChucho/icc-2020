### funciones para el script de Consumo e ICC

# funcion de cambios porcentuales
cpm <- function(x, q = 1){
  
  index <- matrix(0, length(x) - q, 2)
  index[,1] <- 1:nrow(index)
  index[,2] <- (q+1):(nrow(index)+q)
  
  xm <- rep(0, nrow(index))
  for(i in 1 : length(xm))
    xm[i] <- x[index[i,2]]/x[index[i,1]]*100-100
  
  return(xm)
}

# lag
y_lag <- function(x, l = 1){ # x <- icc_inegi
  x_trans <- c(rep(NA,l), data.matrix(x))
  return(x_trans)  
}
