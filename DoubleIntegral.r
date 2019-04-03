rm(list = setdiff(ls(), lsf.str()))

################################## AUXILIARY ##################################

#---------
# STEP   \____________________________________________________________________

getHT <- function(a, b, N)
{
  h = (b-a)/(2*N)
  return(h)
}

getHQ <- function(a, b, N)
{
  h = (b-a)/2
  return(h)
}

getHS13 <- function(a, b, N)
{
  h = (b-a)/(3*N)
  return(h)
}

getHS38 <- function(a, b, N)
{
  h = (3*(b-a))/(8*N)
  return(h)
}

#---------
# ROOT  \____________________________________________________________________

getXYi <- function(a, b, N)
{
  h = (b-a)/N
  Xi = matrix(nrow = N+1)
  
  for (n in 1:(N+1)){
    Xi[n] = a+((n-1)*h)
  }
  
  return(Xi)
}

getRXYi <- function(a, b, N)
{
  
  h = (b-a)/2
  h2 = (b+a)/2
  
  if(N==2) u = matrix(c(-0.57735, 0.57735), nrow=2)
  else if(N==4) u = matrix(c(-0.8611361, -0.33998104, 0.33998104, 0.8611361), nrow=4)
  
  rows = dim(u)[1]
  Ri = matrix(nrow = rows)
  
  for(i in 1:rows){
    Ri[i] = (h*u[i]) + h2
  }
  
  return(Ri)
}

#---------
# WEIGHT   \____________________________________________________________________

getWT <- function(N)
{
  w = matrix(2, nrow = N+1)
  w[1] = 1
  w[N+1] = 1
  
  return(w)
}

getWS13 <- function(N)
{
  w = matrix(2, nrow = N+1)
  
  for(row in 2:N){
    if(row%%2==0){
      w[row]=4
    }
    else
      w[row]=2
  }
  
  w[1] = 1
  w[N+1] = 1
  
  return(w)
}

getWS38 <- function(N)
{
  w = matrix(2, nrow = N+1)
  
  for(row in 2:N){
    
    if(row%%4==0){
      w[row]=2
    }
    else
      w[row]=3
  }
  
  w[1] = 1
  w[N+1] = 1
  
  return(w)
}

getWQ <- function(N)
{
  if(N==2) return(matrix(1, nrow=2))
  else if(N==4) return(matrix(c(0.34785485, 0.65214515, 0.65214515, 0.34785485), nrow=4))
}

#-----------
# TABLES   \____________________________________________________________________

tableFxy <- function(xi, yi, exp)
{
  
  # Preparing table of f(x,y)
  
  fxy_table = matrix(nrow = dim(xi)[1], ncol = dim(yi)[1])
  names = matrix(nrow = dim(xi)[1])
  
  for(i in 1:dim(fxy_table)[1])
  {
    names[i] = paste(c("X[", i, "] = ", xi[i]), collapse = "")
  }
  
  rownames(fxy_table) = names
  
  names = matrix(nrow = dim(yi)[1])
  
  for(i in 1:dim(fxy_table)[2])
  {
    names[i] = paste(c("Y[", i, "] = ", yi[i]), collapse = "")
  }
  colnames(fxy_table) = names
  
  # Populating values
  
  for(i in 1:dim(fxy_table)[1]){
    for(j in 1:dim(fxy_table)[2]){
      
      x = xi[i]
      y = yi[j]
      fxy_table[i,j] = eval(exp)
      
    }
  }  
  
  
  return(fxy_table)
  
}

tableWeights <- function(wxi, wyi)
{
  # Preparing weight chart
  
  weights_table = matrix(nrow = dim(wxi)[1], ncol = dim(wyi)[1])
  names = matrix(nrow = dim(weights_table)[1])
  
  for(i in 1:dim(weights_table)[1])
  {
    names[i] = paste(c("Wx[", i, "]"), collapse = "")
  }
  rownames(weights_table) = names
  
  for(i in 1:dim(weights_table)[2])
  {
    names[i] = paste(c("Wy[", i, "]"), collapse = "")
  }
  colnames(weights_table) = names
  
  # Preparing values
  
  for(i in 1:dim(weights_table)[1]){
    for(j in 1:dim(weights_table)[2]){
      
      weights_table[i,j] = (wxi[i])*(wyi[j])
      
    }
  }
  
  return(weights_table)
}

tableWFxy <- function(fxy, weights)
{
  # Preparing final table of w*f(x,y)
  
  WFxy = matrix(nrow = dim(fxy)[1], ncol = dim(fxy)[2])
  names = matrix(nrow = dim(WFxy)[1])
  
  for(i in 1:dim(WFxy)[1])
  {
    names[i] = paste(c("Wxy[", i, "]"), collapse = "")
  }
  rownames(WFxy) = names
  
  for(i in 1:dim(WFxy)[2])
  {
    names[i] = paste(c("F(x,y)[", i, "]"), collapse = "")
  }
  colnames(WFxy) = names
  
  # Preparing values
  
  for(i in 1:dim(WFxy)[1]){
    for(j in 1:dim(WFxy)[2]){
      
      WFxy[i,j] = fxy[i,j]*weights[i,j]
      
    }
  }
  
  return(WFxy)
}

#-------------
# INTEGRATION  \____________________________________________________________________

sumWFxy <- function(WFxy)
{
  sum = 0
  for(i in 1:dim(WFxy)[1]){
    for(j in 1:dim(WFxy)[2]){
      sum = sum + WFxy[i,j]
    }
  }
  
  return(sum)
}

################################## MAIN FUNCTION ##############################

main <- function()
{
  
  #-----------------------------------------------------------------------#|
  #                           DATAS OF INTEGRATION                        #|
  #_______________________________________________________________________#|
  
  ax = 1
  bx = 2
  Nx = 4
  
  ay = 2
  by = 3
  Ny = 4
  
  exp = parse(text="(x^2)*y")
  methodx = "Q"
  methody = "S13"
  
  #-----------------------------------------------------------------------#|
  #                       EXECUTION OF CALCULATION                        #|
  #_______________________________________________________________________#|
                                                                          
  # Gets roots and weights in x                                             
  
  xi =   switch(methodx, 
          "T" = getXYi(ax, bx, Nx), 
          "Q" = getRXYi(ax, bx, Nx), 
          "S13" = getXYi(ax, bx, Nx), 
          "S38" = getXYi(ax, bx, Nx)
         )
                                                                            
  wxi =   switch(methodx, 
            "T" = getWT(Nx), 
            "Q" = getWQ(Nx), 
            "S13" = getWS13(Nx), 
            "S38" = getWS38(Nx)
          ) 
                                                                            
  # Gets roots and weights in y

  yi =   switch(methody, 
            "T" = getXYi(ay, by, Ny), 
            "Q" = getRXYi(ay, by, Ny), 
            "S13" = getXYi(ay, by, Ny), 
            "S38" = getXYi(ay, by, Ny)
         )
  
  wyi =   switch(methody, 
            "T" = getWT(Ny), 
            "Q" = getWQ(Ny), 
            "S13" = getWS13(Ny), 
            "S38" = getWS38(Ny)
          )
  

  # Auxiliary Tables
  
  fxy = tableFxy(xi, yi, exp)
  weights = tableWeights(wxi, wyi)
  
  # Final tables
  
  WFxy = tableWFxy(fxy, weights)
  
  # Sum and result
  
  sum = sumWFxy(WFxy)
  
  hx =  switch(methodx, 
                "T" = getHT(ax, bx, Nx), 
                "Q" = getHQ(ax, bx, Nx), 
                "S13" = getHS13(ax, bx, Nx), 
                "S38" = getHS38(ax, bx, Nx)
        )
  
  hy =  switch(methody, 
                "T" = getHT(ay, by, Ny), 
                "Q" = getHQ(ay, by, Ny), 
                "S13" = getHS13(ay, by, Ny), 
                "S38" = getHS38(ay, by, Ny)
        )
  
  h = hx*hy
  result = sum*h
  
  
  #-----------------------------------------------------------------------#|
  #                    PRINT OF RESULTS                                   #|
  #_______________________________________________________________________#|
  
  cat("\n       ")
  
  cat("#-----------------------------------------------------------------------#|
       #                               f(x,y)                                  #|
       #_______________________________________________________________________#|\n\n")
  
  print(fxy)
  
  cat("\n       ")
  
  cat("#-----------------------------------------------------------------------#|
       #                               Weights                                   #|
       #_______________________________________________________________________#|\n\n")
  
  print(weights)
  
  cat("\n       ")
  
  cat("#-----------------------------------------------------------------------#|
       #                           Weights x f(x,y)                              #|
       #_______________________________________________________________________#|\n\n")
  
  print(WFxy)
  
  cat("\n       ")
  
  cat("#-----------------------------------------------------------------------#|
       #                      Results of integration                           #|
       #_______________________________________________________________________#|\n\n")
  
  cat("\n\tSum: ")
  print(sum)
  cat("\n\thx * hy: ")
  print(h)
  cat("\n\tResult: ")
  sprintf(result, fmt = '%#.8f')
  
}

main()
