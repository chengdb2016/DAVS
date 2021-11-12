Cal_bias<-function(A,True_A){
  # Calculate bias, the proportion of 
  
  bias<-abs((A-True_A)/True_A)*100
  
  
  bias
}
Cal_ATEerror<-function(A,True_A){
  # Calculate bias, the proportion of 
  
  ATEerror<-abs(A-True_A)
  
  
  ATEerror
}
