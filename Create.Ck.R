create_Ck<-function(Lksub1,k){
  countn<-1
  Ck<-list()
  len_Lksub1<-length(Lksub1)
  for (i in 1:(len_Lksub1-1)) {
    for (j in (i+1):len_Lksub1) {
      L1<-Lksub1[[i]]
      L2<-Lksub1[[j]]
      if(all((L1==L2)[-k])==TRUE & (L1==L2)[k]==FALSE){
        Ck_item<-union(L1,L2)
        temps <- combn(Ck_item, k,simplify = FALSE)
        require("fastmatch")
        if(all(temps %fin% Lksub1) == TRUE){
          Ck[[countn]]<- Ck_item
          countn<-countn+1
        }
      }
    }
  }
  return(Ck)
}




