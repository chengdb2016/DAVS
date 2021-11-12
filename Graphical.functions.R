###### Functions needed #######################
## Linear regression ##########
lm.cov <- function (C, y, x) {
  solve(C[x, x], C[x, y, drop = FALSE])[1, ]
}

## Function that computes the set Possible-D-SEP(X,Y), modified by DMalinsky from version done by Spirtes
pdsepset.reach <- function(a,b,c,adjacency)
{
  ## reachable      set of vertices;
  ## edgeslist      array[1..maxvertex] of list of edges
  ## numvertex      integer
  ## labeled        array (by depth) of list of edges that have been labeled
  
  makeedge <- function(x,y) list(list(x,y))
  
  legal.dsep <- function(r,s) {
    ## Modifying global 'edgeslist'
    if (((adjacency[r[[1]],r[[2]]] == 2 &&
          adjacency[s,     r[[2]]] == 2 && r[[1]] != s) || ((adjacency[r[[1]],s] != 0 && r[[1]] != s))) &&  (is.poss.ancestor(s,a,adjacency) || is.poss.ancestor(s,b,adjacency))    && (is.poss.ancestor(r[[2]],a,adjacency) || is.poss.ancestor(r[[2]],b,adjacency))) {
      edgeslist[[r[[2]]]] <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)
    }
  }
  
  
  initialize.pdsep <- function(x,y) mapply(makeedge, x=x, y=y)
  
  labeled <- list()
  numvertex <- dim(adjacency)[1]
  edgeslist <- list()
  
  for (i in 1:numvertex)
    edgeslist <- c(edgeslist,list(which(adjacency[,i] != 0)))
  labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  edgeslist[[a]] <- list()
  depth <- 2
  repeat {
    labeled[[depth]] <- list()
    for (i in seq_along(labeled[[depth-1]])) {
      lab.i <- labeled[[depth-1]][[i]]
      edgestemp <- edgeslist[[lab.i[[2]]]]
      if (length(edgestemp) == 0) break
      for (j in seq_along(edgestemp))
        labeled[[depth]] <- union(legal.dsep(lab.i, edgestemp[[j]]),
                                  labeled[[depth]])
    }
    if (length(labeled[[depth]]) == 0)
      break
    ## else :
    depth <- depth  + 1
  }
  dsep <- unique(unlist(labeled))
  dsep
} # end function

reach <- function(a,b,c,adjacency)
{
  ## reachable      set of vertices;
  ## edgeslist      array[1..maxvertex] of list of edges
  ## numvertex      integer
  ## labeled        array (by depth) of list of edges that have been labeled
  
  makeedge <- function(x,y) list(list(x,y))
  
  legal.pdsep <- function(r,s) {
    ## Modifying global 'edgeslist'
    if ((adjacency[r[[1]],r[[2]]] == 2 &&
         adjacency[s,     r[[2]]] == 2 && r[[1]] != s) ||
        (adjacency[r[[1]],s] != 0 && r[[1]] != s)) {
      edgeslist[[r[[2]]]] <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)
    }
  }
  
  initialize.pdsep <- function(x,y) mapply(makeedge, x=x, y=y)
  
  labeled <- list()
  numvertex <- dim(adjacency)[1]
  edgeslist <- list()
  for (i in 1:numvertex)
    edgeslist <- c(edgeslist,list(which(adjacency[,i] != 0)))
  labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  edgeslist[[a]] <- list()
  depth <- 2
  repeat {
    labeled[[depth]] <- list()
    for (i in seq_along(labeled[[depth-1]])) {
      lab.i <- labeled[[depth-1]][[i]]
      edgestemp <- edgeslist[[lab.i[[2]]]]
      if (length(edgestemp) == 0) break
      for (j in seq_along(edgestemp))
        labeled[[depth]] <- union(legal.pdsep(lab.i, edgestemp[[j]]),
                                  labeled[[depth]])
    }
    if (length(labeled[[depth]]) == 0)
      break
    ## else :
    depth <- depth  + 1
  }
  unique(unlist(labeled))
}
### taken from pcalg package 11/27/2016 ###
possibleDe <- function(amat,x)
{
  ## Purpose: in a DAG, CPDAG, MAG, or PAG determine which nodes are
  ##          possible descendants of x on definite status paths
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: matrix corresponding to the DAG, CPDAG, MAG, or PAG
  ## - x: node of interest
  ## ----------------------------------------------------------------------
  ## Value:
  ## - de.list: array containing the possible descendants of x
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 26 Apr 2012, 16:58
  
  stopifnot(is.matrix(amat))
  p <- nrow(amat)
  is.de <- rep.int(FALSE, p) ##
  ## 1. case: x is a possible child of itself
  is.de[x] <- TRUE
  ## 2. case: find all the possible children of x
  indD <- which(amat[x,] != 0  & amat[,x] != 2 & !is.de) ## x (o,-)-* d
  i.pr <- rep(x,length(indD))
  while (length(indD) > 0) {
    ##next element in the queue
    d <- indD[1]
    indD <- indD[-1]
    pred <- i.pr[1]
    i.pr <- i.pr[-1]
    is.de[d] <- TRUE
    a.d <- amat[,d]
    a.d.p <- a.d[pred]
    ## find all possible children of d not visited yet
    indR <- which(amat[d,] != 0 & a.d != 2 & !is.de) ## d (o,-)-* r
    for(j in seq_along(indR)) {
      ## check that the triple <pred,d,r> is of a definite status
      ## 1. d is a collider on this subpath; this is impossible
      ##    because the edge between d and r cannot be into d
      ## 2. d is a definite non-collider
      r <- indR[j]
      if (a.d.p == 3 || a.d[r] == 3 ||
          (a.d.p == 1 && a.d[r] == 1 && amat[pred,r] == 0)) {
        ## update the queues
        indD <- c(indD, r)
        i.pr <- c(i.pr, d)
      }
    }
  }
  ## return 'de.list' :
  which(is.de)
  
} ## {possibleDe}

### is a an ancestor of b in graph g?
#fff.2 <- memoize(is.poss.ancestor)
is.poss.ancestor <- function(a, b, g,visited=NULL){
  if(a==b) return(TRUE)
  foundpath <- c()
  ind1 <- which(g==3, arr.ind=TRUE, useNames=FALSE) #tails
  ind11 <- which(g==1, arr.ind=TRUE, useNames=FALSE) #circles
  ind1 <- rbind(ind1,ind11) ## tails and circles
  ind1 <- subset(ind1,ind1[,2]==a) # array of indices for tails and circles out of A
  if(nrow(ind1)==0) return(FALSE) # addition: if there are no tails or circles at A
  for(x in 1:nrow(ind1)){ # loop through tails and circles out of A
    if(ind1[x,1] %in% visited) next
    if(g[ind1[x,2],ind1[x,1]]==2 || g[ind1[x,2],ind1[x,1]]==1){ # if there is an arrowhead or circle at the other end of the x-th tail (call this C)
      if(ind1[x,1]==b){
        foundpath <- append(foundpath,TRUE)
        break
      }
      if(any(g[,ind1[x,1]]==3 | g[,ind1[x,1]]==1)){ # if there are any tails or circles out of C
        a_old <- a
        a2 <- ind1[x,1]
        if(a2==a_old) next
        foundpath <- append(foundpath,is.poss.ancestor(a2,b,g,visited=c(visited,a_old)))
        if(any(foundpath)==TRUE) break
      }
    } # if there isn't an arrowhead at C - !(A-->C) - don't return anything
  } # for x in 1:nrow(ind1)
  if(any(foundpath)==TRUE) return(TRUE)
  else return(FALSE)
} # end function
