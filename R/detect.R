#' Classfy patients
#' @description This function reads a patient's SNV patterns and return the classfication result
#' @param x Patient's SNV patterns, should be a list of six matrices
#' @param cancer The list of cancers that processed in getMatrix.R
#' @param k The k of KNN
#' @param ignoreSame Dicedes whether ignore the exactly same matrices in all six dimensions
#' @return The result and param parameter same will be TRUE if there was a same patient
detect = function(x , cancer , k = 6 , ignoreSame = TRUE){

  maxId = 0
  for(i in 1:length(cancer)){
    if(length(cancer[[i]][[7]]) > maxId){
      maxId = length(cancer[[i]][[7]])
    }
  }
  distance = matrix(nrow = maxId , ncol = length(cancer))
  same = FALSE
  for(i in 1:length(cancer)){
    for(j in 1:length(cancer[[i]][[7]])){
      temp = dist(rbind(x[[1]],cancer[[i]][[1]][j,]),"euclidean") + dist(rbind(x[[2]],cancer[[i]][[2]][j,]),"euclidean") + dist(rbind(x[[3]],cancer[[i]][[3]][j,]),"euclidean") + dist(rbind(x[[4]],cancer[[i]][[4]][j,]),"euclidean") + dist(rbind(x[[5]],cancer[[i]][[5]][j,]),"euclidean") + dist(rbind(x[[6]],cancer[[i]][[6]][j,]),"euclidean")
      if(temp == 0){
        same = TRUE
        if(ignoreSame){
          next
        }else{
          distance[j,i] = temp
        }
      }else{
        distance[j,i] = temp
      }
    }
  }
  minK_yIndex = arrayInd(order(distance,decreasing=FALSE)[1:k],dim(distance))[,1]
  minK_xIndex = arrayInd(order(distance,decreasing=FALSE)[1:k],dim(distance))[,2]
  cat(minK_yIndex,"\n")
  cat(minK_xIndex,"\n")
  classfication = c()
  for(i in 1:length(cancer)){
    classfication[i] = 0
  }
  for(i in minK_xIndex){
    classfication[i] = classfication[i] + 1
  }
  cat(classfication,"\n")

  k_dis = c()
  for(i in 1:length(cancer)){
    k_dis[i] = 0
  }
  for(i in 1:k){
    k_dis[minK_xIndex[i]] = k_dis[minK_xIndex[i]] + distance[ minK_yIndex[i],minK_xIndex[i]]
  }
  k_backup = which(classfication == max(classfication),arr.ind=TRUE)
  k_min = as.numeric(max(k_dis))
  for(i in k_backup){
    if(k_dis[i] <= k_min){
      k_min = k_dis[i]
    }
  }
  k_result = which(k_dis == k_min)
  return(list(k_result,same))
}
