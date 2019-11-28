#' Process MAF files
#' @description This function reads MAF files and extracts six mutli-dimensional SNV matrices
#' @param x MAF file
#' @return Six matrices and all the id of MAF file
getMatrix = function(x){

  # read maf
 #x = "lgg.maf"
 maf = read.csv(file = x, header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = "#")


  #Test for adjacent base data

  if("Genome_Plus_Minus_10_Bp" %in% colnames(maf)){
    gene_13 =maf[,"Genome_Plus_Minus_10_Bp"]
  }else if("ref_context" %in% colnames(maf)){
    gene_21 =maf[,"ref_context"]
    gene_13 = c()
    for(i in 1:length(gene_21)){
      gene_13[i] = toupper(substring(gene_21[i] ,5 ,17))
    }
  }else{
    return("no data")
  }

  # Read data

  gene_re =maf[,"Reference_Allele"]
  type =maf[,"Variant_Type"]
  sta =maf[,"Mutation_Status"]
  alle =maf[,"Tumor_Seq_Allele2"]
  gene_id =maf[,"Matched_Norm_Sample_Barcode"]
  id = unique(gene_id)

  # Replace the same base substitution type

  sig1_orig = c("G>A","G>T","G>C","A>G","A>T","A>C")
  sig1_repl = c("C>T","C>A","C>G","T>C","T>A","T>G")

  #Count subscripts that satisfy the conditions

  count = 0
  index = array()
  for( i in 1:length(gene_13)){
    if(type[i] == "SNP" & sta[i] == "Somatic" & gene_id[i] %in% id){
      count <- count + 1
      index[count] <- i
    }
  }
  errorNum = 0
  errorIndex = array()
  for(i in 1:count){
    if(substring(gene_13[index[i]],7,7) != gene_re[index[i]]){
      errorNum = errorNum+1
      errorIndex[errorNum] = index[i]
    }
  }
  index = setdiff(index , errorIndex)
  count = count - errorNum

  # Building a one-dimensional feature matrix

  mut1 = array()
  for(i in 1:count){
    mut1[i] = paste0(substr(gene_13[index[i]],7,7),">",alle[index[i]])
  }
  for(j in 1:count){
    for(i in 1:6){
      if(mut1[j] == sig1_orig[i]){
        mut1[j] = sig1_repl[i]
      }
    }
  }
  nu1 = matrix(0,nrow = length(id) ,ncol = 6 ,dimnames = list(id ,sig1_repl))
  for(i in 1:count){
    nu1[gene_id[index[i]],mut1[i]] = nu1[gene_id[index[i]],mut1[i]] + 1 #下标映射，计数加一
  }

  # Building a two-dimensional left-adjacent feature matrix

  mut2_l = array()
  for(i in 1:count){
    mut2_l[i] = paste0(substr(gene_13[index[i]],6,6),".",substr(gene_13[index[i]],7,7),">",alle[index[i]])
  }
  for(j in 1:count){
    for(i in 1:6){
      if(substr(mut2_l[j],3,5) == sig1_orig[i]){
        mut2_l[j] = paste0(substr(mut2_l[j] ,1,2),sig1_repl[i])
      }
    }
  }
  sig2_repl_l = unique(mut2_l)
  nu2_l = matrix(0,nrow = length(id) ,ncol = 24 ,dimnames = list(id ,sig2_repl_l))
  for(i in 1:count){
    nu2_l[gene_id[index[i]],mut2_l[i]] = nu2_l[gene_id[index[i]],mut2_l[i]] + 1
  }

  # Building a two-dimensional right-adjacent feature matrix

  mut2_r = array()
  for(i in 1:count){
    mut2_r[i] = paste0(substr(gene_13[index[i]],7,7),">",alle[index[i]],".",substr(gene_13[index[i]],8,8))
  }
  for(j in 1:count){
    for(i in 1:6){
      if(substr(mut2_r[j],1,3) == sig1_orig[i]){
        mut2_r[j] = paste0(sig1_repl[i] , substr(mut2_r[j] ,4,5) )
      }
    }
  }
  sig2_repl_r = unique(mut2_r)
  nu2_r = matrix(0,nrow = length(id) ,ncol = 24 ,dimnames = list(id ,sig2_repl_r))
  for(i in 1:count){
    nu2_r[gene_id[index[i]],mut2_r[i]] = nu2_r[gene_id[index[i]],mut2_r[i]] + 1
  }

  # Building a three-dimensional left-adjacent feature matrix

  mut3_l = array()
  for(i in 1:count){
    mut3_l[i] = paste0(substr(gene_13[index[i]],5,5),".",substr(gene_13[index[i]],6,6),".",substr(gene_13[index[i]],7,7),">",alle[index[i]])
  }
  for(j in 1:count){
    for(i in 1:6){
      if(substr(mut3_l[j],5,7) == sig1_orig[i]){
        mut3_l[j] = paste0(substr(mut3_l[j],1,4) , sig1_repl[i] )
      }
    }
  }
  sig3_repl_l = unique(mut3_l)
  nu3_l = matrix(0,nrow = length(id) ,ncol = 96 ,dimnames = list(id ,sig3_repl_l))
  for(i in 1:count){
    nu3_l[gene_id[index[i]] , mut3_l[i]] = nu3_l[gene_id[index[i]] , mut3_l[i]] + 1
  }

  # Building a three-dimensional both sides adjacent feature matrix

  mut3_m = array()
  for(i in 1:count){
    mut3_m[i] = paste0(substr(gene_13[index[i]],6,6),".",substr(gene_13[index[i]],7,7),">",alle[index[i]],".",substr(gene_13[index[i]],8,8))
  }
  for(j in 1:count){
    for(i in 1:6){
      if(substr(mut3_m[j],3,5) == sig1_orig[i]){
        mut3_m[j] = paste0(substr(mut3_m[j],1,2) ,sig1_repl[i], substr(mut3_m[j],6,7))
      }
    }
  }
  sig3_repl_m = unique(mut3_m)
  nu3_m = matrix(0,nrow = length(id) ,ncol = 96 ,dimnames = list(id ,sig3_repl_m))
  for(i in 1:count){
    nu3_m[gene_id[index[i]] , mut3_m[i]] = nu3_m[gene_id[index[i]] , mut3_m[i]] + 1
  }

  # Building a three-dimensional right-adjacent feature matrix

  mut3_r = array()
  for(i in 1:count){
    mut3_r[i] = paste0(substr(gene_13[index[i]],7,7),">",alle[index[i]],".",substr(gene_13[index[i]],8,8),".",substr(gene_13[index[i]],9,9))
  }
  for(j in 1:count){
    for(i in 1:6){
      if(substr(mut3_r[j],1,3) == sig1_orig[i]){
        mut3_r[j] = paste0(sig1_repl[i], substr(mut3_r[j],4,7) )
      }
    }
  }
  sig3_repl_r = unique(mut3_r)
  nu3_r = matrix(0,nrow = length(id) ,ncol = 96 ,dimnames = list(id ,sig3_repl_r))
  for(i in 1:count){
    nu3_r[gene_id[index[i]] , mut3_r[i]] = nu3_r[gene_id[index[i]] , mut3_r[i]] + 1
  }
  return(list(nu1,nu2_l,nu2_r,nu3_l,nu3_m,nu3_r,id))
}
