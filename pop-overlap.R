########This new intersect function calculates the intersection of two vectors considering duplicates
new_intersect<- function(a,b){
  rep(sort(intersect(a, b)), pmin(table(a[a %in% b]), table(b[b %in% a]))) 
}

### This function is able to expand the population from the original format to a vector with 10000 elements
expansion<-function(filename,line,bias,supply){
  expanded<-c()
  for (i in 1:(length(scan(paste(filename,"_",supply,"_",bias,".dat",collapse=NULL,sep=""), nlines = 1, what = "character",comment.char = "#",skip = line,quiet=TRUE))-1)) {
    test<-strsplit(scan(paste(filename,"_",supply,"_",bias,".dat",collapse=NULL,sep=""), nlines = 1, what = "character",comment.char = "#",skip = line,quiet=TRUE)[i+1],split = "")
    expanded<-c(expanded,rep(as.numeric(paste(test[[1]][1:(which(test[[1]]==":")-1)],sep = "",collapse = "")),times=as.numeric(paste(test[[1]][(which(test[[1]]==":")+1):(length(test[[1]]))],sep = "",collapse = ""))))
  }
  return(expanded)
}

############################################################# MAIN
  
args = commandArgs(trailingOnly=TRUE)

  filename<-args[1]
  supply<-args[2]
  bias_jacc<-c()
  lines_jacc<-c()
  count<-0
  record_alpha<-c()
  record_jacc<-c()
  for (k in seq(1,100,10)) {
    newk<-k+9
      temp_overlap<-0
      count<-count+1
      for (i in seq(k,newk-1,1)) {
        for (q in seq(i+1,newk,1)) {
          temp_overlap<-temp_overlap+length(new_intersect(expansion(filename,i,"0.50",supply),expansion(filename,q,"0.50",supply)))/10000
        }
          temp_bias<-0
          for (l in format(seq(0.05,0.90, by=0.05),2)) {
            newl<-format(as.numeric(l)+0.05,2)
            for (p in format(seq(newl,0.95, by=0.05),2)) {
              temp_bias<-temp_bias+length(new_intersect(expansion(filename,i,l,supply),expansion(filename,i,p,supply)))/10000
              record_alpha<-c(record_alpha,abs(as.numeric(p)-as.numeric(l)))
              record_jacc<-c(record_jacc,length(new_intersect(expansion(filename,i,l,supply),expansion(filename,i,p,supply)))/10000)
            }
          }
          temp_bias<-temp_bias/171
          bias_jacc<-c(bias_jacc,temp_bias)
      }
      temp_bias<-0
      for (l in format(seq(0.05,0.90, by=0.05),2)) {
        newl<-format(as.numeric(l)+0.05,2)
        for (p in format(seq(newl,0.95, by=0.05),2)) {
          temp_bias<-temp_bias+length(new_intersect(expansion(filename,i+1,l,supply),expansion(filename,i+1,p,supply)))/10000
          record_alpha<-c(record_alpha,abs(as.numeric(p)-as.numeric(l)))
          record_jacc<-c(record_jacc,length(new_intersect(expansion(filename,i+1,l,supply),expansion(filename,i+1,p,supply)))/10000)
        }
      }
      temp_bias<-temp_bias/171
      bias_jacc<-c(bias_jacc,temp_bias)
      temp_overlap<-temp_overlap/45
      lines_jacc<-c(lines_jacc,temp_overlap)
  }
  
  write.table(bias_jacc, paste(filename,"_",supply,"_biasjacc.dat",sep = "",collapse = ""), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
  
  write.table(lines_jacc, paste(filename,"_",supply,"_linesjacc.dat",sep = "",collapse = ""), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
  
  write.table(record_alpha, paste(filename,"_",supply,"_alpha.dat",sep = "",collapse = ""), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
  
  write.table(record_jacc, paste(filename,"_",supply,"_record.dat",sep = "",collapse = ""), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
