#Like VARX lasso, but splits the data up to be run in parallel
map.VARXlasso = function(y, x, p, b, y.spec=matrix(1,nrow=nrow(y),ncol=nrow(y)), 
  x.spec=matrix(1,nrow=nrow(y),ncol=nrow(x)), computenodes=4, mapFolder) {

  varz = VARXZ(y,x,p,b)
  Z = varz$Z
  write.table(Z, file=paste(mapFolder,'/VARXlasso.Z.',p,'.',b,sep=''))
  write.table(y.spec, file=paste(mapFolder,'/VARXlasso.yspec.',p,'.',b,sep=''))
  write.table(x.spec, file=paste(mapFolder,'/VARXlasso.xspec.',p,'.',b,sep=''))
  
  partition.nrows = floor(ncol(y) / computenodes)
  for(i in 1:computenodes) {
    if(i == computenodes) {
      y.partition = varz$y.p
      write.table(y.partition, file=paste(mapFolder,'/VARXlasso.y',i,'.',p,'.',b,sep=''))
    } else {
      y.partition = varz$y.p[,1:partition.nrows,]
      write.table(y.partition, file=paste(mapFolder,'/VARXlasso.y',i,'.',p,'.',b,sep=''))
      varz$y.p = varz$y.p[,-(1:partition.nrows),]
    }
  }
}

.reduce.VARXlassocore = function(i, y, Z, p, b, y.spec, x.spec, lengthY, lengthX) {
  rowIndex = i[1]
  x.to.remove = which(!x.spec[rowIndex,])
  y.to.remove = which(!y.spec[rowIndex,])
  np = lengthY * p
  z.x.to.remove = c()
  z.y.to.remove = c()
  if(length(x.to.remove) > 0) {
    z.x.to.remove = as.vector(sapply(1:b, function(ii) {
      x.to.remove + np + lengthX*(ii-1)
    }))
  }
  if(length(y.to.remove) > 0) {
    z.y.to.remove = as.vector(sapply(1:p, function(ii) {
      y.to.remove + lengthY*(ii-1)
    }))
  }
  z.to.remove = unlist(c(z.x.to.remove, z.y.to.remove))
  if(length(z.to.remove > 0)) {
    Z.reduced = Z[,-z.to.remove]
  } else {
    Z.reduced = Z
  }
       
  i.lasso = cv.glmnet(Z.reduced, i[-1])
  i.lasso.s = i.lasso$lambda.min

  B = rep(0, ncol(Z)+1)  #make vector of coefficients full size
  B.reduced = rep(0, ncol(Z.reduced) + 1)

  i.coef = coef(i.lasso, i.lasso.s)
  B.reduced[(i.coef@i)+1] = i.coef@x
  if(length(z.to.remove) > 0) {
    B[-(z.to.remove+1)] = B.reduced
    return (B)
  }
  else {
    return (B.reduced)
  }
}

reduce.VARXlasso = function(y, x, Z, p, b, y.spec=matrix(1,nrow=nrow(y),ncol=nrow(y)), 
  x.spec=matrix(1,nrow=nrow(y),ncol=nrow(x)), lengthY, lengthX, id, 
  reduceFolder, numcore=1, ...) {

  y.augmented = rbind(1:ncol(y), y)
  if(numcore==1) {
    var.lasso = apply(y.augmented, 2, .reduce.VARXlassocore, y=y,Z=Z,p=p,b=b,
      y.spec=y.spec, x.spec=x.spec, lengthY, lengthX)
  } else {
    y.augmented.list = c(unname(as.data.frame(y.augmented)))  #converts matrix to list
    var.lasso = mclapply(y.augmented.list, .reduce.VARXlassocore, y=y,Z=Z,p=p,b=b,
      y.spec=y.spec, x.spec=x.spec, lengthY, lengthX, 
      mc.cores=numcore, ...)
  }

  #mclapply returns a list, so unlist it and recreate the B matrix
  var.lasso = matrix(unlist(var.lasso),nrow=(lengthY*p + lengthX*b+1), ncol=ncol(y))
  rownames(var.lasso) = c('intercept', rownames(Z))
  colnames(var.lasso) = colnames(y)
  write.table(var.lasso, file=paste(reduceFolder,'/VARXlasso.B.',id,'.',p,'.',b,sep=''))
}
 
gather.VARXlasso = function(foldername) {
  files = as.list(list.files(path=foldername))
  B = lapply(files, function(x) {
    read.table(paste(foldername,'/',x,sep=''))
  })
  do.call('cbind', B)
}

