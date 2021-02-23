#require(rjags)

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

zcalc = function(samples, credMass=0.95, CI=c(.025,.975), delta=0, filters=NULL) {
  # Number of models passed to the function
  sampleCount = 1
  
  # If only one model is passed make a list of it
  if(!is.mcmc.list(samples)) {
    sampleCount = length(samples)
  }else{
    samples = list(samples)
  }
  
  # Make sure filter is a vector
  if(!is.null(filters) && !is.vector(filters))
  {
    filters = c(filters)
  }
  
  # Pre Calculate the number of expected rows
  rowCount = 0
  drops = list()
  for(k in 1:sampleCount) {
    drops[[k]] = c(0)
    if(is.null(filters)) {
      rowCount = rowCount + nvar(samples[[k]])
    }else{
      nvar = nvar(samples[[k]])
      varnames = varnames(samples[[k]])
      for(l in 1:nvar) {
        varname = varnames[[l]]
        isOK = FALSE
        for(f in 1:length(filters)) {
          isOK = isOK || regexpr(filters[[f]], varname)[1] > 0
        }
        if(isOK) {
          rowCount = rowCount + 1
        }else{
          drops[[k]] = c(drops[[k]],l)
        }
      }
    }
  }
  
  columnNames = c()
  # Pre allocate the df
  result = data.frame(
    mean=rep(NaN,rowCount),
    median=rep(NaN,rowCount),
    mode=rep(NaN,rowCount),
    sd=rep(NaN,rowCount),
    hdiLow=rep(NaN,rowCount),
    hdiHigh=rep(NaN,rowCount),
    quantileLow=rep(NaN,rowCount),
    quantileHigh=rep(NaN,rowCount),
    SS=rep(NaN,rowCount),
    ESS=rep(NaN,rowCount),
    RHAT=rep(NaN,rowCount),
    stringsAsFactors=FALSE
  )
  
  # Keeping track of the currently edited row
  currentRow = 0
  
  # Process the models
  for(k in 1:sampleCount) {
    # Make the name prefix if multiple models are present
    prefix = ""
    if( sampleCount > 1 ) {
      prefix = paste(k,".",sep="")
    }
    
    # Get the sample
    sample = samples[[k]]
    
    # Some common values
    variables = nvar(sample)
    varnames = varnames(sample)
    iterations = niter(sample)
    chains = nchain(sample)
    
    for(j in 1:variables) {
      if(!(j %in% drops[[k]])) {
        currentRow = currentRow + 1
        
        uvalue = unlist(sample[,j])
        value = sample[,j]
        
        columnNames = c(columnNames, paste(prefix,varnames[[j]],sep=""))
        
        result[currentRow,"SS"] <- iterations * chains
        result[currentRow,"ESS"] <- as.integer(round(effectiveSize(uvalue),1))
        result[currentRow,"mean"] <- mean(uvalue)
        result[currentRow,"median"] <- median(uvalue)
        
        mcmcDensity = density(uvalue)
        result[currentRow,"mode"] <- mcmcDensity$x[which.max(mcmcDensity$y)]
        
        HDI = HDIofMCMC( uvalue , credMass )
        result[currentRow,"hdiLow"] <- HDI[1]
        result[currentRow,"hdiHigh"] <- HDI[2]
        
        resultCI = quantile(uvalue, CI)
        
        result[currentRow,"quantileLow"] <- resultCI[1]
        result[currentRow,"quantileHigh"] <- resultCI[2]
        
        result[currentRow,"sd"] <- sd(uvalue)
        
        # RHAT calc
        
        # Get chain stats
        chainmeans = c()
        chainvars = c()
        for(i in 1:chains) {
          sum = sum(value[[i]])
          var = var(value[[i]])
          mean = sum / iterations
          
          chainmeans = c(chainmeans,mean)
          chainvars = c(chainvars,var)
        }
        globalmean = sum(chainmeans) / chains;
        
        #w in gelmanrubin with code === value
        #w = sum(var(coda)) / nchains;
        globalvar = sum(chainvars) / chains;
        
        # Compute between- and within-variances and MPV
        b = sum((chainmeans - globalmean)^2) * iterations / (chains - 1);
        
        varplus = (iterations - 1) * globalvar / iterations + b / iterations;
        
        # Gelman-Rubin statistic
        rhat = sqrt(varplus / globalvar);
        
        result[currentRow,"RHAT"] <- rhat
      }
    }
  }
  
  # Round a bit
  result = data.frame(apply(result, 2, function(x) round(x,4)))
  
  # Rename columns
  if(length(result) > 0) {
    names(result)[names(result) == 'hdiLow'] <- paste(sprintf("%.0f", round(credMass*100, digits = 2)),"HDI_L",sep="% ")
    names(result)[names(result) == 'hdiHigh'] <- paste(sprintf("%.0f", round(credMass*100, digits = 2)),"HDI_H",sep="% ")
    names(result)[names(result) == 'quantileLow'] <- paste("CrI",sprintf("%.1f%%", round(CI[1]*100, digits = 3)),sep=" ")
    names(result)[names(result) == 'quantileHigh'] <- paste("CrI", sprintf("%.1f%%", round(CI[2]*100, digits = 3)),sep=" ")
  }
  
  # Set the row names
  row.names(result) <- columnNames
  
  # Return
  result
}
