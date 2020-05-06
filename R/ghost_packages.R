require('stats')
require('R6')
require('pracma')
require('devtools')
{
  HTCls <- R6Class("HTCls",
                   lock_objects = FALSE,
                   portable = TRUE,
                   public = list(
                     contentASkey = data.frame(),
                     offsets = as.vector(0),
                     #class constructor
                     initialize = function (inContentASkey, inNewoffset) {
                       self$contentASkey = inContentASkey
                       self$offsets = c(self$offsets, inNewoffset)
                       #self$offsets = c(inNewoffset)
                     },
                     addOffset = function () {
                       append(self$offsets, newoffset)
                     },

                     findWindow = function (indf) {
                       return (exactidentical.norowname(indf,self$contentASkey))
                     }
                   )
  )
}
{
  MissCls <- R6Class("MissCls",
                     public = list(
                       post = data.frame(),
                       prior = data.frame(),
                       miss = data.frame(),
                       resolver = data.frame(),
                       #missrn =  as.integer(0),
                       missrn = as.vector(0),
                       confidence =  as.integer(0),
                       initialize = function (priorin,missin,postin,cnts) {
                         self$post = postin
                         self$prior = priorin
                         self$miss = missin
                         self$missrn = cnts
                       },
                       addConfidence = function () {
                         self$confidence = self$confidence+1
                       }
                     )
  )
}
reconstruct<-function(data_frame,constraintCol, wSize,direction_save,epsilon) {

  data_frame[data_frame == ""] <- NA

  ### converting the constraint name to number
  if(is.numeric(constraintCol)) {}
  else{
    constraintCol<-constraint_check(constraintCol,data_frame)
  }

  #### checking the insert dataframe
  if(is.numeric(data_frame[,-constraintCol])) { stop('The inserted dataframe is not discreted!!') }
  else{
    df1<-data_frame
    ht = slidewindow(data_frame, constraintCol, wSize)
    missingSN <- data.frame() # this dataframe holds missing data + their posterior and prior data
    cnt =1 # counter for while loop

    while (cnt < nrow(data_frame)) {
      cur = data_frame[cnt,]
      curNoCST = c(cur[,-c(constraintCol)])
      priordf <- data.frame()
      postdf <- data.frame()
      missingdf <- data.frame()
      if (all(is.na(curNoCST)) && (cnt > wSize) && (cnt <= (nrow(data_frame)-wSize)) ) { # an NA row has been found
        missingdf <- rbind(missingdf,data_frame[cnt,]) # add it to the missing column
        cntPri = 1
        cntPost = 1
        wPost = wSize
        wPri = wSize

        #### construct the prior part
        while (cntPri <= wPri) {
          xPri = data_frame[cnt-cntPri,]
          xnoCSTPri = c(xPri[,-c(constraintCol)]) # remove the Constraint column
          if (!all(is.na(xnoCSTPri))) {
          priordf <- rbind(xPri,priordf) # add prior elements that has data (only if it is not NA row)
          }
          else{
            missingdf <- rbind(missingdf,xPri) # prior one is missing and thus added to the missing list
            wPri <- wPri + 1 # this has been added to make another check and neglect the window check.
          }
          cntPri <- cntPri +1
        }
        nextMiss = FALSE
        #### construct the postrior part
        while (cntPost <= wPost) {

          xPos = data_frame[cnt + cntPost,]
          xnoCSTPos = c(xPos[,-c(constraintCol)]) # remove the constraint column
          if (!all(is.na(xnoCSTPos))) {
            postdf <- rbind(postdf, xPos)
            # add prior elements that has data (obviously it is not NA row)
          }else{
            wPost <- wPost +1 # increase the post window size
            missingdf <- rbind(missingdf,xPos) # prior one is missing
            nextMiss = TRUE
          }
          cntPost <- cntPost +1
        }

        if (!nextMiss) {
          priordf <- priordf[ order(row.names(priordf)), ]
          missingSN = rbind(missingSN,priordf)
          missingSN = rbind(missingSN,missingdf)
          missingSN = rbind(missingSN,postdf)
          missingSN <- missingSN[order(row.names(missingSN)),]
          obj <- MissCls$new(priordf, missingdf, postdf, rownames(missingdf))
          #--------------------------------------------------------------------------
          df1<-resolveNoSignal_revised(obj,df1,constraintCol,wSize,epsilon,cnt)
          #--------------------------------------------------------------------------
          missingSN = NULL
          priordf <- NULL
          postdf <- NULL
          missingdf <- NULL
        }
        nextMiss = FALSE
      }
      cnt <- cnt+1
    }


    if(check(direction_save)==0){
      output_file<-df1
      output_file
      return(output_file)
    }
    if(check(direction_save)==1){
      out_csv(df1,direction_save)
      message("The .csv file is saved in ",direction_save," direction")}
  }
}
saxTransform <- function(data_frame, buckets, skipColumnVec,constraint_row) {

  #Will take the data in

  if (buckets < 2 || buckets > 26) {
    message("Error: Buckets must be between 2 and 26!")
    return(NULL)
  }
  data <- as.matrix(data_frame)


  # if(nrow(data)>ncol(data)) {data<-t(data)}

  minData <- min(data,na.rm=TRUE)
  maxData <- max(data,na.rm=TRUE)
  width <- maxData - minData
  bucketWidth <- floor(width / buckets)


  ####-------------------convert num to char
  i=1
  j=1
  for (i in 1:ncol(data)){

    for(j in 1:nrow(data)){

      if(!is.na(data[j,i])){data[j,i]<-chr(97+(floor((as.numeric(data[j,i]) - as.numeric(minData)) / bucketWidth)))}
      j=j+1
    }
    i=i+1
  }
  ##-----------------removing the constraint rows and coloumns
  if(check(constraint_row)==1) { data[,constraint_row]<-data_frame[,constraint_row]  }
  if(check(skipColumnVec)==1) { data[skipColumnVec,]<-as.numeric(data_frame[skipColumnVec,] ) }
  data
}
###################---------------------------------------------------------------
resolveNoSignal_revised<-function(obj, df1, constraintCol,wSize, epsilon,cnt_first) {

  prior <- obj$prior
  prior<-as.data.frame(prior)
  posterior <- obj$post
  posterior<-as.data.frame(posterior)
  cnt_first1<-cnt_first
  miss <- obj$miss
  miss<-as.data.frame(miss)
  miss <- miss[order(row.names(miss)),]

  cnt <- 1
  curr <- data.frame()
  while (cnt < nrow(df1)){
    cnt2 <- cnt + (wSize-1)
    curr <- data.frame(df1[cnt:cnt2,])

    if (identical.norowname(prior,curr,epsilon)) {

      missEidx <- cnt2+(nrow(miss))
      missBidx <- cnt2+1
      potentialResolver <- data.frame(df1[missBidx:missEidx,])

      if (exactidentical.norowname(potentialResolver[,constraintCol],miss[,constraintCol])) {
        postBIdx <- missEidx+1
        postEIdx <- missEidx+nrow(posterior)
        dfPost <- data.frame(df1[postBIdx:postEIdx,])
        if (identical.norowname(posterior,dfPost,epsilon)) {
          if (countNulls(potentialResolver) < countNulls(miss)) {
            obj$resolver = potentialResolver
            obj$addConfidence()
            df1<-write2file_revised(obj,obj$missrn,wSize,cnt_first1,df1)
          }
        }
      }
      cnt = cnt + wSize
    }else {
      cnt = cnt+1
    }
    curr<- NULL
  }
  return(df1)
}

write2file_revised<-function(inObj, missidx,wSize,cnt,df1) {
  dfprior = inObj$prior
  dfresolver = inObj$resolver
  dfposterior = inObj$post

  dfall <- data.frame()
  dfall <- rbind(dfall,dfprior)
  dfall <- rbind(dfall,dfresolver)
  dfall <- rbind(dfall,dfposterior)

  L<-cnt+1-nrow(dfresolver)
  H<-L+nrow(dfresolver)-1
  df1[L:H,]<-dfresolver
  return(df1)
}

out_csv<-function(data_frame,direction_save) {
  name_file<-paste(direction_save,format(Sys.time(),"%d-%b-%Y %H.%M"),".csv")
  data_frame[is.na(data_frame)]<-""
  write.table(data_frame,name_file, append = FALSE, sep = ",", dec = ".",
              row.names = FALSE, col.names = TRUE)
}

ht = list() # this is the collection of all Hash objects within their offest

constraint_check<-function(constraintCol,data_frame){

  names_co<-colnames(data_frame)
  j=0
  for(i in 1:length(names_co))
  {
    if(names_co[i]==constraintCol){j<-i}
  }
  if(j==0){stop('The inserted constraint name is not valid!!')}

  else{return(j)  }}

countNulls <- function(indf) {
  a = 0
  for (i in 1:nrow(indf)) {
  a = a+sum(is.na(indf[i,]))
  }
  return (a)
}

hasNullRow <- function(indf) {
  for (i in 1:nrow(indf)) {
    if  ( all(is.na((indf[i,]))) ) {
      return (TRUE)
    }
  }
  return (FALSE)
}

identical.norowname <- function(df1,df2,epsilon) {
  rownames(df1) <- NULL
  rownames(df2) <- NULL
  sim = 0
  for (i in 1:ncol(df1))
  {
    if (identical(df1[i], df2[i]))
    {
      sim <- sim+1
    }
  }
  sim=sim/ncol(df1)
  if (sim>=epsilon)
    return(TRUE)
  else
    return(FALSE)
}

exactidentical.norowname <- function(df1,df2) {
  rownames(df1) <- NULL
  rownames(df2) <- NULL
  return (identical(df1,df2,1))
}

check=function(x) tryCatch(if(class(x) == 'logical') 1 else 1, error=function(e) 0)


slidewindow <- function(indf, constraintCol, windowS) {
  counter = 1
  tmpDF.rownames = NULL

  for (i in 1:(windowS+1)){ # i in 1:(windowS+1

    counter = i
    while (counter < nrow(indf)) {
      tmpDF = data.frame( indf[counter:((counter+windowS)-1 ), ])
      tmpDF = subset(tmpDF, select = -c(constraintCol)) # remove the constraint col for now, I don't know why I did it ???
      rownames(tmpDF) = NULL
      tmpOffsets = c(counter) # only the offset of the first window record will be used
      remove = c(0) # remove 0 and repeats from offests
      tmpOffsets = tmpOffsets[!tmpOffsets %in% remove]

      if (hasNullRow(tmpDF)==TRUE) {
      }else {
        objHT = HTCls$new(tmpDF, tmpOffsets)
        ht = addtoList(objHT, ht)
      }
      counter = counter + windowS
    }
  }
  return (ht)
}

addtoList <- function(inHTCls, inHT) {
  if (length(inHT) == 0) {
    tmpOffsets = c(inHTCls$offsets) # only the offset of the first window record will be used
    remove = c(0)
    tmpOffsets = tmpOffsets[!tmpOffsets %in% remove]
    inHTCls$offsets = tmpOffsets
    inHT[[1]] = inHTCls

  }
  else {
    for (k in 1: (length(inHT)) )  {

      a = inHT[[(k)]]$contentASkey
      b = inHTCls$contentASkey
      comparision = exactidentical.norowname (a, b)
      if (comparision) { # the windoe found in the Hash Table, thus offset of the founded window in Hashtable should be updated.
        curroffset = inHT[[(k)]]$offsets
        inHT[[(k)]]$offsets = appendOffset(curroffset, inHTCls$offsets)

        break; # the elemnt has been found so there is no need to check the rest of Hash Table.
      }else { #not existed added into the Hash Table and list size should be increased too

        if (k == (length(inHT))) { # we are at the end of the list and now it should be added into the list
          inHT[[(length(inHT)+1)]] = inHTCls

          break;  #because the new window has been added and there is no need for comparision anymore
        }
      }
    }
  }
  return (inHT)
}

appendOffset <- function(curoffsets, newoffsets ){
  result = as.vector(rbind(curoffsets,newoffsets))
  result = unique(result)
  remove = c(0)
  result = result[!result %in% remove]
  return (result)
}

chr <- function(n) { rawToChar(as.raw(n)) }






