###############################################################################
### Report generator for targeted HRMS (thermo) analysis
### Author Benoit Carreres
### Email:   Benoit.carreres@rd.nestle.com
versionNumber="1.0.0"
###############################################################################
Sys.setenv(LANG = "en")### to get messages in english (else it is mixed languages)
suppressWarnings(suppressMessages(require(openxlsx)))
suppressWarnings(suppressMessages(require(tidyr)))
suppressWarnings(suppressMessages(require(plyr)))
suppressWarnings(suppressMessages(require(dplyr)))
suppressWarnings(suppressMessages(require(rsq)))
suppressWarnings(suppressMessages(require(Rcpp)))
#############################################
############# Global variables   ############
#############################################
# ndigits.RSQ = 3
# ndigits.RTmin = 2
ndigits.Cc = 3
ccUnit = ""
### Preset styles
headerStyle <- createStyle(fontColour="#000000",textDecoration ="bold",halign="center",valign="center",fgFill="#EAEAEA",border="Bottom",borderStyle="thin",borderColour="#888888")
headerStyle.leftArea <- createStyle(halign="left")## to be applyed after normal header style, with stack.
leftArea <- createStyle(fgFill="#e1e1ff", halign="left", numFmt="TEXT")
resArea <- createStyle(fgFill="#f0ffe1", halign="center")
commentStyle <- createStyle(fgFill="#ffffe1", border=NULL, halign="center")
headerStyleMain <- createStyle(fontSize = 11, fontColour="#000000",textDecoration ="bold",halign="center",valign="center",fgFill="#EAEAEA",border="LeftRight",borderStyle="thin",borderColour="#888888")
headerStyleSimple <-  createStyle(fontColour="#000000",textDecoration ="bold",fgFill="#EAEAEA",border="TopBottom",borderStyle="medium",borderColour="#000000")
headerStyleSimple.meta <-  createStyle(fontColour="#000000",textDecoration ="bold",halign="left",valign="center",fgFill="#EAEAEA")
groupFirstCol <- createStyle(border = "Left", borderColour = "#888888")
greenText <- createStyle(fontColour = "#119911", fgFill = "#DEFFEF")
# greenBoldText <- createStyle(fontColour = "#119911",textDecoration ="bold" , fgFill = "#DEFFEF")
orangeText <- createStyle(fontColour = "#e06000", fgFill = "#DEFFEF")
redText <- createStyle(fontColour = "#dd1111", fgFill = "#DEFFEF")
blackText <- createStyle(fontColour = "#000000", fgFill = "#DEFFEF")
########################################################################################################################
########################################################################################################################
##################################################    Functions   ######################################################
########################################################################################################################
printManual <- function(){
  cat(readLines("README.txt"), sep="\n")
}
quoteNonNumeric <- function(v) {
  v[suppressWarnings(is.na(as.numeric(v)))] = paste0("\"",v[suppressWarnings(is.na(as.numeric(v)))],"\"")
  return(v)
}
smartConcat <- function(a,b,sep=" ") {
  sp1 = strsplit(a, "[ _]")[[1]]
  sp2 = strsplit(b, "[ _]")[[1]]
  ch.eq = which(!sp1 == sp2)[1]
  if(!is.na(ch.eq)) {
    a = paste(a,paste(sp2[ch.eq:length(sp2)], collapse=" "), collapse="_", sep=sep)
  } else {
    a = paste(unique(c(a,b)), collapse="\n")
  }
  return(a)
}
setDfCols2Num <- function(df, colNums){
  for(colNum in colNums) {
    df[which(df[,colNum] %in% c("N/A","N/F","NaN")),colNum] = NA
    df[,colNum] = as.numeric(df[,colNum])
  }
  return(df)
}
testUniformity <- function(df) {
  globalCompoundCounts = as.numeric(names(which.max(table(table(df$Compound.Name)))))
  outstandingCompound = names(which(table(df$Compound.Name) != globalCompoundCounts))
  if(length(outstandingCompound)>0){
    outstandingCompoundCounts = table(as.vector(df$Compound.Name[df$Compound.Name == outstandingCompound]))
    printableCounts = paste(names(outstandingCompoundCounts), outstandingCompoundCounts, sep = "=", collapse = "; ")
    warning("\n/!\\\n/!\\\tTable:",deparse(substitute(df)),"\n/!\\\tThe majority of compounds counts ",globalCompoundCounts," measurements.\n/!\\\tSome compounds have a different number of measurements than others: ",printableCounts,"\n/!\\")
  }
}
testLevelsConsistency <- function(df, type){### Works only on QC
  if(class(df) == "list") {
    invisible(lapply(df, function(df.notlist) testLevelsConsistency(df.notlist, type)))
  } else{
    ## Double filtration:
    ### Filter by sample Sample.name, which correspond to the part before the first "_".
    ### Filter by sample type, to avoid issues like QC overlapping revoveries.
    declaredMd = metaData2[grepl(paste0("^",paste(unique(df$Sample.name), collapse="|^")),metaData2[,"Raw.File.Name"]) & metaData2$Sample.Type %in% type,]
    ## Identify recovery sample and set sample level as pretended 0
    recov.sampleName = declaredMd$"Raw.File.Name"[which(declaredMd$Level=="QC_Rec")]
    df$Sample.Level[df$Sample.Raw.File.Name == recov.sampleName] <- 0
    #
    declaredLevels = unique(declaredMd$Level)
    declaredLevels = suppressWarnings(declaredLevels[order(as.numeric(gsub("^([0-9]+).*$","\\1",declaredLevels)))])### Warning expected if no number is in the level
    levelCounts = df %>% group_by(Compound.Name,Sample.name) %>% count(Sample.Level) %>% count
    w.counts = which(levelCounts$n!=length(declaredLevels))
    if(length(w.counts)>0){
      message("/!\\ Inconsistency in number of injection levels for the following compound data.")
      message("/!\\ \t Expected ",length(declaredLevels)," injection levels: ",paste(declaredLevels,collapse=", "))
      rownames(df) = NULL
      lapply(w.counts, function(w) message("/!\\ \t Found ",levelCounts$Compound.Name[w]," with ",levelCounts$n[w]," injection levels."))
      message("/!\\ See extract of data as found in the input:")
      print(lapply(w.counts, function(w) df[df$Compound.Name == levelCounts$Compound.Name[w] & df$Sample.name == levelCounts$Sample.name[w],c("Compound.Name","Sample.Raw.File.Name","Sample.Type","Sample.Level")]))
      stop("\n/!\\\n/!\\Errors in the consistency of injection levels, check above messages.\n/!\\\n\n")
    }
  }
}
## Calculations
QC.CV.byRow <- function(df) {
  injectionLevels = unique(gsub(" ","",metaData2[which(metaData2$Sample.Type %in% c("QC Std","Chk Std")),]$Level))
  injectionLevels = injectionLevels[order(as.numeric(gsub("^([0-9]+).*$","\\1",injectionLevels)))]### catching and sorting levels from batch table
  df.area = df[,grep("^Peak.Area",names(df), value=TRUE)]
  df.lvls = df[,grep("^Sample.Level", names(df), value=TRUE)]
  isnadf = is.na(df.area)
  if(any(isnadf)) {
    warning("\n/!\\\n/!\\\t",length(which(isnadf))," NA(s) found in peak area values. There must be missing measurements in the QC.","\n/!\\")
  }
  CVcols = do.call(rbind,sapply(1:nrow(df), function(i) {
    r = unlist(df.area[i,])
    rg = factor(df.lvls[i,], levels=sort(unique(unlist(df.lvls[i,]))))
    if(length(levels(rg)) != length(injectionLevels)) {
      stop("Inconsistency in number of injection levels for the following compound data.\n\tFound ",length(levels(rg)),": ",paste(levels(rg),collapse=", "),"\n\tExpected ",length(injectionLevels),": ",paste(injectionLevels,collapse=", "))
    }
    levelGroups = sapply(levels(rg), function(l) which(rg == l), simplify=FALSE)
    rg=as.numeric(as.vector(rg))
    model = lm(r~rg)
    CVs = sapply(levelGroups, function(cols) {
      if(length(cols)==1) {
        # warning("only 1 value found for the fortification level, returning CV of 0")### This is fine
        return(0)
      } else if(mean(r[cols])==0) {
        return(NA)
      } else {
        round((sd(r[cols])/mean(r[cols]))*100)
      }
    })
    names(CVs) = paste0("CV.",injectionLevels)
    CVs
  }, simplify=FALSE))
  df = cbind(df, CVcols)
  ### TODO merge both apply funcs
  if(length(injectionLevels) >1) {
    lin.stats = do.call(rbind,lapply(1:nrow(df), function(i) {
      y = unlist(df.area[i,])
      x = as.numeric(as.vector(df.lvls[i,]))
      # rg = as.numeric(as.vector(factor(df.lvls[i,], levels=sort(unique(unlist(df.lvls[i,]))))))
      model = lm(y~x)
      model.s = summary(model)
      r2 = model.s$adj.r.squared
      r2[is.nan(r2)] = NA
      d=((((y-model$coefficients[1])/model$coefficients[2])-x)/x)*100
      if(all(is.na(d))) {
        # if(is.nan(md)|is.na(md)|length(md)==0) message("max deviation is not a number or na: ", md, "\nCalculation: y=",paste0(y, collapse="_"), ", coef1=", model$coefficients[1], ", coef2=", model$coefficients[2], ", x=",paste0(x, collapse="_"),"\nx=",x,", y=",y,", df ",paste0(df[i,], collapse="; "))
        message("no calculable deviation: ",paste0(d, collapse=", "))
        message("Name: ", df$Compound.Name[i])
        md=NA
      } else {
        md = max(abs(d[!is.infinite(d)]), na.rm=TRUE)
      }
      # print(df[i,])
      # md[is.nan(md)] = NA
      cbind(RSQ=r2,max.deviation=md)
    }))
    df = cbind(df,lin.stats)
    df$QC.linearity = QC.lin(df)
  }
  return(df)
}
QC.check <- function(df) {
  apply(df,1,function(r){
    CV.num = as.numeric(r[grep("^CV",names(r))])
    r.lvls = r[grep("^Sample.Level", names(r), value=TRUE)]
    Inj.zero = gsub("^.+Level_","",names(r.lvls)[which(r.lvls==0)])
    if(any(r[grep("^Confirm", grep(paste0(Inj.zero,"$"),names(r), invert=TRUE, value=TRUE), value=TRUE)]=="NotFound")) {return("Not analyzable")}
    else if(all(is.na(CV.num))){return("Not analyzable")}
    else if(any(is.na(CV.num))){return("To check")}
    else if(any(CV.num > as.numeric(md.params$QCRepeatability.thres["value"]))) {return("To check")}
    else if(length(which(r[grep("^Confirm", names(r), value=TRUE)]=="Identified")) > length(r[grep("^Confirm", names(r), value=TRUE)])/4)  {return("Identified")}
    else (return("OK"))## Else it should be fine. Either all confirmed, few identified, or other citerias are low enough.
  })
}
checkAndModDataCol.injAreaRatios <- function(dataTable, colName, dataTable.withArea, is.screening){
  compColName = ifelse("ID" %in% names(dataTable) & "ID" %in% names(dataTable.withArea), "ID", "Compound.Name")
  m=match(dataTable[,compColName],dataTable.withArea[,compColName])
  df=do.call(rbind,lapply(split(dataTable.withArea[m,], dataTable.withArea[m,"Sample.name"]), function(dataTable.withArea.spl){
    ## Select the highest level that is not NA
    areaInj2plus = grep("^Peak.Area_Inj_[2-9][0-9]?$",colnames(dataTable.withArea.spl), value=TRUE)
    areaMin = head(areaInj2plus,1)## injections 2+, take the first (smallest)
    if(!all(is.na(dataTable.withArea.spl[areaInj2plus]))) {
      areaInj2plus = areaInj2plus[which(apply(dataTable.withArea.spl[areaInj2plus],2,function(x) all(!is.na(x))))]
    } else{
      message("Seems like only one injection level was provided for the sample. Check sample names for spelling mistakes.\n\t Table header:\n")
      print(head(dataTable.withArea.spl))
      stop("Seems like only one injection level was provided for the sample. Check sample names for spelling mistakes.\n")
    }
    areaMax = tail(areaInj2plus,1)## injections 2+, take the latest (biggest)
    levelMax = gsub("^Peak.Area_(Inj_[2-9][0-9]?)$","Sample.Level_\\1",areaMax)
    w.QN = dataTable.withArea.spl$"QN/QL" == "QN" | (is.na(dataTable.withArea.spl$"QN/QL") & !is.screening)
    areaRef = as.numeric(dataTable.withArea.spl$Peak.Area_Inj_1)
    w.isNum = suppressWarnings(!is.na(as.numeric(dataTable.withArea.spl[,colName])))
    ### testing if unpike area is <66% of max spiked area. If that is the case, the value is too large.
    inf.tooHigh = w.QN & w.isNum & areaRef > (0.66 * as.numeric(dataTable.withArea.spl[,areaMax]))
    w.tooHigh = which(inf.tooHigh)
    #### This is correct, the CC value given is not interpretable, hence simply wrong and can be <2*levelMax. It is possible to see a value that is lower than that reported with this ">" value.
    inf.tooHigh[w.tooHigh] = paste0(">",dataTable.withArea.spl[w.tooHigh,levelMax]*2)
    is.tooLow = w.QN & w.isNum & areaRef < (0.05 * as.numeric(dataTable.withArea.spl[,areaMin]))
    data.frame(ID=dataTable.withArea.spl$ID,inf.tooHigh=inf.tooHigh, is.tooLow=is.tooLow)
  }))
  m2 = match(dataTable$ID,df$ID)
  df = df[m2,]
  ### takes simple table as input and return the content of that column to be replaced
  w.tooHigh = which(df$inf.tooHigh!="FALSE")## because there are >x values
  dataTable[w.tooHigh,colName] <- df$inf.tooHigh[w.tooHigh]
  w.tooLow = which(df$is.tooLow)
  dataTable[w.tooLow,colName] <- ifelse(is.screening,"<STC","<LOQ")### WARNING this is not the right place to do that!!
  return(dataTable)
}
QC.lin <- function(df) {
  apply(df,1,function(r){
    if(as.numeric(r[grep("^RSQ$",names(r))]) < as.numeric(md.params$RSQthreshold["value"]) || is.na(r[grep("^RSQ$",names(r))])) {return("Check RSQ")}
    else if(all(!is.na(r[grep("^max.deviation$",names(r))])) & any(as.numeric(r[grep("^max.deviation$",names(r))]) > as.numeric(md.params$QCLinearity.thres["value"]))) {return("Check Deviation")} ###TODO give a dedicated threshold ### linearity to RSQ
    else (return("OK"))## Else it should be fine. Either all confirmed, few identified, or other citerias are low enough.
  })
}
Blank.check <- function(df.blk, df.qc) {### TODO The two sets of rules overlap totaly and maybe a combinaison needs to be done.
  df.blk$QC.median = apply(df.blk,1,function(r){
    median(as.numeric(unlist(df.qc[which(as.character(df.qc$Compound.Name) == r["Compound.Name"]),grep("^Peak\\.Area_Inj",names(df.qc))])))
  })
  df.blk$max.CarryOver = apply(df.blk,1,function(r){
    max(as.numeric(r[grep("^Peak\\.Area_Inj",names(df.blk))])) / as.numeric(r["QC.median"]) * 100
  })
  df.blk$CarryOver = ifelse(df.blk$max.CarryOver <= as.numeric(md.params$CarryOver.thres["value"]), "OK", "To check")
  df.blk$max.CarryOver[is.nan(df.blk$max.CarryOver)] = NA ### deviding by 0 results in NAN, hence an error message. Making it NA make the cell empty.
  df.blk$CarryOver[is.na(df.blk$CarryOver)] = "-"
  return(df.blk)
}
sample.linCc <- function(df, is.recovery=FALSE){
  if(is.recovery) {
    df$recov.Level.theoretical = df$Sample.Level_Inj_1
    df$Sample.Level_Inj_1 = 0
  }
  newDF = t(apply(df, 1, function(row) {
    x = as.numeric(as.vector(row[grep("^Sample.Level",names(row), value=TRUE)]))
    y = as.numeric(as.vector(row[grep("^Peak.Area",names(row), value=TRUE)]))
    model <- lm(y ~ x)
    Concentration = signif(unname(model$coefficients[1]/model$coefficients[2]), digits = ndigits.Cc)
    r2 = as.vector(rsq(model))
    Concentration[is.nan(Concentration)] = NA
    r2[is.nan(r2)] = NA
    ### Check validity of values
    conf = row[grep("^Confirm_",names(row), value=TRUE)]
    conf.unspiked = row[names(conf)[which(x == 0)]]
    conf.spiked = row[names(conf)[which(x > 0)]]
    ### Recovery
    if(is.recovery) {
      QC.recovery = (Concentration / as.numeric(row["recov.Level.theoretical"])) * 100
      QC.recovery[is.nan(QC.recovery)] = NA
      if(any(conf.spiked == "NotFound" | conf.unspiked == "NotFound")) {## setting values to NA because not valid
        Concentration = NA
        QC.recovery = NA
      }
      return(c(RSQ=r2,recov.Level.calculated=Concentration, QC.recovery=QC.recovery))
    }
    ##### For quantification mode
    if(any(conf.spiked == "NotFound")) {
      Concentration = "Not analyzable"
      Comment = "Not found in spiked sample(s)"
      ccCol = "#dd1111"#"red"
    } else if(all(conf.unspiked == "Identified")){
      Comment = "Identified in unspiked sample(s)"
      ccCol = "#e06000"#"orange"
    } else if(all(conf.unspiked == "Confirmed")){
      Comment = ""
      ccCol = "#119911"#"green"
    } else if(any(conf.unspiked == "Confirmed")){
      Comment = "Identified in unspiked sample(s)"
      ccCol = "#119911"#"green"
    } else if(any(conf.spiked == "Confirmed")){
      Concentration = "Not detected"
      Comment = ""
      ccCol = "#119911"#"green"
    } else {
      Concentration = "Not detected"
      Comment = "Identified in spiked sample(s)"
      ccCol = "#e06000"#"orange"
    }
    if(!is.na(r2) & r2 < as.numeric(md.params$RSQthreshold["value"])) {Comment = ifelse(Comment=="", "Check RSQ", paste0(Comment, "; Check RSQ"))}
    return(c(RSQ=r2,Concentration=Concentration,ccCol=ccCol,Comment=Comment))
  }))
  return(cbind(df,newDF))
}
sample.screeningInfo <- function(df, is.recovery=FALSE){
  w.sl = grep("Sample\\.Level_",colnames(df))
  w.pa = grep("Peak\\.Area_",colnames(df))
  sl.cs = colSums(df[,w.sl])### Could also base on the first row
  w.z  = as.vector(which(sl.cs == 0))### Could also base on the first row
  w.min  = as.vector(which(sl.cs == min(sl.cs[sl.cs != 0])))
  confCols = df[,grep("^Confirm_",colnames(df), value=TRUE)]
  if(length(w.z)>1) message("More than one unspiked samples detected, reporting the worst available value")
  if(length(w.min)>1) message("More than one spiked samples detected, reporting using mean value")
  area.zmax = matrixStats::rowMaxs(as.matrix(df[,w.pa][,w.z]))### as.matrix is to avoid making conditions when only 1 col -> vector
  area.mm = rowMeans(as.matrix(df[,w.pa][,w.min]))
  area.ratio = area.zmax / area.mm * 100### Percentage prefered by operators, threshold also described as percentage
  area.ratio[is.nan(area.ratio) | is.infinite(area.ratio)] = NA## Set NAN to NA
  res = ifelse(area.ratio > df$"Cut-off", "Suspect", "Negative")
  res[is.na(res)] = "Not analyzable"
  if(is.recovery){
    res[res=="Suspect"] = "OK"
    res[res=="Negative"] = "To check"
    qcsCol = rep("#119911", nrow(df))##green
    qcsCol[res=="To check"] = "#e06000"#"orange"
    qcsCol[res=="Not analyzable"] = "#dd1111"#"red"
    qcsCol[rowSums(confCols == "Identified")/ncol(confCols) >1/3] =  "#e06000"#"orange"
    newCols = data.frame(QC.Area.ratio = area.ratio, QC.Screening = res,qcsCol = qcsCol)
  } else {
    newCols = data.frame(Area.ratio = area.ratio, Screening = res)
  }
  df = cbind(df,newCols)
  return(df)
}
samples.analyze <- function(data.sample, is.recovery=FALSE) {
  if(any(colnames(data.sample) == "Test.portion.factor")){
    w = !is.na(data.sample$Test.portion.factor)
    data.sample$Sample.Level[w] = data.sample$Sample.Level[w] * data.sample$Test.portion.factor[w]
    data.sample$LOQ[w] = data.sample$LOQ[w] * data.sample$Test.portion.factor[w]
  }
  df = pivot_wider(data.sample, id_cols=c(Compound.Name,CAS.Number,Family,Sample.name,LOQ,"Cut-off","QN/QL"), names_from = ReplicateName, values_from = c(Sample.Level , Peak.Area, Confirm, "RT.(min)", Detected.Mass, "m/z.Delta","F1.m/z","F2.m/z","F3.m/z","F4.m/z","F5.m/z",F1.Delta.MZ,F2.Delta.MZ,F3.Delta.MZ,F4.Delta.MZ,F5.Delta.MZ,NbFrag))
  if(nrow(df)==0) {
    warning("Empty sample, most likely an empty row in the original CSV.")
    return()
  }
  if(all(is.na(df$"Cut-off")) | all(df$"Cut-off"==0)) df = df[,colnames(df)!="Cut-off"]
  df = merge(df, QC.w[,c("Compound.Name","RT.(min)_Mean")], by="Compound.Name")
  df$RT_QC = df$"RT.(min)_Mean"
  df$"RT.(min)_Mean" = NULL
  df = suppressWarnings(sample.linCc(df, is.recovery=is.recovery))
  df$RSQ = as.numeric(as.vector(df$RSQ))### Forcing here because it is systematically changed back to factor
  ### Screening mode
  if(is.screening) {
    df = sample.screeningInfo(df, is.recovery=is.recovery)
  }
  return(df)
}
###
groupsColOrder = function(df) {
  allColNames = names(df)
  groupedCols = list()
  ### Static columns, manually set
  preferedOrder = list(
    "Static_information" = c(
      grep("^ID$",allColNames),
      grep("^Batch$",allColNames),
      grep("^Sample\\.name$",allColNames),
      grep("^Compound\\.Name$",allColNames),
      grep("^CAS\\.Number$",allColNames),
      grep("^Family$",allColNames),
      grep("^Cut-off$",allColNames),
      grep("^LOQ$",allColNames),
      grep("^QN/QL$",allColNames)
    ),
    "Results" = c(
      grep("^Screening$",allColNames),
      grep("^Area\\.ratio$",allColNames),
      grep("^Concentration$",allColNames),
      grep("^SAP$",allColNames)
    ),
    "Data control" = c(
      grep("^QC\\.Screening$",allColNames),
      grep("^QC\\.Area\\.ratio$",allColNames),
      grep("^QC\\.recovery$",allColNames),
      grep("recov\\.Level\\.calculated",allColNames),
      grep("recov\\.Level\\.theoretical",allColNames),
      grep("^QC\\.repeatability$",allColNames),
      grep("^QC\\.linearity$",allColNames),
      grep("^CV",allColNames),
      grep("^RSQ$",allColNames),
      grep("^max\\.deviation$",allColNames),
      grep("^CarryOver$",allColNames),
      grep("^max\\.CarryOver$",allColNames)
    ),
    "Identification_criteria" = grep("^Confirm",allColNames),
    c(grep("^Comment$",allColNames)),### No category name for that
    "Fortification_level" = grep("^Sample.Level",allColNames),
    "Area" = c(grep("^QC\\.median$",allColNames),which((!grepl("^CV",allColNames)) & grepl("^Peak\\.Area_[^(cv)]",allColNames))),
    "RT" = which(grepl("^RT",allColNames)),
    "Number of fragments" = grep("^NbFrag",allColNames)
  )
  ### Dynamic columns, to be automatically set: MZs
  mzCols = grep("Detected.Mass|m\\/z|MZ",allColNames)
  spikingGroupNames = gsub(".*?_(.+)$","\\1",allColNames)
  mzColsGroupNames = spikingGroupNames[mzCols]# selecting the cols
  mzColsGroupList = sapply(unique(mzColsGroupNames),function(g) mzCols[which(mzColsGroupNames == g)], USE.NAMES=TRUE, simplify=FALSE)## listing each col IDs per naming groups
  mzColsGroupList = sapply(names(mzColsGroupList), function(g) {
    cs = allColNames[mzColsGroupList[[g]]]
    ct = gsub(paste0("_",g,"$"),"",cs)
    mzColsGroupList[[g]][c(which(ct=="Detected.Mass"),which(ct=="m/z.Delta"),which(ct=="NbFrag"),which(!ct%in%c("m/z.Delta","NbFrag","Detected.Mass")))]
  }, USE.NAMES=TRUE, simplify=FALSE)
  ###
  groupedCols = c(preferedOrder,"mz"=mzColsGroupList)### Will add mz in front of the mz groups
  ### If anything else
  leftOverCols = which(!1:length(allColNames) %in% unlist(groupedCols))
  if(length(leftOverCols)>0) {
    warning("Found unspecified columns, adding to the end of the table: ", paste0(leftOverCols, collapse=", ", "---"))
    groupedCols = c(groupedCols,unspeciFiedCols = leftOverCols)
  }
  groupedCols
}
addReplicateName <- function(df, splitVector=1) {
  ### Base to rename injections in the order of appearance
  df$ReplicateName = ""
  df = do.call(rbind,lapply(split(df, splitVector), function(spl) {
    spl$ReplicateName = paste0("Inj_",as.numeric(factor(spl$Sample.Raw.File.Name, levels=unique(spl$Sample.Raw.File.Name))))
    return(spl)
  }))
  return(df)
}
renameHeader <- function(header) {#### Last modifications of headers to better fit operators likings
  ### Important adjustments
  header = gsub("^Fortification_level$",paste0("Fortification level (",md.params$ccUnit["value"],")"),header)
  header = gsub("^RT$","RT (min)",header)
  header = gsub("^Concentration$",paste0("Concentration (",md.params$ccUnit["value"],")"),header)
  header = gsub("^(CV.*)$","\\1 (%)",header)
  header = gsub("^max\\.CarryOver$", "max.CarryOver (%)",header)
  header = gsub("^max\\.deviation$", "max.deviation (%)",header)
  header = gsub("Area\\.ratio$", "Area.ratio (%)",header)
  header = gsub("^QC.recovery$", "QC.recovery (%)",header)
  header = gsub("^recov.Level.calculated$",paste0("Calculated (",md.params$ccUnit["value"],")"),header)
  header = gsub("^recov.Level.theoretical$",paste0("Theoretical (",md.params$ccUnit["value"],")"),header)
  ### Special cases
  if(is.screening) {### if screening mode
    header = gsub("^LOQ$",paste0("STC (",md.params$ccUnit["value"],")"),header)
    header = gsub("^Cut-off$","Cut-off (%)",header)
  } else {### if not screening mode, then quantitative mode (default)
    header = gsub("^LOQ$",paste0("LOQ (",md.params$ccUnit["value"],")"),header)
  }
  ### Aesthetics change of header
  header = gsub("m\\.z","m/z",header)
  header = gsub("mz","m/z",header)
  header = gsub("\\."," ",header)
  header = gsub("_"," ",header)
  return(header)
}
### splitting header in two lines and merging grouped cols
headerSplitRename <- function(df){### Function to generate splitted headers in in the excel sheet. Columns should already be ordered using the function "groupsColOrder"
  df.colOrder = groupsColOrder(df)
  mainHeader = unlist(lapply(1:length(df.colOrder),function(gn) {
    colN = NULL
    if(length(df.colOrder[[gn]])==1) {
      colN = names(df.colOrder)[gn]
    } else if(length(df.colOrder[[gn]])>1) {
      colN = names(df.colOrder)[gn]
      colN = c(colN,rep("",length(df.colOrder[[gn]])-1))
    }
    colN
  }))
  secHeader = unlist(sapply(1:length(df.colOrder),function(gn) {
    gName = names(df.colOrder)[gn]
    cols = df.colOrder[[gn]]
    if(gName == ""){### if no group name, keep the colname
      colN=names(df)[cols]
    } else {### if group name, remove either prefix or sufix
      colN = names(df)[cols]
      if(grepl("^mz",gName)) {
        colN = gsub("(.*?)_.+$","\\1",colN)
      } else {
        colN = gsub(".*?_(.+)$","\\1",colN)
      }
    }
    colN
  }))
  mainHeader = renameHeader(mainHeader)
  secHeader = renameHeader(secHeader)
  rbind(mainHeader,secHeader)
}
buildSheetbyGroups = function(df, wb, sht, metaDataTable=NULL) {
  message("Sheet: ", sht)
  sr=2
  if(sht=="Report") sr=sr+1
  ### Preliminary step, remove special columns that should not be reported
  if(any(names(df) == "ccCol")){
    df.ccCol = df$ccCol
    df.ccCol = factor(df.ccCol)
    df = df[,!names(df) %in% "ccCol"]
  }
  if(any(names(df) == "qcsCol")){
    df.qcsCol = df$qcsCol
    df.qcsCol = factor(df.qcsCol)
    df = df[,!names(df) %in% "qcsCol"]
  }
  ### First, order columns according to the group col ordering
  df = df[,as.vector(unlist(groupsColOrder(df)))]
  df = setDfCols2Num(df,grep("F[0-9]+\\.m/z_|Detected\\.Mass|Delta", names(df)))

  addWorksheet(wb, sheet = sht)
  writeData(wb, sheet = sht, x = df, headerStyle = headerStyle, startRow = sr)
  addStyle(wb, sheet = sht, createStyle(halign="center"), rows=sr:(sr+nrow(df)), cols=1:ncol(df), gridExpand = TRUE, stack = TRUE)### defualt style for everything
  #############
  ## Columns ##
  #############
  ### Groups
  df.colOrder = groupsColOrder(df)### at this point it should be already ordered and this is only useful to group cols
  ### Group separation lines
  df.groupFirstCol = unlist(sapply(df.colOrder, function(x) if(length(x)>0) min(x)))### Logically, cols should have been ordered prior this call
  addStyle(wb, sheet = sht, style = groupFirstCol, cols = df.groupFirstCol, rows = sr, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet = sht, style = groupFirstCol, cols = df.groupFirstCol, rows = (1:nrow(df)+sr), gridExpand = TRUE, stack = TRUE)
  ### Areas style
  addStyle(wb, sheet = sht, style = headerStyle.leftArea, rows = sr, cols = df.colOrder[[1]], stack = TRUE)
  addStyle(wb, sheet = sht, style = leftArea, cols = df.colOrder[[1]], rows = (1:nrow(df)+sr), gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet = sht, style = resArea, cols = unlist(df.colOrder[2:5]), rows = (1:nrow(df)+sr), gridExpand = TRUE, stack = TRUE)
  #################
  ### Coloring text
  #################
  ###### Complex Color conditions, use integrated column
  if(exists("df.ccCol")){### based on preset color as a data column
    ccColN = which(names(df) == "Concentration")
    for (col in levels(df.ccCol)) {### +sr on the row to skip the headers
      addStyle(wb, sheet = sht, createStyle(fontColour = col, halign="center"), rows = which(df.ccCol==col)+sr, cols = ccColN, gridExpand = FALSE, stack = TRUE)
    }
    lapply(which(names(df) == "Screening"), function(colN){
      for (col in levels(df.ccCol)) {### +sr on the row to skip the headers
        addStyle(wb, sheet = sht, createStyle(fontColour = col, halign="center"), rows = which(df.ccCol==col)+sr, cols = colN, gridExpand = FALSE, stack = TRUE)
      }
    })
  }
  if(exists("df.qcsCol")){### based on preset color as a data column
    qcsColN = which(names(df) == "QC.Screening")
    for (col in levels(df.qcsCol)) {### +sr on the row to skip the headers
      addStyle(wb, sheet = sht, createStyle(fontColour = col, halign="center"), rows = which(df.qcsCol==col)+sr, cols = qcsColN, gridExpand = FALSE, stack = TRUE)
    }
  }

  ###### Basic value/word based coloring
  lapply(which(names(df) == "QC.recovery"), function(col){
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0("> ",as.numeric(md.params$recov.Thres["value2"])), style = redText)###condition
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0("< ",as.numeric(md.params$recov.Thres["value1"])), style = redText)###condition
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "between", rule = c(as.numeric(md.params$recov.Thres["value1"]), as.numeric(md.params$recov.Thres["value2"])), style = greenText)###condition
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "contains", rule="Not tested", style = blackText)
  })
  lapply(which(names(df) == "RSQ"), function(col){
    col = which(names(df) == "RSQ")
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0(">= ",as.numeric(md.params$RSQthreshold["value"])), style = greenText)###condition
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0("< ",as.numeric(md.params$RSQthreshold["value"])), style = redText)###condition
    addStyle(wb, sheet = sht, createStyle(numFmt = "0.000"), rows = (1:nrow(df))+sr, cols = col, stack = TRUE)
  })
  ### RT
  lapply(df.colOrder$RT, function(col){
    addStyle(wb, sheet = sht, createStyle(numFmt = "0.00"), rows = (1:nrow(df))+sr, cols = col, gridExpand = TRUE, stack = TRUE)
  })
  lapply(which(names(df) == "max.CarryOver"), function(col){
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0(">= ",as.numeric(md.params$CarryOver.thres["value"])), style = redText, gridExpand = FALSE)###condition
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0("< ",as.numeric(md.params$CarryOver.thres["value"])), style = greenText, gridExpand = FALSE)###condition
    addStyle(wb, sheet = sht, createStyle(numFmt = "0"), rows = (1:nrow(df))+sr, cols = col, stack = TRUE)
  })
  lapply(which(names(df) == "max.deviation"), function(col){
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0(">= ",as.numeric(md.params$QCRepeatability.thres["value"])), style = redText, gridExpand = FALSE)###condition
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0("< ",as.numeric(md.params$QCRepeatability.thres["value"])), style = greenText, gridExpand = FALSE)###condition
    addStyle(wb, sheet = sht, createStyle(numFmt = "0.0"), rows = (1:nrow(df))+sr, cols = col, stack = TRUE)
  })
  lapply(grep("^CV\\.",names(df)), function(col){
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0("> ",as.numeric(md.params$QCLinearity.thres["value"])), style = redText, gridExpand = FALSE)###condition
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0("<= ",as.numeric(md.params$QCLinearity.thres["value"])), style = greenText, gridExpand = FALSE)###condition
    addStyle(wb, sheet = sht, createStyle(numFmt = "0"), rows = (1:nrow(df))+sr, cols = col, gridExpand = FALSE, stack = TRUE)
  })
  lapply(grep("^Peak\\.Area|QC\\.median",names(df)), function(col){
    addStyle(wb, sheet = sht, createStyle(numFmt = "SCIENTIFIC"), rows = (1:nrow(df))+sr, cols = col, gridExpand = TRUE, stack = TRUE)
  })
  lapply(which(names(df) == "QC.repeatability"), function(col){
    conditionalFormatting(wb, sht, rows = (1:nrow(df))+sr, cols = col, type = "contains", rule="OK", style = greenText, gridExpand = FALSE, stack = TRUE)
    conditionalFormatting(wb, sht, rows = (1:nrow(df))+sr, cols = col, type = "contains", rule="To check", style = orangeText, gridExpand = FALSE, stack = TRUE)
    conditionalFormatting(wb, sht, rows = (1:nrow(df))+sr, cols = col, type = "contains", rule="Identified", style = orangeText, gridExpand = FALSE, stack = TRUE)
    conditionalFormatting(wb, sht, rows = (1:nrow(df))+sr, cols = col, type = "contains", rule="Not analyzable", style = redText, gridExpand = FALSE, stack = TRUE)
  })
  lapply(which(names(df) == "QC.linearity"), function(col){
    conditionalFormatting(wb, sht, rows = (1:nrow(df))+sr, cols = col, type = "contains", rule="OK", style = greenText, gridExpand = FALSE, stack = TRUE)
    conditionalFormatting(wb, sht, rows = (1:nrow(df))+sr, cols = col, type = "contains", rule="Check RSQ", style = orangeText, gridExpand = FALSE, stack = TRUE)
    conditionalFormatting(wb, sht, rows = (1:nrow(df))+sr, cols = col, type = "contains", rule="Check Deviation", style = orangeText, gridExpand = FALSE, stack = TRUE)
  })
  lapply(which(names(df) == "CarryOver"), function(col){
    conditionalFormatting(wb, sht, rows = (1:nrow(df))+sr, cols = col, type = "contains", rule="OK", style = greenText, gridExpand = FALSE, stack = TRUE)
    conditionalFormatting(wb, sht, rows = (1:nrow(df))+sr, cols = col, type = "contains", rule="To check", style = orangeText, gridExpand = FALSE, stack = TRUE)
  })
  ############
  ### peaksInfo formating
  lapply(grep("F[0-9]+\\.m/z_|Detected\\.Mass", names(df)), function(col){
    addStyle(wb, sheet = sht, createStyle(numFmt = "0.0000"), rows = (1:nrow(df))+sr, cols = col, gridExpand = TRUE, stack = TRUE)
  })
  lapply(grep("Delta", names(df)), function(col){
    addStyle(wb, sheet = sht, createStyle(numFmt = "0.0"), rows = (1:nrow(df))+sr, cols = col, gridExpand = TRUE, stack = TRUE)
  })
  ############ SCREENING
  lapply(which(names(df) %in% c("Area.ratio")), function(col){
    cutoff.colLetter = toupper(letters[which(names(df) == "Cut-off")])
    ar.colLetter = toupper(letters[col])
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0(ar.colLetter,1+sr," > ",cutoff.colLetter,1+sr), style = redText, gridExpand = FALSE)###condition
    conditionalFormatting(wb, sht, cols = col, rows = (1:nrow(df))+sr, type = "expression", rule=paste0(ar.colLetter,1+sr," <= ",cutoff.colLetter,1+sr), style = greenText, gridExpand = FALSE)###condition
  })
  lapply(which(names(df) %in% c("Area.ratio","QC.Area.ratio","QC.recovery")), function(col){
    ar.colLetter = toupper(letters[col])
    addStyle(wb, sheet = sht, createStyle(numFmt = "0"), rows = (1:nrow(df))+sr, cols = col, gridExpand = FALSE, stack = TRUE)
  })
  #############
  ## Headers ##
  #############
  df.newHeaders = headerSplitRename(df)### this one will get it's own group order to make it independent and avoid mistakes
  ###
  #### Build the main header, above the normal header.
  if(!all(df.newHeaders[1,]=="")){
    df.newHeaders[1,]
    writeData(wb, sheet = sht, x = t(df.newHeaders[1,]), startRow=sr-1, colNames = FALSE)
    addStyle(wb, sheet = sht, style = headerStyleMain, rows = sr-1, cols = 1:ncol(df), gridExpand = TRUE, stack = TRUE)
    setRowHeights(wb, sht, rows = sr-1, heights = 25)
    ignore = lapply(df.colOrder[lengths(df.colOrder)>1], function(x) mergeCells(wb, sheet=sht, cols=x, rows=sr-1))
    writeData(wb, sheet = sht, x = t(df.newHeaders[2,]), startRow=sr, colNames = FALSE)
  }
  setRowHeights(wb, sht, rows = sr, heights = 15)
  ## merge group cols
  setColWidths(wb, sht, cols=1:ncol(df), widths = "auto", ignoreMergedCells = TRUE, hidden = colnames(df) == "ID")
  ###################
  ############ Setting cold width for few annoying ones
  w.1 = colnames(df) %in% "Compound.Name"
  w.2 = colnames(df) %in% c("CAS.Number","Family")
  w.3 = colnames(df) %in% "Batch"
  w.4 = colnames(df) %in% "QN/QL"
  if(any(w.1)) setColWidths(wb, sht, cols=which(w.1), widths = 32)
  if(any(w.2)) setColWidths(wb, sht, cols=which(w.2), widths = 15)
  if(any(w.3)) setColWidths(wb, sht, cols=which(w.3), widths = ifelse(length(args)==1,0,2))
  if(any(w.4)) setColWidths(wb, sht, cols=which(w.4), widths = 3)
  ## Add column filters, at the location of the right header (StartRow)
  if(sht != "Report") {
    addFilter(wb, sht, row = sr, cols = 1:ncol(df))
    freezePane(wb, sht, firstActiveRow = sr+1, firstActiveCol = max(df.colOrder[[1]])+1)
  } else {
    addFilter(wb, sht, row = sr, cols = 1:(ncol(df)+1))### Considering the future addition of the comment column.
    freezePane(wb, sht, firstActiveRow = sr+1)
  }
  ### set static info with numbers centered to middle
  addStyle(wb, sheet = sht, style=createStyle(halign="center"), cols=grep("^Cut-off|^STC|^LOQ|QN\\/QL",names(df)), rows=sr:(nrow(df)+sr), stack=TRUE,gridExpand=TRUE)
  return(sr)
}
addBatchColorStyle <- function(wb, sht, df, sr=2) {
  for(b in seq_along(unique(df$Batch))) {
    w = which(df$Batch == b)
    addStyle(wb, sheet = sht, createStyle(fgFill=batch.color[b]), rows = w+(sr-1), cols = which(colnames(df)=="Batch"),stack=TRUE)
  }
}
########################################################################################################################
########################################################################################################################
#################################################   Data Analysis   ####################################################
########################################################################################################################
message("-------------------------\n---- Begining report ----\n-------------------------\n")
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n\tYou may drag and drop the correctly formated xlsx file onto this executable file.\n\n",printManual(), call.=FALSE)
}
message("Input file(s) = ", paste(args, collapse=", "))
file=args[1]
if(length(args)==2) file=smartConcat(gsub("\\.xlsx$","",args[1]),gsub("\\.xlsx$","",args[2]), sep="_")
outFile = gsub("\\.xlsx$","",file)
outFile = paste0(outFile,"_report",if(length(args)>1) "_m",".xlsx")
message("Output file = ", outFile,"\n")
m.metaData1 = data.frame()
m.metaData2 = data.frame()
m.is.screening = NULL
m.Blank.w = data.frame()
m.QC.w = data.frame()
m.Samples.w = data.frame()
m.data.samplesRecov.w = data.frame()
m.Simple.table = data.frame()
m.allData = data.frame()
l.metaData1 = list()

for(i in 1:length(args)) {
  inFile = args[i]
  message("----- Input & Output ----")
  message("---- Processing file ----\n")
  allData = read.xlsx(inFile, sheet="allData")
  ##################    EXTRACTING METADATA INFORMATION     ######################
  message("Processing metadata")
  metaData1 = read.xlsx(inFile, sheet="metadata")
  ### Metadata 1
  ##### Testing conformity of required infor from metaData1
  md.params = list(
    ccUnit = c(name = "Concentration Unit", value = metaData1[grep("^Concentration Unit",metaData1[,1]),2]),
    CarryOver.thres = c(name = "Carry over max", value = as.numeric(metaData1[grep("^Carry over max",metaData1[,1]),2])),
    QCRepeatability.thres = c(name = "QC Solvent Repeatability", value = as.numeric(metaData1[grep("^QC Solvent Repeatability",metaData1[,1]),2])),
    QCLinearity.thres = c(name = "QC Solvent Linearity", value = as.numeric(metaData1[grep("^QC Solvent Linearity",metaData1[,1]),2])),
    RSQthreshold = c(name = "Linearity", value = as.numeric(metaData1[grep("^Linearity",metaData1[,1]),2])),
    # Scrg.thres = c(name = "Screening Cut-off", value = as.numeric(metaData1[grep("^Screening Cut-off",metaData1[,1]),2])),
    QtOrScrg = c(name = "Screening or Quantitative mode", value = metaData1[grep("Screening or Quantitative mode",metaData1[,1], fixed=TRUE),2]),
    recov.Thres = c(name = "QC sample recovery", value = as.numeric(unlist(strsplit(metaData1[grep("^QC sample recovery",metaData1[,1]),2], split="-"))))### because 2 values, this one has elecments called "value1" and "value2"
  )
  metaData1 = metaData1[metaData1[,1]!="Screening Cut-Off Unit",]
  metaData1[grep("^Linearity",metaData1[,1]),2]=as.numeric(metaData1[grep("^Linearity",metaData1[,1]),2])### force numeric to avoid formating issue: 0.97999999999999998

  is.screening = grepl("screening",md.params$QtOrScrg["value"], ignore.case=TRUE)
  if(any(lengths(md.params)<2)) {message(printManual(), "\n\nCannot identify the required metadata: \"", md.params[[which(lengths(md.params)<2)]]["name"], "\"\nPlease check section \"Metadata requirements\" of the manual to resolve the issue.\n\n")}
  ##### Now extrating needed info
  metaData1 = rbind(c("Reporting version",versionNumber), metaData1)
  ###### OUTPUT for merging files
  m.metaData1 = merge(m.metaData1, metaData1)
  l.metaData1[[inFile]] = metaData1
  ###### Reformat from 2 to 4 column format
  metaData1.1 = metaData1[1:which(metaData1=="Operator"),]
  metaData1.2 = metaData1[(which(metaData1=="Operator")+1):nrow(metaData1),]
  metaData1 = data.frame(matrix(NA, nrow = max(nrow(metaData1.1),nrow(metaData1.2)), ncol = 4))
  metaData1[1:nrow(metaData1.1),1:2] = metaData1.1
  metaData1[1:nrow(metaData1.2),3:4] = metaData1.2
  ### Metadata2
  metaData2 = read.xlsx(inFile, sheet="Batch")
  metaData2$Sample.Type[metaData2$Level  %in% c("0xLOQ","0xSTC") & metaData2$Sample.Type %in% c("QC Std", "Chk Std")] <- "Solvent"

  ##################    EXTRACTING METADATA INFORMATION     ######################
  message("Processing data")
  allData$Compound.Name = factor(allData$Compound.Name)
  allData$Sample.name = gsub("^(.+?)_.+$","\\1",allData$Sample.Raw.File.Name)
  allData$Sample.Type[grep("QC_0xLOQ",allData$Sample.name)] = "Solvent"
  allData$Sample.Level = as.vector(allData$Sample.Level)
  allData$Sample.Level[allData$Sample.Level=="N/A"] = 0
  allData$Sample.Level = as.numeric(allData$Sample.Level)
  #
  allData$"RT.(min)" = as.vector(allData$"RT.(min)")
  allData$"RT.(min)"[allData$"RT.(min)" == "N/F"] <- NA
  allData$"RT.(min)" = as.numeric(allData$"RT.(min)")
  ### MZ delta has (ppm) within the value, to remove. Also, "0.123" is written "-.123"
  allData$"m/z.Delta" = gsub("-\\.","-0.",allData$"m/z.Delta")
  allData$"m/z.Delta" = gsub("^\\.","0.",allData$"m/z.Delta")
  allData$"m/z.Delta" = gsub(" \\(ppm\\)","",allData$"m/z.Delta")
  ############################       Data split        ###########################
  message("Adding replicate name to Blanks")
  data.blanks = addReplicateName(allData[allData$Sample.Type == "Solvent",])
  ### Separating samples, qc and blanks based on type. Then the replicates will be based on naming pattern ("ReplicateName")
  message("Adding replicate name to QC")
  data.QC = addReplicateName(allData[allData$Sample.Type %in% c("QC Std","Chk Std"),])
  message("Adding replicate name to samples")
  data.samples = allData[!allData$Sample.Type %in% c("Solvent","QC Std","Chk Std") ,]
  ### Samples need extra separation, Sample name needs to be generated properly.
  data.samples$Sample.name = gsub("^(.+?)_.+$","\\1",data.samples$Sample.name)
  data.samples = addReplicateName(data.samples, data.samples$Sample.name)
  message("---- Verifyring data ----\n")
  message("-- uniformity\n")
  testUniformity(data.QC)
  testUniformity(data.blanks)
  testUniformity(data.samples)
  message("-- level consistency\n")
  message("- QC\n")
  testLevelsConsistency(data.QC, type=c("QC Std","Chk Std"))
  message("- Blank\n")
  testLevelsConsistency(data.blanks, type="Solvent")
  ## samples are tested bellow
  message("-------------------------\n")
  message("---- Processing data ----\n")
  message("-      QC table dims: ",paste0(dim(data.QC),collapse=", "),"\n")
  message("-  Blanks table dims: ",paste0(dim(data.blanks),collapse=", "),"\n")
  message("- Samples table dims: ",paste0(dim(data.samples),collapse=", "),"\n")
  #################################
  ##########     QC     ###########
  #################################
  message("-- Processing QC\n")
  QC.w = pivot_wider(data.QC, id_cols=c(Compound.Name,Family,CAS.Number), names_from = ReplicateName, values_from = c(Sample.Level , "RT.(min)", Peak.Area, Confirm, Detected.Mass, "m/z.Delta","F1.m/z","F2.m/z","F3.m/z","F4.m/z","F5.m/z",F1.Delta.MZ,F2.Delta.MZ,F3.Delta.MZ,F4.Delta.MZ,F5.Delta.MZ,NbFrag))
  QC.w = QC.CV.byRow(QC.w)
  QC.w$"RT.(min)_Mean" = rowMeans(QC.w[,grep("^RT.\\(min\\)",names(QC.w), value=TRUE)], na.rm=TRUE)
  QC.w$QC.repeatability = QC.check(QC.w)
  QC.w = QC.w[order(factor(QC.w$QC.repeatability, levels=c("To check","NotFound","Identified","OK"))),]
  QC.w$Batch = i
  #################################
  ##########   Blanks   ###########
  #################################
  message("-- Processing Blanks\n")
  Blank.w = pivot_wider(data.blanks, id_cols=c(Compound.Name,Family,CAS.Number), names_from = ReplicateName, values_from = c(Sample.Level , "RT.(min)", Peak.Area, Confirm, Detected.Mass, "m/z.Delta","F1.m/z","F2.m/z","F3.m/z","F4.m/z","F5.m/z",F1.Delta.MZ,F2.Delta.MZ,F3.Delta.MZ,F4.Delta.MZ,F5.Delta.MZ,NbFrag))
  Blank.w = Blank.check(Blank.w, QC.w)
  Blank.w = Blank.w[order(factor(Blank.w$CarryOver, levels=c("To check","OK"))),]
  Blank.w$Batch = i
  #################################
  ##########   Sample   ###########
  #################################
  message("-- Processing Samples\n")
  data.samples.l = split(data.samples, data.samples$Sample.name)
  data.samples.l = data.samples.l[!names(data.samples.l) == ""]
  message("- level consistency\n")
  testLevelsConsistency(data.samples.l, type="Cal Std")
  ## Recovery sample
  data.samplesRecov.w=NULL
  Samples.whichRecovery = gsub("_.+$","", metaData2$Raw.File.Name[which(metaData2$Level == "QC_Rec")])
  if(length(Samples.whichRecovery)>0){
    data.samplesRecov.l = data.samples.l[names(data.samples.l) %in% Samples.whichRecovery]
    data.samplesRecov.w.l = lapply(data.samplesRecov.l,function(x) samples.analyze(x, is.recovery=TRUE))
    data.samplesRecov.w = do.call(rbind.fill,data.samplesRecov.w.l)
    data.samplesRecov.w$Batch = i
  }
  ## Normal Sample
  #### Check if the sample has only 1 injection level, which should be a problem.
  samples.singleInjection = names(which(lapply(data.samples.l, function(x) length(unique(x$Sample.Raw.File.Name)))==1))
  if(length(samples.singleInjection)>0) {
    stop("\n\n\nThe following sample(s) contain only one injection level:\n\t",paste0(samples.singleInjection, collapse=", "),"\n\n\n")
  }
  data.samples.l = data.samples.l[!names(data.samples.l) %in% Samples.whichRecovery]
  data.samples.w.l = lapply(data.samples.l,samples.analyze)
  Samples.w = do.call(rbind.fill,data.samples.w.l)
  Samples.w$Batch = i
  Samples.w = cbind(ID=apply(Samples.w[,c("Batch","Sample.name","Compound.Name")],1,function(row) paste0(row,collapse="_")), Samples.w)
  Samples.w$ID = gsub("-|\\.| ","_",Samples.w$ID)
  #################################
  #######   Simple Report   #######
  #################################
  message("-- Building simple report tale\n")
  Simple.table = Samples.w[,c("ID","Batch","Sample.name","Compound.Name","CAS.Number","Family","LOQ","QN/QL","Concentration","ccCol")]
  if(!is.null(Samples.w$Screening)) Simple.table$Screening = Samples.w$Screening
  if(!is.null(data.samplesRecov.w$QC.Screening)) Simple.table = cbind(Simple.table, QC.Screening = data.samplesRecov.w$QC.Screening[match(Simple.table$Compound.Name,data.samplesRecov.w$Compound.Name)])
  if(!is.null(data.samplesRecov.w$qcsCol)) Simple.table$qcsCol = data.samplesRecov.w$qcsCol
  if(!is.null(data.samplesRecov.w)) Simple.table = cbind(Simple.table, QC.recovery = data.samplesRecov.w$QC.recovery[match(Simple.table$Compound.Name,data.samplesRecov.w$Compound.Name)])
  ## Merge while keeping order (merge sort false does not work)
  Simple.table = cbind(Simple.table, QC.repeatability = QC.w$QC.repeatability[match(Simple.table$Compound.Name,QC.w$Compound.Name)])
  if("QC.linearity" %in% names(QC.w)) Simple.table = cbind(Simple.table, QC.linearity = QC.w$QC.linearity[match(Simple.table$Compound.Name,QC.w$Compound.Name)])
  Simple.table = cbind(Simple.table, CarryOver = Blank.w$CarryOver[match(Simple.table$Compound.Name,Blank.w$Compound.Name)])
  #
  Simple.table$Concentration = as.vector(Simple.table$Concentration)
  Concentration.num = suppressWarnings(as.numeric(Simple.table$Concentration))
  whichBLOQ = !is.na(Concentration.num) & Concentration.num < Simple.table$LOQ ## Supressing warnings because there are non numeric values that I want to ignore anyway
  whichBBLOQ = !is.na(Concentration.num) & Concentration.num < (Simple.table$LOQ*0.3)
  whichPos = !is.na(Concentration.num) & Concentration.num >= 0
  Simple.table$Concentration[which(whichBLOQ & whichPos)] <- paste0(ifelse(is.screening, "<STC", "<LOQ"), " (",Simple.table$Concentration[whichBLOQ & whichPos],")")
  Simple.table$Concentration[which(whichBBLOQ)] <- ifelse(is.screening, "<STC", "<LOQ")## also allow negatives.
  w.QL = Simple.table$"QN/QL" == "QL" | (is.na(Simple.table$"QN/QL") & is.screening)
  Simple.table$Concentration[which(w.QL & whichPos)] <- "Detected"## TODO, if QN/QL has no value, should be QL for screening. Else, block analysis if QN/QL is not fully filled.
  Simple.table$Concentration[which(w.QL & whichBLOQ)] <- "Detected"
  ### If screening, the decision is defined that the detection is <STC and is validated though it may not be a quantitative measure.
  if(is.screening) Simple.table$Concentration[which(Simple.table$Screening=="Negative" & Simple.table$Concentration=="Detected")] <- "<STC"
  ########
  metaData2 = cbind(Batch = i, metaData2)
  allData = cbind(Batch = i, allData)
  ##################################################
  #######  MERGING data for multiple inputs  #######
  ##################################################
  m.is.screening = all(c(m.is.screening, is.screening))
  m.metaData2 = dplyr::bind_rows(m.metaData2, metaData2)

  m.Blank.w = dplyr::bind_rows(m.Blank.w,Blank.w)
  m.QC.w = dplyr::bind_rows(m.QC.w,QC.w)
  m.Samples.w = dplyr::bind_rows(m.Samples.w,Samples.w)
  m.data.samplesRecov.w = dplyr::bind_rows(m.data.samplesRecov.w,data.samplesRecov.w)
  m.Simple.table = dplyr::bind_rows(m.Simple.table,Simple.table)

  m.allData = dplyr::bind_rows(m.allData, allData)
}

### Now removing the QC linearity and repeatability columns
m.QC.w$QC.linearity = NULL
m.QC.w$QC.repeatability = NULL
### Logically, also removing the Screening column as well
# Samples.w$Screening = NULL### Operator preference to keep it in
metaDataTable = data.frame(matrix(ncol=4))
for(i in 1:length(l.metaData1)){
  metaData1 = l.metaData1[[i]]
  metaData1.1 = metaData1[1:which(metaData1=="Operator"),]
  metaData1.2 = metaData1[(which(metaData1=="Operator")+1):nrow(metaData1),]
  metaData1 = data.frame(matrix(NA, nrow = max(nrow(metaData1.1),nrow(metaData1.2)), ncol = 4))
  metaData1[1:nrow(metaData1.1),1:2] = metaData1.1
  metaData1[1:nrow(metaData1.2),3:4] = metaData1.2
  metaDataTable = rbind(metaDataTable, c(names(l.metaData1)[i],"","",""))
  metaDataTable = rbind(metaDataTable, metaData1)
}
metaDataTable = metaDataTable[-1,]

is.overlapingBatches = length(unique(m.allData$Compound.Name)) != nrow(unique(m.allData[,c("Compound.Name","Batch")]))
metaData2 = m.metaData2
allData = m.allData
is.screening = m.is.screening
Blank.w = m.Blank.w
QC.w = m.QC.w
Samples.w = m.Samples.w
data.samplesRecov.w = m.data.samplesRecov.w
Simple.table = m.Simple.table
if(nrow(data.samplesRecov.w)==0) data.samplesRecov.w=NULL

#########################
#########################
##### Build Report ######
#########################
#########################
message("---- Building report ----\n")
wb <- createWorkbook()
options("openxlsx.minWidth" = 14)### works, 15 is too big
modifyBaseFont(wb, fontSize = 9, fontColour = "black", fontName = "Calibri")
#########################
### New Info tab
addWorksheet(wb, sheet = "Info")
nmd = nrow(metaDataTable)
addStyle(wb, sheet = "Info", createStyle(fontSize=16,textDecoration ="bold",halign="center",valign="center",fgFill="#EAEAEA",border="Bottom",borderStyle="medium",borderColour="#444444"), rows = 1, cols = 1:6)## header of metadata
writeData(wb, sheet = "Info", x = c("Analytical Report - UHPLC-HRMS Analysis"), startRow = 1, startCol=1, colNames = FALSE)
mergeCells(wb, sheet = "Info", cols=1:(ncol(metaDataTable)+2), rows=1)
setRowHeights(wb, "Info", rows = 1, heights = 35)
sr=2
for(i in 1:length(l.metaData1)){
  name = names(l.metaData1)[i]
  df = l.metaData1[[name]]
  df.1 = df[1:which(df=="Operator"),]
  df.2 = df[(which(df=="Operator")+1):nrow(df),]
  df = data.frame(matrix(NA, nrow = max(nrow(df.1),nrow(df.2)), ncol = 4))
  df[1:nrow(df.1),1:2] = df.1
  df[1:nrow(df.2),3:4] = df.2
  ### add batch and name rows
  sr=sr+1
  # batch.color[name]
  writeData(wb, sheet = "Info", x = paste0("Batch ",which(names(l.metaData1)==name)), startRow = sr, startCol=2, colNames = FALSE)
  addStyle(wb, sheet = "Info", createStyle(fontSize=12,textDecoration ="bold",valign="center",halign="center",fgFill="#EAEAEA",border="Top",borderStyle="medium",borderColour="#444444"), rows = sr, cols = 2:5)## header of metadata
  writeData(wb, sheet = "Info", x = c(name), startRow = sr+1, startCol=2, colNames = FALSE)
  addStyle(wb, sheet = "Info", createStyle(fontSize=12,textDecoration ="bold",valign="center",fgFill="#EAEAEA",border="Bottom",borderStyle="medium",borderColour="#444444"), rows = sr+1, cols = 2:5)## header of metadata
  mergeCells(wb, sheet = "Info", cols=2:(ncol(metaDataTable)+1), rows=sr)
  mergeCells(wb, sheet = "Info", cols=2:(ncol(metaDataTable)+1), rows=sr+1)
  ### add data of batch
  sr=sr+2
  writeData(wb, sheet = "Info", x = df, startRow = sr, startCol=2, colNames = FALSE)
  addStyle(wb, sheet = "Info", createStyle(fontSize=10,textDecoration ="bold",valign="center",fgFill="#EAEAEA"), rows = (sr:(sr+nrow(df)-1)), cols = 2)## header of metadata
  addStyle(wb, sheet = "Info", createStyle(fontSize=10,textDecoration ="bold",valign="center",fgFill="#EAEAEA"), rows = (sr:(sr+nrow(df)-1)), cols = 4)## header of metadata
  #### whole surrounding
  addStyle(wb, sheet = "Info", createStyle(border="top",borderStyle="medium",borderColour="#444444"), rows = sr-2, cols = 2:5, stack=TRUE)
  addStyle(wb, sheet = "Info", createStyle(border="bottom",borderStyle="medium",borderColour="#444444"), rows = sr+(nrow(df)-1), cols = 2:5, stack=TRUE)
  addStyle(wb, sheet = "Info", createStyle(border="left",borderStyle="medium",borderColour="#444444"), rows = ((sr-2):(sr+(nrow(df)-1))), cols = 2, stack=TRUE)
  addStyle(wb, sheet = "Info", createStyle(border="right",borderStyle="medium",borderColour="#444444"), rows = ((sr-2):(sr+(nrow(df)-1))), cols = 5, stack=TRUE)
  sr=sr+nrow(df)
}
setColWidths(wb, "Info", cols=2:5, widths = 40, ignoreMergedCells = TRUE)



## data tables
simpSR = buildSheetbyGroups(Simple.table, wb, sht="Report")
sampSr = buildSheetbyGroups(Samples.w, wb, sht="Samples")
# saveWorkbook(wb,outFile, overwrite=TRUE)
if(!is.null(data.samplesRecov.w)) recovSr = buildSheetbyGroups(data.samplesRecov.w, wb, sht="QC_sample")
qcSr = buildSheetbyGroups(QC.w, wb, sht="QC_solvent") ### Sr not needed, catching the return to avoid a print
blkSr = buildSheetbyGroups(Blank.w, wb, sht="Blank") ### Sr not needed, catching the return to avoid a print
########  Add hyperlinks for easy access to information
STn =  names(Simple.table)[!names(Simple.table) %in% c("ccCol","qcsCol")][unlist(groupsColOrder(Simple.table[,!names(Simple.table) %in% c("ccCol","qcsCol")]))]### reordering as done in the buildSheetbyGroups
### Sample
Simple.table = checkAndModDataCol.injAreaRatios(Simple.table, "Concentration", Samples.w, is.screening)
if(is.overlapingBatches) {### new problematic function
  samples.cell.rowNum = paste0('CELL("row",INDEX(Samples!$A$',sampSr+1,':$A$',sampSr+nrow(Samples.w),',MATCH(1,(Report!$B',(1:nrow(Simple.table))+simpSR,'=Samples!$B$',sampSr+1,':$B$',sampSr+nrow(Samples.w),') *(Report!$C',(1:nrow(Simple.table))+simpSR,'=Samples!$C$',sampSr+1,':$C$',sampSr+nrow(Samples.w),') *(Report!$D',(1:nrow(Simple.table))+simpSR,'=Samples!$D$',sampSr+1,':$D$',sampSr+nrow(Samples.w),'),0)))')
} else {### old functionning well method
  samples.cell.rowNum = paste0('CELL("row",INDEX(Samples!$A$',sampSr+1,':$A$',sampSr+nrow(Samples.w),',MATCH($A',(1:nrow(Simple.table))+simpSR,',Samples!$A$',sampSr+1,':$A$',sampSr+nrow(Samples.w),',0)))')
}

linkReport2sample = paste0('HYPERLINK("#Samples!"&',samples.cell.rowNum,'&":"&',samples.cell.rowNum,',"',Simple.table$Concentration,'")')

class(linkReport2sample) = c(class(linkReport2sample),"formula")
writeData(wb, "Report", x = linkReport2sample, startCol = which(STn=="Concentration"), startRow = simpSR+1, colNames=FALSE)
### Sample QC recovery
if("QC.recovery" %in% STn && !is.null(data.samplesRecov.w) ) {
  ## making sure we are matching the rightfull data using compound names
  m=match(Simple.table$Compound.Name,data.samplesRecov.w$Compound.Name)
  ## If theoretical value is 0, means no recovery and no test needed.
  Simple.table$QC.recovery[data.samplesRecov.w$recov.Level.theoretical[m]==0] <- "Not tested"
  ## If compound was not found for any injection, the compound is not analyzable
  confCols = grep("^Confirm_Inj_[1-9][0-9]?$",colnames(data.samplesRecov.w))
  for(confCol in confCols){
    Simple.table$QC.recovery[data.samplesRecov.w[m,confCol]=="NotFound"] <- "Not analyzable"
  }
  ## quotes non numeric cells
  Simple.table$QC.recovery = quoteNonNumeric(Simple.table$QC.recovery) ### make links
  recovNr = nrow(data.samplesRecov.w)
  if(is.overlapingBatches) {### new problematic function
    qcsamples.cell.rowNum = paste0('MATCH(1,(Report!$D$',(1:nrow(Simple.table))+recovSr,'=QC_sample!$B$',recovSr+1,':$B$',recovSr+nrow(data.samplesRecov.w),')* (Report!$B$',(1:nrow(Simple.table))+recovSr,'=QC_sample!$A$',recovSr+1,':$A$',recovSr+nrow(data.samplesRecov.w),'),0)')
  } else {### old functionning well method
    qcsamples.cell.rowNum = paste0('CELL("row",INDEX(QC_sample!$C$',recovSr+1,':$C$',recovSr+recovNr,',MATCH($D',(1:nrow(Simple.table))+simpSR,',QC_sample!$C$',recovSr+1,':$C$',recovSr+recovNr,',0)))')
  }
  linkReport2qcsample = paste0('HYPERLINK("#QC_sample!"&',qcsamples.cell.rowNum,'&":"&',qcsamples.cell.rowNum,',',Simple.table$QC.recovery,')')## this one is not quoted, need to quote characters
  class(linkReport2qcsample) = c(class(linkReport2qcsample),"formula")
  writeData(wb, "Report", x = linkReport2qcsample, startCol = which(STn=="QC.recovery"), startRow = simpSR+1, colNames=FALSE)##
}
### Screening
if("Screening" %in% STn) {
  ### samples.cell.rowNum is technically the same as for the Concentration link, but because it may be removed later, I prefer to keep them independent
  if(is.overlapingBatches) {### new problematic function
    samples.cell.rowNum = paste0('MATCH(1,(Report!$D$',(1:nrow(Simple.table))+sampSr,'=Samples!$D$',sampSr+1,':$D$',sampSr+nrow(Samples.w),')* (Report!$C$',(1:nrow(Simple.table))+sampSr,'=Samples!$C$',sampSr+1,':$C$',sampSr+nrow(Samples.w),')* (Report!$B$',(1:nrow(Simple.table))+sampSr,'=Samples!$B$',sampSr+1,':$B$',sampSr+nrow(Samples.w),'),0)')
  } else{### old functionning well method
    samples.cell.rowNum = paste0('CELL("row",INDEX(Samples!$A$',sampSr+1,':$A$',sampSr+nrow(Samples.w),',MATCH($A',(1:nrow(Simple.table))+simpSR,',Samples!$A$',sampSr+1,':$A$',sampSr+nrow(Samples.w),',0)))')
  }
  linkReport2sample = paste0('HYPERLINK("#Samples!"&',samples.cell.rowNum,'&":"&',samples.cell.rowNum,',"',Simple.table$Screening,'")')
  class(linkReport2sample) = c(class(linkReport2sample),"formula")
  writeData(wb, "Report", x = linkReport2sample, startCol = which(STn=="Screening"), startRow = simpSR+1, colNames=FALSE)##
}
if("QC.Screening" %in% STn) {
  recovNr = nrow(data.samplesRecov.w)
  if(is.overlapingBatches) {### new problematic function
    qcsamples.cell.rowNum = paste0('MATCH(1,(Report!$D$',(1:nrow(Simple.table))+recovSr,'=QC_sample!$B$',recovSr+1,':$B$',recovSr+nrow(data.samplesRecov.w),')* (Report!$B$',(1:nrow(Simple.table))+recovSr,'=QC_sample!$A$',recovSr+1,':$A$',recovSr+nrow(data.samplesRecov.w),'),0)')
  } else {### old functionning well method
    qcsamples.cell.rowNum = paste0('CELL("row",INDEX(QC_sample!$C$',recovSr+1,':$C$',recovSr+recovNr,',MATCH($D',(1:nrow(Simple.table))+simpSR,',QC_sample!$C$',recovSr+1,':$C$',recovSr+recovNr,',0)))')
  }
  linkReport2qcsample = paste0('HYPERLINK("#QC_sample!"&',qcsamples.cell.rowNum,'&":"&',qcsamples.cell.rowNum,',"',Simple.table$QC.Screening,'")')
  class(linkReport2qcsample) = c(class(linkReport2qcsample),"formula")
  writeData(wb, "Report", x = linkReport2qcsample, startCol = which(STn=="QC.Screening"), startRow = simpSR+1, colNames=FALSE)##
}
### QC
if(is.overlapingBatches) {### new problematic function
  qc.cell.rowNum = paste0('MATCH(1,(Report!$D$',(1:nrow(Simple.table))+qcSr,'=QC_solvent!$B$',qcSr+1,':$B$',qcSr+nrow(QC.w),')* (Report!$B$',(1:nrow(Simple.table))+qcSr,'=QC_solvent!$A$',qcSr+1,':$A$',qcSr+nrow(QC.w),'),0)')
} else {### old functionning well method
  qc.cell.rowNum = paste0('CELL("row",INDEX(QC_solvent!$B$',qcSr+1,':$B$',qcSr+nrow(QC.w),',MATCH($D',(1:nrow(Simple.table))+simpSR,',QC_solvent!$B$',qcSr+1,':$B$',qcSr+nrow(QC.w),',0)))')
}
linkReport2qc = paste0('HYPERLINK("#QC_solvent!"&',qc.cell.rowNum,'&":"&',qc.cell.rowNum,',"',Simple.table$QC.repeatability,'")')
class(linkReport2qc) = c(class(linkReport2qc),"formula")
writeData(wb, "Report", x = linkReport2qc, startCol = which(STn=="QC.repeatability"), startRow = simpSR+1, colNames=FALSE)

if("QC.linearity" %in% STn) {
  if(is.overlapingBatches) {### new problematic function
    #
  } else {### old functionning well method
    qc.cell.rowNum = paste0('CELL("row",INDEX(QC_solvent!$B$',qcSr+1,':$B$',qcSr+nrow(QC.w),',MATCH($D$',(1:nrow(Simple.table))+simpSR,',QC_solvent!$B$',qcSr+1,':$B$',qcSr+nrow(QC.w),',0)))')
  }
  linkReport2qc = paste0('HYPERLINK("#QC_solvent!"&',qc.cell.rowNum,'&":"&',qc.cell.rowNum,',"',Simple.table$QC.linearity,'")')
  class(linkReport2qc) = c(class(linkReport2qc),"formula")
  writeData(wb, "Report", x = linkReport2qc, startCol = which(STn=="QC.linearity"), startRow = simpSR+1, colNames=FALSE)
}
### Blank
if(is.overlapingBatches) {### new problematic function
  warning("Batches mix have overlapping names of compounds, using adapted matching function that is bugged by EXCEL.")
  blk.cell.rowNum = paste0('MATCH(1,(Blank!$A$',blkSr+1,':$A$',blkSr+nrow(Blank.w),'=Report!$B',(1:nrow(Simple.table))+simpSR,')*(Blank!$B$',blkSr+1,':$B$',blkSr+nrow(Blank.w),'=Report!$D',(1:nrow(Simple.table))+simpSR,'),0)')
} else {### old functionning well method
  blk.cell.rowNum = paste0('CELL("row",INDEX(Blank!$B$',blkSr+1,':$B$',blkSr+nrow(Blank.w),',MATCH(Report!$D',(1:nrow(Simple.table))+simpSR,',Blank!$B$',blkSr+1,':$B$',blkSr+nrow(Blank.w),',0)))')### working old way
}
#### This one should work, but excel adds a damn @ before the blank in the match formula. replace all @ with nothing works. File is correctly generated, office bug: https://answers.microsoft.com/en-us/msoffice/forum/all/excel-365-functions-preceded-by-in-formulas/4c0cc0f1-ed2b-4d37-83ff-a1aca3b3fac2
linkReport2blk = paste0('HYPERLINK("#Blank!"&',blk.cell.rowNum,'&":"&',blk.cell.rowNum,',"',Simple.table$CarryOver,'")')


class(linkReport2blk) = c(class(linkReport2blk),"formula")
writeData(wb, "Report", x = linkReport2blk, startCol = which(STn=="CarryOver"), startRow = simpSR+1, colNames=FALSE)
########
### Adding SAP column
sapColN = length(STn)+1
sapCol = Simple.table$Concentration
w.sapCol = which(grepl("^<STC|^<LOQ|^Not detected$",sapCol) & Simple.table$"QN/QL"=="QN")
sapCol[w.sapCol] <- paste0("<",Simple.table[w.sapCol,grep("^STC|^LOQ",colnames(Simple.table))])
writeData(wb, sheet = "Report", x = data.frame("SAP"=sapCol), headerStyle = headerStyleSimple, startCol = sapColN, startRow = simpSR)
addStyle(wb, sheet = "Report", style=createStyle(halign="center", valign="center"), rows = simpSR, cols = sapColN)
addStyle(wb, sheet = "Report", style=resArea, rows = (1:nrow(Simple.table)+simpSR), cols = sapColN)
setColWidths(wb, "Report", cols=sapColN, widths = 20,hidden=TRUE, ignoreMergedCells = TRUE) ## Auto seems to be problematic here. probably because of formula

########
### Adding & linking comment
Samples.w.c = Samples.w[,!names(Samples.w) %in% "ccCol"]
commentColN = letters[which(readWorkbook(wb, "Samples", rows=2, colNames=FALSE)=="Comment")]
commentLink = paste0('INDEX(Samples!$',commentColN,'$',sampSr+1,':$',commentColN,'$',sampSr+nrow(Samples.w),',MATCH($A',(1:nrow(Simple.table))+simpSR,',Samples!$A$',sampSr+1,':$A$',sampSr+nrow(Samples.w),',0))')
class(commentLink) = c(class(commentLink),"formula")
cmtColN = length(STn)+2
writeData(wb, sheet = "Report", x = data.frame("Comment"=commentLink), headerStyle = headerStyleSimple, startCol = cmtColN, startRow = simpSR)
addStyle(wb, sheet = "Report", style=createStyle(halign="center", valign="center"), rows = simpSR, cols = cmtColN)
addStyle(wb, sheet = "Report", style=commentStyle, rows = (1:nrow(Simple.table)+simpSR), cols = cmtColN)
setColWidths(wb, "Report", cols=cmtColN, widths = 50, ignoreMergedCells = TRUE) ## Auto seems to be problematic here. probably because of formula
#### Modify header of report
addStyle(wb, sheet = "Report", createStyle(fontSize=16,textDecoration ="bold",halign="center",valign="center",fgFill="#EAEAEA",border="Bottom",borderStyle="medium",borderColour="#444444"), rows = 1, cols = 1:cmtColN)## header of metadata
writeData(wb, sheet = "Report", x = c("Simplified Results Report"), startRow = 1, startCol=1, colNames = FALSE)
mergeCells(wb, sheet = "Report", cols=1:cmtColN, rows=1)
setRowHeights(wb, "Report", rows = 1, heights = 35)
## remove formating from the 1st of the 2 rows header of results table
removeCellMerge(wb, sheet = "Report", cols=1:(ncol(Simple.table)), rows=2)
addStyle(wb, sheet = "Report", style=createStyle(), cols=1:cmtColN, rows=2, stack=FALSE)
writeData(wb, sheet = "Report",x=t(rep("",20)), startCol=1, startRow=2, colNames=FALSE)
## Set wide line at top of the table
addStyle(wb, sheet = "Report", style=headerStyleSimple, cols=1:cmtColN, rows=simpSR, stack=TRUE)


############
### AllSample
addWorksheet(wb, sheet = "allData")
writeData(wb, sheet = "allData", x = allData, headerStyle = headerStyle)
setRowHeights(wb, "allData", rows = 1, heights = 25)
setColWidths(wb, "allData", cols=1:ncol(allData), widths = "auto")
addFilter(wb, "allData", row = 1, cols = 1:ncol(allData))
freezePane(wb, "allData", firstRow = TRUE)
############
### Batch
addWorksheet(wb, sheet = "Batch")
writeData(wb, sheet = "Batch", x = metaData2, headerStyle = headerStyle)
setRowHeights(wb, "Batch", rows = 1, heights = 25)
setColWidths(wb, "Batch", cols=1:ncol(metaData2), widths = "auto")
addFilter(wb, "Batch", row = 1, cols = 1:ncol(metaData2))
freezePane(wb, "Batch", firstRow = TRUE)
############
#### Protect certain worksheet for accidental editing
protectWorksheet(wb,sheet = "allData", protect = TRUE, lockAutoFilter = FALSE, lockSelectingLockedCells = FALSE, lockSorting = FALSE, lockFormattingColumns = FALSE)
protectWorksheet(wb,sheet = "Batch", protect = TRUE, lockAutoFilter = FALSE, lockSelectingLockedCells = FALSE, lockSorting = FALSE, lockFormattingColumns = FALSE)
########################
###  Write WorkBook  ###
message("---- Writing to file ----\n")
# saveWorkbook(wb,outFile, overwrite=TRUE)
tryCatch(
  expr = {saveWorkbook(wb,outFile, overwrite=TRUE)},
  error = function(e){stop(e)},
  warning = function(w){warning(w)}
)
########################





#
