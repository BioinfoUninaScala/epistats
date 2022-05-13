#' Extracting epiallele information from customised target regions.
#' @description
#' This is the function which takes as input the user-defined regions of interest and it returns for each of them epiallele composition and summary statistics, such as average DNA methylation, Shannon Entropy, etc...
#' @importFrom BSgenome getSeq
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom purrr map_df
#' @importFrom Biostrings vmatchPattern matchPattern reverseComplement DNAString strsplit
#' @importFrom XVector subseq
#' @importFrom GenomicAlignments stackStringsFromGAlignments
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom magrittr %>%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom tidyr as_tibble
#' @param align A GAlignments object obtained trough the loading and filtering steps.
#' @param bin A GRanges object that is the output obtained from one of the functions makeWindows() and makeBins().
#' @param aname A character indicating the analysis name ("SampleName")
#' @param threshold An integer indicating the minimum coverage required to analyse the interval. This is specified here since it could be possible that DNA methylation status cannot be assessed for some reads, thus decreasing altering the coverage of the previously selected regions.
#' @param pathDir A string indicating the path of the directory desired to save the output.
#' @param genome A DNAStringSet object containing the fasta genome
#' @param myfuns A list object containing the functions to be applied to the binary epiallele matrix (to compute average DNA methylation, Shannon Entropy, etc...).
#' @param bisu.Thresh An integer indicating the bisulfite efficiency threshold to be used to discard reads or not (Default = 0).
#' @param stranded Logical indicating whether the analysis should be performed stranded or not (Default = FALSE).
#' @param mode A character indicating the pattern to be used to perform the methylation analysis (one of "CG", "CA", "CC", "CT").
#' @param remove.Amb Logical indicating whether ambiguous reads (containing not expected nucleotides) should be discarded or not (Default = TRUE).
#' @param cores An integer indicating the number of cores to be used for the computation.
#' @param retain.reads Logical indicating whether reads that have a low bisulfite conversion efficiency should be kept or not (Default = TRUE).
#' @param get.cPos Logical indicating whether the coordinates of the CpGs displayed through the analysis should be returned as a text output or not (Default = FALSE).
#' @return A list object containing different dataframes.
#'
#' intervals = Dataframe having all the analysed regions as rows and the different summary statistics as columns.
#'
#' epi = Dataframe that is the compressed form of the binary matrix containing the epiallele composition for each analysed region.
#'
#' log = Dataframe containing the coordinates of the genomic regions that have been discarded for the analysis.
#' @export
#' @examples
#' data <- loadInput(bamfile, genomefile)
#' algn <- data[[1]]
#' Genome <- data [[2]]
#'
#' ## Keeping only standard chromosomes
#' filtered <- filterBAM(algn,
#'                       keepStChr =TRUE,
#'                       keepM = FALSE,
#'                       retainChr = NULL)
#'
#' ## Selecting regions with a minimum coverage
#' covered <- filterByCoverage(algn, threshold = 50, minsize = 8)
#'
#' ## Choose one of the two available methods to design target regions.
#' windows <- makeWindows(gr = covered,
#'                        Genome = Genome,
#'                        window = 50,
#'                        step = 1,
#'                        mode = "CG",
#'                        min.C = 2,
#'                        max.C = 50,
#'                        cores = 40)
#'
#' ## Extract the GRanges object containing the target regions
#' regions <- windows$windows
#'
#' ## Select the functions that work on the binary epiallele matrix
#' myfuns <- list("dist" = cdist,
#'                "epi" = epi,
#'                "singleton" = singleton,
#'                "maxfreq" = maxfreq,
#'                "shannon" = shannon,
#'                "mean_met" = meanMeth,
#'                "num_cg" = numCG,
#'                "num_reads" = num_reads)
#'
#' ## Extract epiallele information
#' aname = "Sample_1"
#'
#' out <- epiAnalysis(align = algn,
#'                    bin = regions,
#'                    aname = aname,
#'                    threshold = 50,
#'                    bisu.Thresh = 0,
#'                    stranded = FALSE,
#'                    mode = "CG",
#'                    remove.Amb = TRUE,
#'                    genome = Genome,
#'                    pathDir = "/home/sarnataro/bam/out_epistats/",
#'                    cores = cores,
#'                    retain.reads = TRUE,
#'                    get.cPos = FALSE,
#'                    myfuns = myfuns)

epiAnalysis= function(align,
                      bin,
                      aname,
                      threshold,
                      pathDir,
                      genome,
                      myfuns,
                      bisu.Thresh = 0,
                      stranded=FALSE,
                      mode="CG",
                      remove.Amb =TRUE,
                      cores = 1,
                      retain.reads = FALSE,
                      get.cPos=FALSE)
{
  ###controllo sui parametri
  if((!mode=="CG") & (stranded == FALSE))
  {
    print(paste("unstranded mode not supported for",mode, "mode",sep=" "))
  } else {
    ## Takes all the refseq in my intervals coordinates
    refseq=BSgenome::getSeq(genome, bin)
    ## Do epiallele_analyse for all the intervals
    df= cbind(data.frame(bin)[,1:3], refseq=as.character(refseq))

    df_split <- split(df, (seq(nrow(df))-1) %/% (nrow(df)/cores))

    regs <- foreach::foreach (dfi = df_split) %do% {
      region = GenomicRanges::GRanges(dfi$seqnames, IRanges::IRanges(dfi$start,dfi$end))
      IRanges::subsetByOverlaps(align, region)
    }
    ## Plan parallel
    cl <- parallel::makeCluster(cores, type = 'PSOCK')
    doParallel::registerDoParallel(cl)
    Mapi = foreach::foreach(a = df_split, b = regs, .verbose = TRUE) %dopar% {
      block = epiallele_analyse_Block(a, b, threshold, bisu.Thresh, stranded, mode, remove.Amb, retain.reads, get.cPos, myfuns)
    }
    parallel::stopCluster(cl)
    ## Unlist dataframes
    Mapi <- unlist(Mapi, recursive = FALSE)
    ##
    intervals <- Mapi %>% purrr::map_df(~ .$intervals)
    epi <- Mapi %>% purrr::map_df(~ .$epi)
    log <- Mapi %>% purrr::map_df(~ .$log)
    #### Costruisce i dataframe e li scrive nei rispettivi file
    #out=do.call(Mapi, c(rbind, out))
    #setta i file e i nomi di colonna
    bedFile=paste(pathDir,paste(aname,"intervals.bed", sep="_"),sep="/")
    write.table(intervals,bedFile, row.names = F, col.names = T, quote= F, sep = "\t")
    epiFile=paste(pathDir,paste(aname,"epiAnalysis.txt", sep="_"),sep="/")
    write.table(epi, epiFile, row.names = F, col.names = T, quote= F, sep ="\t")
    logFile=paste(pathDir,paste(aname,"log.txt", sep="_"),sep="/")
    write.table(log, logFile, row.names = F, col.names = F, quote= F, sep ="\t")
    if(get.cPos==TRUE)
    {
      cPosFile=paste(pathDir,paste(aname,"cPos.txt", sep="_"), sep="/")
      write.table(out[["CPos"]], cPosFile, row.names = F, col.names = F, quote = F, sep="\t")
    }
  }
  return(list(intervals,epi,log))
}


read_filter=function(alignObj)
{
  align_trunk= alignObj[Biostrings::vmatchPattern("+", alignObj)]
  alignObj= alignObj[which(align_trunk@ranges@width == 0)]
  alignObj= alignObj[!duplicated(names(alignObj))]
  return(alignObj)
}


get_cPos=function(rseq, mode, strand, bisu.Thresh)
{
  if (strand == "plus")
  {
    cMode = Biostrings::matchPattern(mode, rseq)@ranges@start
    ###solo strand plus
    if (bisu.Thresh==0)
    {
      ###se non deve fare il controllo del bisulfito
      cPos=list(as.list(cMode),list(),list(),list())
    }else{
      ###se deve fare il controllo del bisulfito
      if (mode == "CG")
      {
        #sulle cpg
        cNotMode <- Biostrings::matchPattern("C", rseq)@ranges@start
        cNotMode= cNotMode[!cNotMode %in% cMode]
      } else {
        #sulle non cpg
        cNotMode <- Biostrings::matchPattern("C", rseq)@ranges@start
        cg_pos = Biostrings::matchPattern("CG", rseq)@ranges@start
        cNotMode = cNotMode[!cNotMode %in% c(cMode,cg_pos)]
      }
      cPos=list(as.list(cMode), as.list(cNotMode),list(),list())
    }

  } else if (strand == "minus"){
    ###solo strand minus
    if(mode=="CG")
    {
      cMode=Biostrings::matchPattern(mode, rseq)@ranges@start
      cMode=cMode
    }else{
      cMode=Biostrings::matchPattern(Biostrings::reverseComplement(Biostrings::DNAString(mode)), rseq)@ranges@start
      cMode=cMode
    }
    if (bisu.Thresh==0)
    {
      ###se non deve fare il controllo del bisulfito
      cPos=list(list(),list(),as.list(cMode),list())
    }else{
      ###se deve fare il controllo del bisulfito
      if (mode == "CG")
      {
        #sulle cpg
        cNotMode <- Biostrings::matchPattern("G", rseq)@ranges@start
        cNotMode= cNotMode[!cNotMode %in% cMode]
      } else {
        #sulle non cpg
        cNotMode <- Biostrings::matchPattern("G", rseq)@ranges@start
        cg_pos = (Biostrings::matchPattern("CG", rseq)@ranges@start)+1
        cNotMode = cNotMode[!cNotMode %in% c(cMode,cg_pos)]
      }
      cPos=list(list(),list(),as.list(cMode), as.list(cNotMode))
    }

  }else{
    ###entrambi gli strand
    cMode_plus = Biostrings::matchPattern(mode, rseq)@ranges@start
    if(mode=="CG")
    {
      cMode_minus=cMode_plus
    }else{
      cMode_minus = Biostrings::matchPattern(Biostrings::reverseComplement(Biostrings::DNAString(mode)), rseq)@ranges@start
    }
    if (bisu.Thresh==0)
    {
      ###se non deve fare il controllo del bisulfito
      cPos=list(as.list(cMode_plus),list(),as.list(cMode_minus),list())
    }else{
      ###se deve fare il controllo del bisulfito
      if (mode == "CG")
      {
        #sulle cpg
        cNotMode_plus=Biostrings::matchPattern("C",rseq)@ranges@start
        cNotMode_plus= cNotMode_plus[!cNotMode_plus %in% cMode_plus]
        cNotMode_minus=Biostrings::matchPattern("G",rseq)@ranges@start
        cNotMode_minus=cNotMode_minus[!cNotMode_minus %in% as.integer(cMode_minus+1)]
      } else {
        #sulle non cpg
        cNotMode_plus <- Biostrings::matchPattern("C", rseq)@ranges@start
        cg_pos = Biostrings::matchPattern("CG", rseq)@ranges@start
        cNotMode_plus = cNotMode_plus[!cNotMode_plus %in% c(cMode_plus,cg_pos)]
        cNotMode_minus <- Biostrings::matchPattern("G", rseq)@ranges@start
        cg_pos = (Biostrings::matchPattern("CG", rseq)@ranges@start)+1
        cNotMode_minus = cNotMode_minus[!cNotMode_minus %in% c(cMode_minus+1,cg_pos)]
      }
      cPos=list(as.list(cMode_plus),as.list(cNotMode_plus),as.list(cMode_minus),as.list(cNotMode_minus))
    }
  }
  return(cPos)
}


extract_matrix=function(alignObj, pos, strand, mode)
{
  if(mode=="bisu"){
    data=lapply(pos, function(x) {y = XVector::subseq(alignObj, x, x); return(unname(as.vector(y)))})
    data=as.data.frame(do.call(cbind,data))
    if(strand == "plus")
    {
      data <- as.data.frame(lapply(data, function(x) ifelse(x == "C", 1, ifelse(x == "T", 0, NA))))
    }else{
      data <- as.data.frame(lapply(data, function(x) ifelse(x == "G", 1, ifelse(x == "A", 0, NA))))
    }
  } else{
    char= unlist(Biostrings::strsplit(mode, split=""))
    data=lapply(pos, function(x) {y = XVector::subseq(alignObj, x, x+1); return(unname(as.vector(y)))})
    data=as.data.frame(do.call(cbind,data))
    if(strand == "plus")
    {
      char= paste(c("C", "T"), char[2], sep="")
      data <- as.data.frame(lapply(data, function(x) ifelse(x == char[1], 1, ifelse(x == char[2], 0, NA))))
    }else{ ### reverse complement
      char= paste(Biostrings::reverseComplement(Biostrings::DNAString(char[2])), c("G", "A"), sep="")
      data <- as.data.frame(lapply(data, function(x) ifelse(x == char[1], 1, ifelse(x == char[2], 0, NA))))
    }
  }
  # rownames(data)=names(alignObj)
  return(data)
}


get_epiMatrix=function(alignObj, bisu.Thresh, remove.Amb, c.Mode, c.NotMode, strand, retain.reads, mode)
{
  data_mode=extract_matrix(alignObj, c.Mode, strand, mode)
  if(remove.Amb==TRUE)
  {
    index=which(rowSums(is.na(data_mode))>0)
    if(length(index)>0)
    {
      data_mode=data_mode[-index,]
    }
  }
  ######passa al controllo del bisulfito
  ######se il cutoff e' zero oppure se non e' zero ma non ci sono le c su cui effettuarlo e l'utente ha specificato di consevare le reads
  if (bisu.Thresh==0 | (bisu.Thresh >0 & (length(c.NotMode)==0 & retain.reads==TRUE)))
  {
    ####se la soglia e' zero
    return(data_mode)
  }else{
    ####se la soglia non e' zero
    if(length(c.NotMode) == 0 & retain.reads== FALSE)
      ######se non ci sono c su cui effettuare il controllo e l'utente ha specificato di non conservare le reads
    {
      return(data.frame())
    }else{
      #####se ci sono le c su cui effettuare il controllo
      data_bisu=extract_matrix(alignObj, c.NotMode, strand, mode= "bisu")
      ###eventualmente va a rimuovere le reads eliminate con amb thresh.
      if(length(index)>0)
      {
        data_bisu=data_bisu[-index,]
      }
      ###esegue filtro del bisulfito
      eff=apply(as.matrix(data_bisu), 1, function(x) 1-(sum(x, na.rm = TRUE)/length(x[!is.na(x)])))
      ########azzera l'efficienza delle reads in cui le non cg sono tutte na
      eff[is.na(eff)]=0
      data_mode=data_mode[eff >=bisu.Thresh,]
      return(data_mode)
    }
  }
}


get_out=function(matrix, out, strand, coord, get.cPos, myfuns)
{
  parameters=lapply(myfuns, function(x) x(matrix))
  mydata= data.frame(as.data.frame(coord),
                     as.data.frame(parameters))
  mydata$strand=strand
  out[["intervals"]]=mydata
  id=apply(mydata[1,c(1:3)],1,function(x) paste(x,collapse="_"))
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix=as.data.frame(table(matrix$epi))
  matrix$id=id
  matrix$strand=strand
  out[["epi"]]=matrix
  if(get.cPos == TRUE)
  {
    Cpos=as.numeric(names(matrix))
    out[["CPos"]]=data.frame('chr'= rep(as.character(seqnames(coord)),length(Cpos)),"start"= Cpos + coord@ranges@start - 1,"end"= Cpos + coord@ranges@start, "strand"=rep(strand, length(Cpos)), "id"=rep(id,length(Cpos)))
  }
  return(out)
}

epiallele_analyse=function(align,
                           bin,
                           threshold,
                           bisu.Thresh,
                           stranded,
                           mode,
                           remove.Amb,
                           rseq,
                           retain.reads,
                           get.cPos,
                           myfuns){
  if(get.cPos==FALSE)
  {
    out=list("intervals"=data.frame(),"epi"=data.frame(),"log"=data.frame())
  }else{
    out=list("intervals"=data.frame(),"epi"=data.frame(),"log"=data.frame(),"CPos"=data.frame())
  }
  ##############################restituisce gli allineamenti delle reads per la regione di interess
  reads <- GenomicAlignments::stackStringsFromGAlignments(align, bin)
  reads <- reads[!duplicated(names(reads))]
  reads_trunk <- reads[Biostrings::vmatchPattern("+", reads)]
  reads <- reads[which(reads_trunk@ranges@width == 0)]
  ## Plus
  align_plus= reads[reads@elementMetadata$strand=="+"]
  ## Minus
  align_minus= reads[reads@elementMetadata$strand=="-"]
  ############################
  ####fa il check dello strand e degli allineamenti, per capire su quali richiamare getMatrix
  if (bin@strand@values =="*" & (length(align_minus)==0 & length(align_plus)==0) |
      bin@strand@values =="-" & length(align_minus)==0|
      bin@strand@values =="+" & length(align_plus)==0 |
      length(reads) < threshold)
  {
    out[["log"]]=as.data.frame(bin)
  } else {
    ##se strand * ed entrambi gli allineamenti sono non vuoti
    if(bin@strand@values == "*" & (length(align_minus)>0 & length(align_plus)>0))
    {
      cPos=get_cPos(rseq, mode,"*",bisu.Thresh)
      data_plus=tidyr::as_tibble(get_epiMatrix(align_plus, bisu.Thresh, remove.Amb, cPos[[1]], cPos[[2]], strand = "plus", retain.reads = retain.reads, mode))
      data_minus=tidyr::as_tibble(get_epiMatrix(align_minus, bisu.Thresh, remove.Amb, cPos[[3]], cPos[[4]], strand = "minus", retain.reads = retain.reads, mode))
      ###se tutte le reads sono filtrate per ambthresh o bisulfito, scrive bin nel log
      if (dim(data_plus)[1]==0 & dim(data_minus)[1]==0 | dim(data_plus)[1]+dim(data_minus)[1] < threshold)
      {
        out[["log"]]=as.data.frame(bin)
      }else{
        if(stranded == FALSE)
        {
          if(dim(data_plus)[1]>0 & dim(data_minus)[1]>0)
          {
            names(data_minus)=names(data_plus)
          }
          data=rbind(data_plus,data_minus)
          colnames(data)=unlist(cPos[[1]])
          out=get_out(data, out, "*", bin, get.cPos, myfuns)
        }else{
          if(dim(data_plus)[1]>0)
          {
            colnames(data_plus)=unlist(cPos[[1]])
            out=get_out(data_plus, out, "+", bin, get.cPos, myfuns)
          }
          if(dim(data_minus)[1]>0)
          {
            colnames(data_minus)=unlist(cPos[[3]])
            out=get_out(data_minus, out, "-", bin, get.cPos, myfuns)
          }
        }
      }
      ##se lo strand * e solo il positivo e' pieno oppure se lo strand e' +
    } else if (bin@strand@values == "*" & (length(align_minus)==0 & length(align_plus)>0)|
             bin@strand@values == "+")
    {
      cPos=get_cPos(rseq,mode,"plus",bisu.Thresh)
      data_plus=tidyr::as_tibble(get_epiMatrix(align_plus, bisu.Thresh, remove.Amb, cPos[[1]], cPos[[2]], "plus", retain.reads, mode))
      if (dim(data_plus)[1] < threshold)
      {
        out[["log"]]=as.data.frame(bin)
      }else{
        colnames(data_plus)=unlist(cPos[[1]])
        out=get_out(data_plus, out, bin@strand@values, bin, get.cPos, myfuns)
      }
    }else {
      ##se lo strand e' * e solo in negativo e' pieno oppure lo strand e' -
      cPos=get_cPos(rseq,mode,"minus",bisu.Thresh)
      data_minus=tidyr::as_tibble(get_epiMatrix(align_minus, bisu.Thresh, remove.Amb, cPos[[3]], cPos[[4]], "minus", retain.reads,mode))
      if (dim(data_minus)[1] < threshold)
      {
        out[["log"]]=as.data.frame(bin)
      }else{
        colnames(data_minus)=unlist(cPos[[3]])
        out=get_out(data_minus, out, bin@strand@values, bin, get.cPos, myfuns)
      }
    }
  }
  return(out)
}


epiallele_analyse_Block <- function(df,
                                    aln_c,
                                    threshold,
                                    bisu.Thresh,
                                    stranded,
                                    mode,
                                    remove.Amb,
                                    retain.reads,
                                    get.cPos,
                                    myfuns){
  ret_list = list()
  for (i in 1:nrow(df)) {
    ret_list[[i]] = epiallele_analyse(align= aln_c,
                                      bin = GenomicRanges::GRanges(df[i,]$seqnames,IRanges::IRanges(df[i,]$start,df[i,]$end), strand = df[i,]$strand),
                                      threshold = threshold,
                                      bisu.Thresh = bisu.Thresh,
                                      stranded = stranded,
                                      mode = mode,
                                      remove.Amb = remove.Amb,
                                      rseq = df[i,]$refseq,
                                      retain.reads = retain.reads,
                                      get.cPos = get.cPos,
                                      myfuns = myfuns)
  }
  return(ret_list)
}










