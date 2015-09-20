################################################################################
################################################################################
#                           gbparseR
#
#                         by Nick Waters 20150920
#                         Version 0_1_0
################################################################################
################################################################################
# This file currently is more of a workflow-  the functions will extract 
#informaton from .gb genbank files, uing a lot of regular expressions and the 
#readlines function.  Everything is base R, so that should make your life easier 
#when trying to manage packages..
#
# get sequences from NCBI ising the dropdown menu for exporting a sequence 
# after ensureing that you are displaying the full sequence at the bottom of the page.


################################################################################
################################################################################
source.file = "CP000253_8325.gb"

################################################################################
####  extract metadata and meat ( the actual features).  spits out a list 
#     containing:
#         the meta info
#         The features
#         the nucleotide sequence

load_gb<-function(source.file){
  input <- readLines(source.file)
  if(length(grep("LOCUS", input))==2){
    print("sorry, single scffolds only for now. Try splitting your .gb. file into
    a file for each scaffold")  
    meta<-""
    meat<-""
    for (i in 1:length(grep("LOCUS", input))){
      j<-input[grep("LOCUS", input)[i]:grep("FEATURES", input)[i]-1]
      meta<-c(meta, j)
    }
    notmeta<-input[!input %in% meta]
    metall<-meta
    meta1<-meta[grep("LOCUS", meta)[1]:grep("//$", meta)]
    meta2<-meta[grep("LOCUS", meta)[2]:length( meta)]
    #  meta<-meta[!duplicated(meta)]
    ####  butcher meat
    meat1<-notmeta[grep("FEATURES", notmeta)[1]:grep("ORIGIN", notmeta)[1]-1]
    meat2<-notmeta[grep("FEATURES", notmeta)[2]:grep("ORIGIN", notmeta)[2]-1]
    ####  get seq   make this a function
    seq1<-notmeta[grep("ORIGIN", notmeta)[1]: grep("FEATURES", notmeta)[2]-1]
    seq2<-notmeta[grep("ORIGIN", notmeta)[2]: length(notmeta)]
    
    if (length(grep("ORIGIN", notmeta))>2){print("sorry mate, can only deal with the first two for now!")}
  } else{
    meta<-input[grep("LOCUS", input)[1]:grep("FEATURES", input)[1]-1]
    meat<-input[!input %in% meta]
    meat<-meat[grep("FEATURES", meat):length(meat)]
    seq<-meat[grep("ORIGIN", meat): length(meat)]
    return(list(meta, meat, seq))
  }
}
#^^^^^^^^^  
#returns<-load_gb("CP000253_8325.gb")
################################################################################
#####  clean up sequence
clean_sequence<-function(seq){#} (seq in grep("seq\\d", ls(), value=T)){#print(seq)}
  fullraw<-seq #get(seq)
  nospace<-gsub("\\s", "",fullraw[2:length(fullraw)])
  nonum<-paste(gsub("\\d", "",nospace), sep="", collapse="")
  nonum<-substr(nonum, 7, nchar(nonum)) #  get rid of word "ORIGIN
  nonum
}

#^^^^^^^^^ 
#returns[[3]]<-clean_sequence(seq = returns[[3]])


################################################################################
###  start making  dataframe containing the location coordinates
get_ranges<-function(x){
  if (grep("FEATURES", x)!=1){stop}   #  make sure its 
  #  x<-meat1  #debug with
  x<-c(x, " 27..27 ")
  # ewxtract ranges 
  featList<-
    gsub("\\)|\\,|\\(|-","", 
         gsub("(\\D*)(\\d*\\.{2}\\d*)(\\D*)", "\\2", 
              
              grep("(\\d*\\.{2}\\d*)",x, value=T)))
  #featList<-gsub("#clean ranges
  indexList<-grep("(\\d*\\.{2}\\d*)",x)
  preZ<-data.frame(index= indexList, loc=featList, stringsAsFactors = F)
  # this next bit gets rid of duplicates;  neat, huh?
  for (i in 2:length(preZ$loc)){
    if (preZ[i-1,"loc"]==preZ[i, "loc"]){preZ[i,"loc"]<-"NA"}
  }
  preZ<-preZ[preZ$loc != "NA",]
  z<-preZ   
  z$id<-1:nrow(z)
  z$loc_start<-as.numeric(gsub("(\\d*)(\\.\\.\\d*)", "\\1", z$loc))
  z$loc_end<-as.numeric(gsub("(\\d*\\.\\.)(\\d*)", "\\2", z$loc))
  z$next_loc<-unlist(c(lapply(z$id[1:nrow(z)-1], function(w){       #  get start of next locus
    as.numeric(z[z$id==w+1, "loc_start"])}), 0))
  z$loc<-NULL
  z
  
}
#^^^^^^^^^ 
#ranges<-get_ranges(returns[[2]])


################################################################################
#  this removes references to locations from the locations lis, and 
#  relevels the id's to get continuous numbering

clean_ranges_rna<-function(gl, data){
  for (i in gl[1,"id"]:gl[nrow(gl) - 1, "id"]){
    j<-i+1
    entry<-data[gl[gl$id==i,"index"]:gl[gl$id==j,"index"]]
    gl[gl$id==i, "keep"]<-ifelse(
      grepl("anticodon=", entry[1]), "NA", "keep"
    )
  }
  gl[gl$id==max(gl$id, na.rm = T), "keep"]<-"keep"    #  take care of the buffer row;  need this to be included
  gl<-gl[gl$keep != "NA",]
  gl[,"keep"]<-NULL
  gl[,"id"]<-1:nrow(gl)#relevel id to get consec numbers
  gl
}

#^^^^^^^^^ 
#ranges<-clean_ranges_rna(ranges, returns)

################################################################################

#  this extracts the features, augmenting the ranges dataframe with gene info, 
#  type, etc
extract_features<-function(data, ranges, locus_tag="locus_tag"){ #y=grep("\\sgene\\s", x)[1]
  
  #   check for buffer row, stop if  absent
  if(tail(ranges,1)$loc_start!=27 | tail(ranges,1)$loc_end != 27){
    stop("where's the buffer row? did you forget to get and clean the ranges?")
  }
  z<-ranges
  for (i in ranges[1,"id"]:(ranges[nrow(ranges),"id"]-1)){ #  print(i)}
    #  can this be made into an apply?
    j<-i+1
    hits<-grep("\\d+\\.\\.\\d+", 
               data[ranges[ranges$id==i,"index"]:ranges[ranges$id==j,"index"]])
    #  this defines the index of where to look for the entry header;  default is 1:3 
    #  for non-gene features;  for genes, it finds the loc pattern \\d..\\d in the 
    #  region and bound of the entry as the second, middle occurance.
    entry_end <-ifelse(length(hits)>=3,
                       hits[2], 3)
    entry<-data[ranges[ranges$id==i,"index"]:ranges[ranges$id==j,"index"]]
    entry_head<-entry[1:entry_end]
    #}#  quick check
    if(!any(grepl("source|assembly_gap|misc_feature|misc_binding|CDS|rRNA|tRNA", 
                  entry_head))){
      if(any(grepl("pseudo", entry))){
        print(paste("caution! id = ", i, "is possibly a pseudogene"))
      } else {
        print(paste(" uh oh...  we got a rogue entry:",i,". Skipping..."))
        next()
      }
      #
    }
    # go time;  stuff is arranged most frequent to least.  probs will need
    #to add lotsa "try" stuff to this
    if(any(grepl("CDS", entry_head))){
      z[z$id==i, "type"]<-"CDS"
      
      pre_locus<-gsub(paste("(.*\\s*\\",locus_tag,"=)(.*)(\")", sep=""),"\\2", 
                      grep(locus_tag,entry, value=T )[1])
      try(z[z$id==i, "gene"]<-
            gsub("(.*\\s*\\gene=\")(.*)(\")","\\2", 
                 grep("gene=",entry, value=T )[1]), silent=T)
      
      z[z$id==i,"locus_tag"]<-
        gsub("(.*?\")(.*)", "\\2",  pre_locus)
      try(z[z$id==i, "inference"]<-
            gsub("(.*\\s*\\inference=)(.*\")(.*)","\\3", 
                 grep("inference",entry, value=T )),silent=T)
      #    z[i, "sequence_ref"]<-
      #      gsub("(.*\\s*\\sequence:)(.*)","\\2", grep("sequence:",entry, value=T ))
      z[z$id==i, "direction"]<-ifelse(!any(grepl("complement",entry_head)), "leading", "compliment")
      try(z[z$id==i, "note"]<-
            gsub("  ","",paste0(entry[grep("note=", entry): 
                                        (grep("codon_start",entry)-1)], collapse="  ")),
          silent=T)
      try(z[z$id==i, "codon_start"]<-
            gsub("(.*\\s*\\codon_start=)(.*)","\\2", 
                 grep("codon_start",entry, value=T )), silent=T)
      try(z[z$id==i, "transl_table"]<-
            gsub("(.*\\s*/transl_table=)(.*)","\\2", 
                 grep("transl_table",entry, value=T )), silent=T)
      try(z[z$id==i, "product"]<-
            gsub("(.*\\s*/product=)(\")(.*)(\")","\\3", grep("product",entry, value=T )),
          silent=T)
      try(z[z$id==i, "protein_id"]<-
            gsub("(.*\\s*/protein_id=)(\")(.*)(\")","\\3", 
                 grep("protein_id",entry, value=T )), silent=T)
      try(z[z$id==i, "db_xref"]<-
            gsub("(.*\\s*/db_xref=)(\")(.*)(\")","\\3", 
                 grep("db_xref",entry, value=T )), silent=T)
      try(z[z$id==i,"translation"]<-
            paste(gsub("([^A-Z]*)([A-Z]*)([^A-Z]*)","\\2",
                       gsub("translation=", "", entry[grep("translation=", entry):length(entry)])), 
                  collapse="",sep=""), silent=T)
      #####  
    } else if(any(grepl("pseudo", entry[length(entry)-1]))){
      z[z$id==i, "type"]<-"pseudo"
      pre_locus<-gsub(paste("(.*\\s*\\",locus_tag,"=)(.*)(\")", sep=""),"\\2", 
                      grep(locus_tag,entry, value=T )[1])
      z[z$id==i,"locus_tag"]<-
        gsub("(.*?\")(.*)", "\\2",  pre_locus)
      z[z$id==i, "direction"]<-
        ifelse(!any(grepl("complement", entry_head)), "leading", "compliment")
      try(z[z$id==i, "note"]<-
            gsub("  ","",paste0(entry[grep("note=", entry): 
                                        (grep("codon_start",entry)-1)], collapse="  ")),
          silent=T)
      
    } else if(any(grepl("tRNA", entry_head))){
      z[z$id==i, "type"]<-"tRNA"
      pre_locus<-gsub(paste("(.*\\s*\\",locus_tag,"=)(.*)(\")", sep=""),"\\2", 
                      grep(locus_tag,entry, value=T )[1])
      z[z$id==i,"locus_tag"]<-
        gsub("(.*?\")(.*)", "\\2",  pre_locus)
      z[z$id==i, "direction"]<-
        ifelse(!any(grepl("complement", entry_head)), "leading", "compliment")
      try(z[z$id==i, "product"]<-
            gsub("(.*\\s*/product=)(\")(.*)(\")","\\3", grep("product",entry, value=T )),
          silent=T)
      try(z[z$id==i, "inference"]<-
            gsub("(.*\\s*\\inference=)(*\")(.*)","\\3",
                 grep("inference",entry, value=T )), silent=T)
      try(z[z$id==i, "anticodon"]<-
            gsub("(.*\\s*/anticodon=)(.*)","\\2", 
                 grep("anticodon",entry, value=T )),silent=T)
      #####
      #   cant use the normal entry_head because occasionally other entries will 
      #   say "rRNA" in header
    } else if(any(grepl("rRNA", entry[c(1,3)]))){
      z[z$id==i, "type"]<-"rRNA"
      pre_locus<-gsub(paste("(.*\\s*\\",locus_tag,"=)(.*)(\")", sep=""),"\\2", 
                      grep(locus_tag,entry, value=T )[1])
      z[z$id==i,"locus_tag"]<-
        gsub("(.*?\")(.*)", "\\2",  pre_locus)
      z[z$id==i, "direction"]<-
        ifelse(!any(grepl("complement", entry_head)), "leading", "compliment")
      try(z[z$id==i, "note"]<-
            gsub("  ","",paste0(entry[grep("note=", entry): 
                                        (grep("codon_start",entry)-1)], collapse="  ")),
          silent=T)
      z[z$id==i, "product"]<-
        gsub("(.*\\s*/product=)(\")(.*)(\")","\\3", grep("product",entry, value=T ))
      #####
    } else if(any(grepl("misc_feature", entry_head))){
      z[z$id==i, "type"]<-"misc_feature"
      z[z$id==i, "direction"]<-
        ifelse(!any(grepl("complement", entry_head)), "leading", "compliment")
      z[z$id==i, "note"]<-
        gsub("  ","",paste0(entry[grep("note=\"", entry): (length(entry)-1)], collapse="  "))
      #####
    } else if(any(grepl("misc_binding", entry_head))){
      z[z$id==i, "type"]<-"misc_binding"
      z[z$id==i, "direction"]<-ifelse(!any(grepl("complement", entry_head)), 
                                      "leading", "compliment")
      try(z[z$id==i, "note"]<-
            gsub("  ","",
                 paste0(entry[grep("note=", entry):
                                (grep("codon_start",entry)-1)],collapse="  ")
            ), silent=T)
      try(z[z$id==i, "bound_moiety"]<-
            gsub("(.*\\s*/bound_moiety=\")(.*)(\")","\\2",
                 grep("bound_moiety",entry, value=T )), silent=T)
      #####
    } else if(any(grepl("assembly_gap", entry_head))){
      z[z$id==i, "type"]<-"assembly_gap"
      try(z[z$id==i, "estimated_length"]<-
            gsub("(.*\\s*/estimated_length=)(.*)","\\2",
                 grep("estimated_length",entry, value=T )), silent=T)
      try(z[z$id==i, "gap_type"]<-
            gsub("(.*\\s*/gap_type=\")(.*)(\")","\\2",
                 grep("gap_type",entry, value=T )), silent=T)
      try(z[z$id==i, "linkage_evidence"]<-
            gsub("(.*\\s*/linkage_evidence=\")(.*)(\")","\\2",
                 grep("linkage_evidence",entry, value=T )), silent=T)
      #####  
    } else if(any(grepl("source", entry_head))){
      z[z$id==i, "type"]<-"source"
      z[z$id==i, "source"]<-
        gsub("  ","",paste0(entry[grep("organism=", entry): 
                                    entry_end], collapse="  "))
    }
  }
  z<-z[1:(nrow(z)-1),]    #  bye bye dummy row
  z
}
##  IMPORTANT!!!  THIS RELIES ON THE RANGES DF HAVING THE LAST ROW AS BUFFER
#       Why is this important?  if you wanna subset, you must include 
#        the last row, either by using tail() or rbind(subset, tail(gl1,1))
#         if you dont, you will get an error.  
# 

#^^^^^^^^^ 
#ranges<-extract_features(data=returns[[2]],ranges =  ranges,locus_tag = "locus_tag")

################################################################################
# this function extracts the nucleotide sequence from the .gb file, along with 
#upstream and downstream regions of interest
get_seqs<-function(gbdf, seq, upstream=500, downstream=500){
  for( i in gbdf[gbdf$type !="source","id"]){ #print(i)}
    #gbdf["id"==i, "seq"]<- 
    gbdf[gbdf$id==i,"dnaseq"]<-
      substr(seq, gbdf[gbdf$id==i, "loc_start"], gbdf[gbdf$id==i, "loc_end"])
    gbdf[gbdf$id==i,"dnaseq_upstream"]<-
      substr(seq, gbdf[gbdf$id==i, "loc_start"]-upstream, 
             gbdf[gbdf$id==i, "loc_start"]-1)
    gbdf[gbdf$id==i,"dnaseq_downstream"]<-
      substr(seq, gbdf[gbdf$id==i, "loc_end"], 
             gbdf[gbdf$id==i, "loc_end"]+downstream)
    #gbdf["id"==i, "upstrea"]
  }
  gbdf
}
#^^^^^^^^^ 
#result<-get_seqs(gbdf = ranges, seq = returns[[3]], upstream = 500, downstream = 500)
