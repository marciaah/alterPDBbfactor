#' @title getPDBinfo.R
#' @description gets the information of the PDB of interest, and identifies those regions that have missing coordinates
#' @param pdbid a PDB identifier (e.g. 2GS6) for only one PDB or the word "all" for all the available PDBs that correspond to the UniProtAC
#' @param uniprotAcc the uniProtKB accession code of your protein (e.g. P00533)
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @importFrom itertools isplitRows
#' @return  A dataframe with PDB information
#' @export
getPDBinfo<-function(pdbid,uniprotAcc){
  #DEBUG
  #uniprotAcc<-"P51587"
  #pdb_id<-"1N0W"
  #chain_id<-"B"
  pdbid<-tolower(pdbid)
  uniprotAcc<-gsub("-\\d+","",uniprotAcc)
  pattern<-gsub("-",".",uniprotAcc)

  #getting the PDBs, ordered by best coverage and resolution for your uniprotAcc of interest
  requestURL <- paste("https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/",uniprotAcc,sep="")
  r <- suppressMessages(httr::GET(requestURL, httr::accept("application/json")))

  #If the call to the server worked well
  if((! suppressMessages(grepl("error|Error|not found|Not Found|exceed|Status: 404|Service Unavailable",r))) &&
     r$status_code!="404" &&
     !purrr::is_empty(r)){
    pdbs <- as.data.frame(jsonlite::fromJSON(jsonlite::toJSON(suppressMessages(httr::content(r))))) ###ACA
    names(pdbs)=gsub(pattern = paste(pattern,".",sep=""),replacement = "",x=names(pdbs))
    pdbs$tax_id<-NULL
    pdbs$resolution[pdbs$experimental_method=="Solution NMR"]<-NA
    if(purrr::is_empty(pdbs$coverage)){pdbs$coverage<-NA}


    pdbs<-data.frame(lapply(pdbs, function(x) unlist(x)))
    pdbcoverage<-data.frame(obs_start=NA, obs_end=NA,  author_start=NA,author_end=NA, pdb_id=NA,chain_id=NA,struct_asym_id_chain_id=NA)

    #If "all" the PDBs were requested, then obtain the list from what is associated to the UniProt
    #Otherwise, the pdb_id to alter will be only the one that was requested
    if(pdbid=="all"){ pdb_list<-unique(pdbs$pdb_id) }else{pdb_list<-pdbid}


    for(pdb_id in pdb_list){

      #WARNING: PDBe GRAPH database Rest API sometimes does not return "entity_id" and "preferred_assembly_id",
      #Hence, two strategies are followed here depending on if this information was present or not
      if( "entity_id" %in% colnames(pdbs) & "preferred_assembly_id" %in% colnames(pdbs))
      {target<-pdbs[pdbs$pdb_id==pdb_id,c("pdb_id","chain_id","entity_id","preferred_assembly_id")]
      with_entity_id=T}else
      { target<-pdbs[pdbs$pdb_id==pdb_id,c("pdb_id","chain_id")]
      with_entity_id=F
      pdbs$preferred_assembly_id<-NA
      }


      for (chain_id in target$chain_id){
        #Observed ranges in a PDB for when there are more than 1 entities, so you get all chains in all entities and can select from there
        requestURL <- paste("https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/",pdb_id,sep="")

        r <- httr::GET(requestURL, httr::accept("application/json"))

        if((! suppressMessages(grepl("error|Error|not found|exceed|Service Unavailable|Status: 404",r))) && !purrr::is_empty(r) && r$status_code!=404){
          coverage <- as.data.frame(jsonlite::fromJSON(jsonlite::toJSON(suppressMessages(httr::content(r)))))
          colnames(coverage)<-c("entity","chains")
          if(with_entity_id) #if with_entity_id==T  (there were entity_id coming from PDBe Graph database)
          {
            coverage<-coverage[coverage$entity %in% target$entity_id[target$chain_id==chain_id],]
            #If >1 entities
            for (entity in unique(target$entity_id[target$chain_id==chain_id])) {
              chains<-as.data.frame(coverage$chains[coverage$entity==entity])
              chains$struct_asym_id<-unlist(chains$struct_asym_id)
              #WARNING: the chain_id that comes from uniProtKB API seems to refer to the chains in "struct_asym_id", assuming this in the following lines
              chains<-chains[chains$struct_asym_id == chain_id,]
              observed<-as.data.frame(chains$observed)
              obs_start<-unlist(observed$start$residue_number)
              #obs_start<-unlist(chains$start$residue_number)
              obs_end<-unlist(observed$end$residue_number)
              #obs_end<-unlist(chains$end$residue_number)
              author_start<-unlist(observed$start$author_residue_number)
              author_end<-unlist(observed$end$author_residue_number)

              #correcting the chain_id, the one that comes from UniProtKB Proteins API corresponds to the chain in the preferred asymmetric unit
              struct_asym_id_chain_id<-chain_id
              chain_id<-chains$chain_id[chains$struct_asym_id==struct_asym_id_chain_id]

              #If 1 entity
              #Sometimes the chain_id for the pdb_id obtained requesting the PDBe API by UniProtAcc, corresponds to struct_asym_id in this table
              #e.g. has chain_id="P" and struct_asym_id="B"
              #entity<-as.data.frame(coverage$chains[coverage$entity==unique(as.numeric(pdbs$entity_id[pdbs$pdb_id==pdb_id]))])
              #obs_start<-unlist(as.data.frame(entity$observed[entity$chain_id==as.character(pdbs$chain_id[pdbs$pdb_id==pdb_id]) || entity$struct_asym_id==as.character(pdbs$chain_id[pdbs$pdb_id==pdb_id])])$start$residue_number)
              #obs_end<-unlist(as.data.frame(entity$observed[entity$chain_id==as.character(pdbs$chain_id[pdbs$pdb_id==pdb_id]) || entity$struct_asym_id==as.character(pdbs$chain_id[pdbs$pdb_id==pdb_id])])$end$residue_number)
              message("rbindinggg")
              pdbcoverage<-as.data.frame(rbind(pdbcoverage,as.data.frame(cbind(obs_start,obs_end,author_start,author_end,pdb_id=rep(pdb_id,length(obs_end)),chain_id=rep(chain_id,length(obs_end)),struct_asym_id_chain_id=rep(struct_asym_id_chain_id,length(obs_end))))))
              message("rbindindeeedd")
            }
          }else{ #else, if with_entity_id==F  (no entities came from PDBe Graph database)
            for(e in 1:length(coverage$chains)){ #going through the possible entities and checking coincidence of chain_id
              chains<-as.data.frame((coverage[coverage$entity==e,]$chains))
              chains<-chains[chains$chain_id==chain_id,]
              chains$struct_asym_id<-unlist(chains$struct_asym_id)
              observed<-as.data.frame(chains$observed)
              obs_start<-unlist(observed$start$residue_number)
              obs_end<-unlist(observed$end$residue_number)
              author_start<-unlist(observed$start$author_residue_number)
              author_end<-unlist(observed$end$author_residue_number)

              #correct the chain_id, the one that comes from UniProtKB Proteins API corresponds to the chain in the preferred asymmetric unit
              struct_asym_id_chain_id<-chain_id
              #chain_id<-chains$chain_id[chains$struct_asym_id==struct_asym_id_chain_id]

              tryCatch(pdbcoverage<-as.data.frame(rbind(pdbcoverage,as.data.frame(cbind(obs_start,obs_end,author_start,author_end,pdb_id=rep(pdb_id,length(obs_end)),chain_id=rep(chain_id,length(obs_end)),struct_asym_id_chain_id=rep(struct_asym_id_chain_id,length(obs_end)))))),
                       error=function(e){
                         message(paste("\nWARNING: Failed obtaining PDB coverage info for ",pdb_id,":",chain_id,". rbind() returned:"))
                         message(e)
                         message("\n\nStructure will not be downloaded and PyMOL script will not be generated")
                       }
              )
            }
          } #end if-else with_entity_id
        } #end if polymer_coverage was successful
      } #end going by chain
    } # end going by pdb
    #remove NAs line

    #THis is if filtering to only one pdb_id of interest
    pdbcoverage<-pdbcoverage[!is.na(pdbcoverage$pdb_id),]
    #pdbs<-pdbs[pdbs$pdb_id==tolower(pdb_id),]
  }else{ #If the call to the PDBe server DID NOT worked well
    stop("EROR!: Failed in retrieving protein PDB structural information from PDBe API\nNo structural coverage available for transferring your scores to the PDB structure :( \n") }


  pdblist<-list(pdbs,pdbcoverage)
  return(pdblist)
}


#' @title alterPDBbfactor.R
#' @description Download and renumber a PDB structure, then alters the B-factor column replacing by scores/numbers of your interest
#' @param pdbid a PDB identifier (e.g. 2GS6)
#' @param chain_id a chain of yoyr PDB of interest (e.g. A)
#' @param uniprotAcc the uniProtKB accession code of your protein (e.g. P00533)
#' @param seq_scores Path to a file where you have two columnsA dataframe with two columns: residue position (aacPos) number and the by-residue scores (score) you want to add in the B-factor field of the PDB (e.g. conservation)
#' aacPos | score
#' 1      |  0.1
#' 2      |  1
#' 3      |  0.5
#'
#' @param outfolder The output folder
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @importFrom itertools isplitRows
#' @return  The PyMol script
#' @export
alterPDBbfactor<-function(pdbid,uniprotAcc,seq_scores,outfolder){

  #DEBUG:------------
  #outfolder<-"./out"
  #pdbid<-"all"
  #uniprotAcc<-"P51587"
  #simulated numberx
  #seq_scores<-data.frame(aacPos=c(1:3418),score=c(round(runif(3418,min=0,max=1),2)))
  #------------------
  seq_scores<-data.table::fread(input=seq_scores,sep = "\t",header = T)
  pdbid<-tolower(pdbid)


  #if outfolder does not exist, create it
  if (!dir.exists(outfolder)){
    dir.create(outfolder)
  }

  #Attempt to obtain the coverage of the PDB over the UniProtKB sequence
  pdb_info<-getPDBinfo(pdb = pdbid,uniprotAcc = uniprotAcc)
  pdb<-as.data.frame(pdb_info[[1]])
  pdbcoverage<-as.data.frame(pdb_info[[2]])




  aux<-subset(unique(as.data.frame(tidyr::separate_rows(pdb,chain_id,sep=", "))),
              select=c("chain_id","pdb_id","start","end","unp_start","unp_end"))

  aux<-unique(merge(x=aux,y=pdbcoverage[,c("pdb_id","chain_id","struct_asym_id_chain_id")],by.x=c("pdb_id","chain_id"),by.y=c("pdb_id","struct_asym_id_chain_id"),all.x=T))

  #THis is to check if the PDB chain has insertion or deletions.
  #If so, the check column will have 1, and this structure will not be straightforward to alter
  #TODO: include this structures, but put different color for this regions
  aux$shift_start<-aux$unp_start-aux$start
  aux$shift_end<-aux$unp_end-aux$end
  aux$check<-ifelse(aux$shift_start==aux$shift_end,0,1)

  #This dataframe will have by-site mapping of the score from uniprot coordinates to PDB coordinates
  pdbinfo<-data.frame(pdb_id=NA,chain_id=NA,pdb_pos=NA,unip_pos=NA,pdb_shift_pos=NA,author_pos=NA,score=NA)

  for(pdb_id in unique(aux$pdb_id))
  {

    #DEBUG:
    #pdb_id<-"8c3n"

    if(unique(aux$check[aux$pdb_id==pdb_id]==0)){

      block<-data.frame(pdb_id=NA,chain_id=NA,pdb_pos=NA,unip_pos=NA,pdb_shift_pos=NA,author_pos=NA,score=NA)
      auxpdb<-aux[aux$pdb_id==pdb_id,]
      for (row in 1:nrow(auxpdb)) {
        #DEBUG:
        #row<-1

        # the chain corresponding to struct_asym_id
        chain_id<-auxpdb$chain_id[as.numeric(row)]

        cat(paste("Doing ",pdb_id," ",chain_id,"\n"))

        # positions present in the PDB
        pos<-NULL
        author_pos<-NULL
        for(region in row.names(pdbcoverage[pdbcoverage$pdb_id==pdb_id & pdbcoverage$struct_asym_id_chain_id==chain_id,])){
          #DEBUG:
          #message(region)

          pos<-c(pos,as.numeric(pdbcoverage$obs_start[as.numeric(region)-1]) : as.numeric(pdbcoverage$obs_end[as.numeric(region)-1]))
          author_pos<-c(author_pos, as.numeric(pdbcoverage$author_start[as.numeric(region)-1]) : as.numeric(pdbcoverage$author_end[as.numeric(region)-1]))
        }
        covered_pos<-sort(pos)
        covered_author<-sort(author_pos)
        covered<-as.data.frame(cbind(covered_pos,covered_author))
        jump<-0

        # This is to map the scores from the uniProt to the amino acids in PDB structure
        for(pospdb in auxpdb$start[row]:auxpdb$end[row]) #positions in pdb numbering
        {
          #DEBUG
          #pospdb<-1

          pdb_pos<-pospdb
          author_pos<-covered$covered_author[covered$covered_pos==pospdb]
          pdb_shift_pos<-auxpdb$unp_start[row]+jump
          unip_pos<-auxpdb$unp_start[row]+jump
          jump<-jump+1
          score<-unique(seq_scores$score[seq_scores$aacPos==unip_pos])

          #Checking whether the position in PDB has coordinates solved
          if(pdb_pos %in% covered$covered_pos)
          {
            #DEBUG:
            #message(paste(pdb_id,chain_id,pdb_pos,unip_pos,pdb_shift_pos,author_pos,colour,group,conservation),sep=",")
            pdbinfo<-rbind(pdbinfo,cbind(pdb_id,chain_id,pdb_pos,unip_pos,pdb_shift_pos,author_pos,score))
            #pdbinfo<-rbind(pdbinfo,cbind(pdb_id,ch,pdb_pos,unip_pos,pdb_shift_pos,author_pos,score))
          }
        }
      }

      #just in case some NAs are thre
      pdbinfo<-unique(pdbinfo[!is.na(pdbinfo$pdb_id),])


      ## Generating a PDB file with conservation score in B-factors column
      # Also the amino acids are renumbered to match the UniProt sequence numbering
      # This PDB structure will be used for coloring by conservation and by CCRs
      skip_to_next <- FALSE
      #attempt to obtain the PDB structure
      tryCatch(pdbstructure <- bio3d::read.pdb(pdb_id,rm.alt=FALSE,verbose = F),
               error = function(e){
                 message(paste("FAILED... ",pdb_id," not downloaded. Maybe the PDB is deprecated? Please, check in https://www.ebi.ac.uk/pdbe/ or https://www.rcsb.org/"))
                 skip_to_next<-TRUE
               } )

      # Skipping to next pdb in loop if an error occurred
      if(skip_to_next) { next }

      #Save the original PDB structure in PDB format. ONLY SAVES ATOM and HETATM records
      bio3d::write.pdb(pdbstructure, file=paste(outfolder,"/",pdb_id,"_original.pdb",sep=""))

      # extract the atom records. Among the columns, "b" has the b-factors
      atom<-pdbstructure$atom

      #adding the score to the atom in structure
      tryCatch( mergeblock<-merge(x=atom,y=pdbinfo[pdbinfo$pdb_id==pdb_id,],by.x=c("chain","resno"),by.y=c("chain_id","author_pos"),all.x=F),
                error = function(e){
                  message(paste("FAILED... I was not able to merge your scores into the PDB:",pdb_id,"\nI'm still not able to process this situations.\nPlease, manually check UniProt vs PDB sequence numbering to find any clues\n"))
                  skip_to_next<-TRUE
                } )

      # Skipping to next pdb in loop if an error occurred
      if(skip_to_next) { next }

      mergeblock$b<-ifelse(is.na(mergeblock$score),-1,mergeblock$score)


      shift<-unique(as.numeric(pdbinfo$pdb_shift_pos[pdbinfo$pdb_id==pdb_id]) - as.numeric(pdbinfo$author_pos[pdbinfo$pdb_id==pdb_id]))
      tryCatch( mergeblock$resno<-mergeblock$resno[mergeblock$pdb_id==pdb_id]+shift,
                error = function(e){
                  cat(paste("FAILED... I was not able to merge your scores into the PDB:",pdb_id,"\nI'm still not able to process this situations.\nPlease, manually check UniProt vs PDB sequence numbering to find any clues\n\n"))
                  skip_to_next<-TRUE
                } )

      # Skipping to next pdb in loop if an error occurred
      if(skip_to_next) { next }

      #order by atom element
      mergeblock<-mergeblock[order(mergeblock[,4]),]
      atom<-atom[order(atom[,2]),]


      if (nrow(atom[atom$eleno %in% mergeblock$eleno,]) == length(mergeblock$eleno)){
        atom[atom$eleno %in% mergeblock$eleno,]<-mergeblock[,c("type","eleno","elety","alt","resid","chain","resno","insert","x","y","z","o","b","segid","elesy","charge")]

        # change the B-factors in the sites that do not have score assigned and other chains other than the one of your protein of interest
        atom$b[atom$b>1 | atom$b<0]<- -1
        pdbstructure$atom<-atom
      }else{
        stop("FAILED...was not able to renumber PDB structure",pdb_id," :(\n")

      }





      # write to PDB file the structure with altered B-factors only in the chain/s of interest
      # WARNING: only the ATOM and HETATM lines are printed out
      bio3d::write.pdb(pdbstructure, file=paste(outfolder,"/",pdb_id,"_Bfact-score.pdb",sep=""))



    }else{
      message(paste(" FAILED with UniProt: ",uniprotAcc," vs PDB: ",pdb_id," residue position correspondence. I'm still not able to process this situations.\nPlease, manually check UniProt vs PDB sequence numbering to find any clues\n" ,sep=""))
    }
  }


  cat("Done! :D (dismiss the Warning messages!) \n")
}
