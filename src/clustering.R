
input_dir = "~/cluster_algorithms/"




# Input (named kegg_kegg)
load("~/input/input_TP_genesets.RData")

load("~/input/input_TP_genesets.RData")


######################################################################################################################
######################################################################################################################
###################################################### MGclus ########################################################



    

    keggclusters=list()
    j=1
    for(i in 1:length(kegg_kegg)){
      onekegg=kegg_kegg[[j]]
      write.table(x=onekegg,paste(input_dir,"onekegg.tsv",sep=""),sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE)
      system(paste("java -jar ",input_dir, "mgclusjar.jar  -f " ,input_dir, "onekegg.tsv -w T -o",input_dir, "/onekeggcluster  ",sep=""))
      
      # Read in the data
      x <- scan("~/cluster_algorithms/onekeggcluster", what="", sep="\n")
      # Separate elements by one or more whitepace
      y <- strsplit(x, "[[:space:]]+")
      
      keggclusters[[j]]=y
      j=j+1
    }
    
    #Naming each list in the list of pathways
    # names(keggclusters)=names(list_g)
    names(keggclusters)=names(kegg_kegg)
    #Each pathway module number
    groups=vector()
    for (i in 1:length(keggclusters)){
      nb=length(keggclusters[[i]])
      groups[i]=nb
    }
    
    
    modules1=data.frame(matrix(NA,nrow=length(unlist(keggclusters)),ncol=2))
    modules1=data.frame(matrix(NA,nrow=0,ncol=2))
    names(modules1)=c("gene","module")
    w=1
    for(i in 1:length(kegg_kegg)){
      if(groups[i]==0){
        
        i=i+0
      }
      else{
        
        for(k in 1:groups[i]){
          
          mod=unlist(keggclusters[[i]][k])
          mod=as.data.frame(mod)
          n=nrow(mod)
          mod1=data.frame(matrix(NA,nrow=n,ncol=2))
          
          names(mod1)<-c("gene","module")
          mod1[,1]=mod
          
          nombre = names(kegg_kegg[i])
          nombre = gsub('([[:punct:]])|\\s+','_',nombre)
          nombre = gsub("[[:blank:]]", "", nombre)
          pathgroup=rep(paste("g_",w,"_#",nombre,sep=""),n)
          pathgroup=trimws(pathgroup, which = "left")
          # pathgroup <- gsub('([[:punct:]])|\\s+','_',pathgroup)
          pathgroup <- sapply(pathgroup, toupper)
          w=w+1
          pathgroup=as.data.frame(pathgroup)
          mod1[,2]=pathgroup
          modules1=rbind(modules1,mod1)
          
        }
        
      }
      
      w=1
    }
    
    modules1 =   modules1[modules1[,2] %in% names(which(table(modules1[,2]) > 2)), ]
    
 
   
    write.table(modules1,"~/mgclus",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

    
  
  ######################################################################################################################
  #######################################################################################################################
  ###################################################### MCL ################################################################
  
    

    mcl_dir="~/cluster_algorithms/mcl/mcl-14-137/"
    kegg_kegg=kegg_kegg
    # kegg_kegg=list_g
    keggclusters=list()
    j=1
    setwd(paste(mcl_dir,"src/shmcl",sep=""))
    
    
    
    for(i in 1:length(kegg_kegg)){
      onekegg=kegg_kegg[[i]]
      write.table(x=onekegg,paste(input_dir,"onekegg.tsv",sep=""),sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE)
      system(paste("./mcl ",input_dir,"onekegg.tsv --abc -o ~/cluster_algorithms/onekeggcluster.mcl",sep=""))
      # Read in the data
      x <- scan("~/cluster_algorithms/onekeggcluster.mcl", what="", sep="\n")
      # Separate elements by one or more whitepace
      y <- strsplit(x, "[[:space:]]+")
      
      keggclusters[[j]]=y
      j=j+1
      
      
    }
    
    #Naming each list in the list of pathways
    # names(keggclusters)=names(list_g)
    names(keggclusters)=names(kegg_kegg)
    #Each pathway module number
    groups=vector()
    for (i in 1:length(keggclusters)){
      nb=length(keggclusters[[i]])
      groups[i]=nb
    }
    
    
    modules1=data.frame(matrix(NA,nrow=length(unlist(keggclusters)),ncol=2))
    modules1=data.frame(matrix(NA,nrow=0,ncol=2))
    names(modules1)=c("gene","module")
    w=1
    for(i in 1:length(kegg_kegg)){
      if(groups[i]==0){
        
        i=i+0
      }
      else{
        
        for(k in 1:groups[i]){
          
          mod=unlist(keggclusters[[i]][k])
          mod=as.data.frame(mod)
          n=nrow(mod)
          mod1=data.frame(matrix(NA,nrow=n,ncol=2))
          
          names(mod1)<-c("gene","module")
          mod1[,1]=mod
          
          nombre = names(kegg_kegg[i])
          nombre = gsub('([[:punct:]])|\\s+','_',nombre)
          nombre = gsub("[[:blank:]]", "", nombre)
          pathgroup=rep(paste("g_",w,"_#",nombre,sep=""),n)
          pathgroup=trimws(pathgroup, which = "left")
          # pathgroup <- gsub('([[:punct:]])|\\s+','_',pathgroup)
          pathgroup <- sapply(pathgroup, toupper)
          w=w+1
          pathgroup=as.data.frame(pathgroup)
          mod1[,2]=pathgroup
          modules1=rbind(modules1,mod1)
          
        }
        
      }
      
      w=1
    }
    
    modules1 =   modules1[modules1[,2] %in% names(which(table(modules1[,2]) > 2)), ]
    
    
    
    write.table(modules1,"~/mcl",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
    
    
    ######################################################################################################################
    #######################################################################################################################
    ###################################################### INFOMAP    ################################################################
    
    
    
    library(igraph)
    # INfomap clustering method
    i = 1
    keggclusters=list()
    j=1
    for(i in 1:length(kegg_kegg)){
      
      genesetd = graph_from_data_frame(kegg_kegg[[i]][1:2], directed = FALSE, vertices = NULL)
      prueba = cluster_infomap(genesetd, e.weights = as.vector(unlist(kegg_kegg[[i]][3])))

      keggclusters[[j]] =  communities(prueba)
      j = j + 1
    }
    
    
    
    #Naming each list in the list of pathways
    # names(keggclusters)=names(list_g)
    names(keggclusters)=names(kegg_kegg)
    #Each pathway module number
    groups=vector()
    for (i in 1:length(keggclusters)){
      nb=length(keggclusters[[i]])
      groups[i]=nb
    }
    
    
    modules1=data.frame(matrix(NA,nrow=length(unlist(keggclusters)),ncol=2))
    modules1=data.frame(matrix(NA,nrow=0,ncol=2))
    names(modules1)=c("gene","module")
    w=1
    for(i in 1:length(kegg_kegg)){
      if(groups[i]==0){
        
        i=i+0
      }
      else{
        
        for(k in 1:groups[i]){
          
          mod=unlist(keggclusters[[i]][k])
          mod=as.data.frame(mod)
          n=nrow(mod)
          mod1=data.frame(matrix(NA,nrow=n,ncol=2))
          
          names(mod1)<-c("gene","module")
          mod1[,1]=mod
          
          nombre = names(kegg_kegg[i])
          nombre = gsub('([[:punct:]])|\\s+','_',nombre)
          nombre = gsub("[[:blank:]]", "", nombre)
          pathgroup=rep(paste("g_",w,"_#",nombre,sep=""),n)
          pathgroup=trimws(pathgroup, which = "left")
          # pathgroup <- gsub('([[:punct:]])|\\s+','_',pathgroup)
          pathgroup <- sapply(pathgroup, toupper)
          w=w+1
          pathgroup=as.data.frame(pathgroup)
          mod1[,2]=pathgroup
          modules1=rbind(modules1,mod1)
          
        }
        
      }
      
      w=1
    }
    
    modules1 =   modules1[modules1[,2] %in% names(which(table(modules1[,2]) > 2)), ]
    
    
    
    write.table(modules1,"~/infomap",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
    
    