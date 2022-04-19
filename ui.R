plotTS<-function(x,fs){
  
  # par(mar=c(1,1,1,1))
  autoplot(x[c(1:(60*fs)),],facet=F,ylim=c(-0.2,0.2))+
    theme_minimal()+ theme(legend.position = "bottom",legend.title = element_blank())
  # autoplot(head(as.xts(as.ts(x,frequency=fs)),60*fs),facet=F,ylim=c(-0.2,0.2))+theme_minimal()+ theme(legend.position = "bottom",legend.title = element_blank())
  
}


downloadButtonRmd <- function (outputId, label = "Download", class = NULL, ...)  {
  tags$a(id = outputId, class = paste("btn btn-default shiny-download-link", 
                                      class), href = "", target = "_blank", download = NA, 
         icon("download"), label, ...)
}

downloadLink<-function(outputId, label = "Download", class = NULL){
  tags$a(id = outputId,
         class = paste("btn btn-default shiny-download-link",class),
         href = "", target = "_blank", download = NA,icon("download"))}

int_to_unit <- function (x, adjustment=2^16) {
  x <- as.numeric(x)
  signs <- sign(x)
  x[signs < 0] <- x[signs < 0] + adjustment
  x
}


findexg<-function(nodes,edges,pos){
  res=as.character(nodes[nodes$id %in% subset(edges,to==pos)$from& nodes$group %in% "physical","label"])
  res<-unique(c(res,as.character(nodes[nodes$id %in% subset(edges,from==pos)$to& nodes$group %in% "physical","label"])))
  return(res)
}

findneighbors<-function(nodes,edges,nom,dir="both",g=NULL){
  if(dir=="both"){
    ids=nodes[nodes$label %in% nom,"id"]
    es=edges[edges$from %in% ids | edges$to %in% ids,]
    if(!is.null(g)){
      gnods=nodes[nodes$group %in% g,"id"]
      gedges=edges$id[(edges$from %in% gnods)|(edges$to %in% gnods)]
      redges=es[es$id %in% gedges,]
      res=unique(nodes[nodes$id %in% redges$from | nodes$id %in% redges$to,"label"])
      res=res[!(res %in% nom)]
    }else{
      res=unique(nodes[nodes$id %in% es$from | nodes$id %in% es$to,"label"])
      res=res[!(res %in% nom)]
    }
    # res=as.character(nodes[nodes$id %in% subset(edges,to==pos)$from & nodes$group %in% "physical","label"])
  }
  if(dir=="to"){
    ids=nodes[nodes$label %in% nom,"id"]
    es=edges[ edges$to %in% ids,]
    if(!is.null(g)){
      gnods=nodes[nodes$group %in% g,"id"]
      gedges=edges$id[(edges$from %in% gnods)]
      redges=es[es$id %in% gedges,]
      res=unique(nodes[ nodes$id %in% redges$from,"label"])
      res=res[!(res %in% nom)]
    }else{
      res=unique(nodes[nodes$id %in% es$from,"label"])
      res=res[!(res %in% nom)]
    }
  }
  if(dir=="from"){
    ids=nodes[nodes$label %in% nom,"id"]
    es=edges[edges$from %in% ids ,]
    if(!is.null(g)){
      gnods=nodes[nodes$group %in% g,"id"]
      gedges=edges$id[(edges$to %in% gnods)]
      redges=es[es$id %in% gedges,]
      res=unique(nodes[nodes$id %in% redges$to ,"label"])
      res=res[!(res %in% nom)]
    }else{
      res=unique(nodes[nodes$id %in% es$to ,"label"])
      res=res[!(res %in% nom)]
    }
  }
  if(length(res)==0){res=""}
  # res=as.character(nodes[nodes$id %in% subset(edges,to==pos)$from& nodes$group %in% "physical","label"])
  # res<-unique(c(res,as.character(nodes[nodes$id %in% subset(edges,from==pos)$to &
  #                                        nodes$group %in% "physical","label"])))
  return(res)
}

filter_value=function(nodes,edges,s){
  fn=nodes[nodes$label==s,"id"]
  fe=edges[edges$from==fn,"title"]
  res=strsplit(fe,":")
  res2=unlist(lapply(res,function(x)x[2]))
  res=unlist(lapply(res,function(x)x[1]))
  ## Making use of to is channel and from is filter
  if("High Pass" %in% res){
    n_f=length(edges[edges$to==fn,"title"])
    hf=rep(paste0("H:",res2[res %in% "High Pass"]),n_f)
  }else{
    n_f=length(edges[edges$to==fn,"title"])
    hf=rep(paste0("H:",0),n_f)
  }
  if("Low Pass" %in% res){
    n_f=length(edges[edges$to==fn,"title"])
    hf2=rep(paste0("L:",res2[res %in% "Low Pass"]),n_f)
  }else{
    n_f=length(edges[edges$to==fn,"title"])
    hf2=rep(paste0("L:",0),n_f)
  }
  if("Notch" %in% res){
    n_f=length(edges[edges$to==fn,"title"])
    hf3=rep(paste0("N:",res2[res %in% "Notch"]),n_f)
  }else{
    n_f=length(edges[edges$to==fn,"title"])
    hf3=rep(paste0("N:",0),n_f)
  }
  
  
  return(paste0(hf3,hf2,hf))
}


findrecons<-function(nodes,edges,pos){
  res=as.character(nodes[nodes$id %in% subset(edges,to==pos)$from& nodes$group %in% "Reconstruction","label"])
  res<-unique(c(res,as.character(nodes[nodes$id %in% subset(edges,from==pos)$to& nodes$group %in% "Reconstruction","label"])))
  return(res)
}

findcluster<-function(nodes,edges,pos){
  res=as.character(nodes[nodes$id %in% subset(edges,to==pos)$from& nodes$group %in% "Cluster","label"])
  res<-unique(c(res,as.character(nodes[nodes$id %in% subset(edges,from==pos)$to& nodes$group %in% "Cluster","label"])))
  return(res)
}

findcomponent<-function(nodes,edges,pos){
  res=as.character(nodes[nodes$id %in% subset(edges,to==pos)$from& nodes$group %in% "Components","label"])
  res<-unique(c(res,as.character(nodes[nodes$id %in% subset(edges,from==pos)$to& nodes$group %in% "Components","label"])))
  return(res)
}
