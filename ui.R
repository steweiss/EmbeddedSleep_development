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