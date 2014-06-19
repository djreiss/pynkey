if ( ! exists( 'organism' ) ) organism <- 'Hpy'
x=read.delim(sprintf("output/%s_stats.tsv", organism))
for(i in 1:ncol(x)){
  if(class(x[[i]])=='factor'){
    if(x[[i]][1]%in%c('true','false')) x[[i]]=as.logical(x[[i]])
    else x[[i]]=as.numeric(gsub('f0','',as.character(x[[i]])))
  }
}
par(mfrow=c(3,3))
try(plot(x$iter,x$RESID,typ='l'))
try(plot(x$iter,x$MEME_PVAL,typ='l'))
try(plot(x$iter,x$STRING_DENS,typ='l'))
try(plot(x$iter,x$ROWS,typ='l'))
try(plot(x$iter,x$COLS,typ='l'))
try(plot(x$iter,x$CLUSTS_PER_ROW,typ='l'))
try(plot(x$iter,x$CLUSTS_PER_COL,typ='l'))
try(plot(x$iter,x$N_MOVES,typ='l'))
try(plot(x$iter,x$N_IMPROVEMENTS,typ='l'))
