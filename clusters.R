#data('hpy',package='cMonkey.data')
require(cMonkey)
#source('~/scratch/biclust/cmonkey.R', chdir=T)
debug.on()

require(parallel)
options(mc.cores=6)

if ( ! exists( 'organism.dir' ) ) {
  if ( exists( 'organism' ) ) organism.dir <- organism
  else organism.dir <- 'Eco'
}
organism <- paste(tolower(substr(organism.dir,1,1)),substr(organism.dir,2,3),sep='')

x=read.delim(sprintf("output/%s_clusters.tsv", organism.dir))
if ( ! exists('ratios') ) ratios <- sprintf('~/python/pynkey/%s/ratios.tsv.gz', organism.dir)
if ( ! exists('n.motifs') ) n.motifs <- 2
e=cmonkey.init(organism=organism, bg.order=3, k.clust=nrow(x), seed.method=c( rows="rnd", cols="rnd" ),
  discard.genome=F, parallel.cores=options('mc.cores')$mc.cores, parallel.cores.motif=options('mc.cores')$mc.cores)
e$cmonkey.re.seed( e )
sys.source("~/scratch/biclust/cmonkey-funcs.R",envir=e,chdir=T)
e$row.score.func='default'
e$iter = 1999
seq.type = 'upstream meme'
e$mast.cmd[1]='./progs//mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 99999 -mev 99999 -mt 0.99 -seqp -remcorr'
names(e$ratios)='ratios' ## for some reason gets set to 'ratios.1'
all.seqs = e$genome.info$all.upstream.seqs[[ seq.type ]]
bg.list = e$genome.info$bg.list[[ seq.type ]]
bg.fname = e$genome.info$bg.fname[ seq.type ]

#tmp <- e$motif.all.clusters( 1:nrow(x), seq.type=names(e$meme.cmd)[1], verbose=T ) ##strsplit( i, " " )[[ 1 ]][ 1 ],
#e$meme.scores[[1]] = tmp
#rm(tmp);gc()

##for (i in 1:nrow(x)) {
tmp = mclapply(1:nrow(x), function(i){
  cat(i," ",x$k[i]," ",x$resid[i],"\n")
  rows=strsplit(as.character(x$rows[i]),',')[[1]]
  cols=strsplit(as.character(x$cols[i]),',')[[1]]
  clust = e$clusterStack[[i]]
  clust$rows = rows
  clust$cols = cols
  clust$resid=x$resid[i]
  meme.out = strsplit(as.character(x$meme_out[i]),"<<<<>>>>")[[1]]
  if ( length( meme.out ) > 0 && ! all(is.na(meme.out)) ) clust$meme.out = e$getMemeMotifInfo( meme.out )
  else clust$meme.out = NULL
  ##print(clust$meme.out)

  mast.out <- list();
  pv.ev = list()
  if ( ! is.null( clust$meme.out ) ) {
    mast.out <- try( e$runMast( meme.out, e$mast.cmd[ seq.type ], names( all.seqs ), all.seqs,
                               verbose=TRUE, seq.type=seq.type, bg.list=bg.list, bgfname=bg.fname ) ) 
    if ( class(mast.out) != 'try-error' ) pv.ev <- e$get.pv.ev.single( mast.out, rows )
  }
  list(clust=clust, pv.ev=pv.ev) ## mast.out=mast.out, ## mast.out is too big!
} )
#rm(x,all.seqs,clust,bg.list,tmp2,bg.fname); gc()

for ( i in 1:length(tmp) ) {
  clust=tmp[[i]]$clust
  ##mast.out=tmp[[i]]$mast.out
  ##pv.ev=tmp[[i]]$pv.ev
  
  e$meme.scores[[seq.type]][[i]] = list( k=i, last.run=FALSE, meme.out=clust$meme.out, pv.ev=tmp[[i]]$pv.ev, prev.run=NULL )

  e$clusterStack[[i]] = clust
  if ( length(clust$meme.out) > 0 ) {
    tmp2 = e$cluster.pclust(i)
    clust$e.val <- tmp2$e.vals
    clust$p.clust <- tmp2$p.clusts
    e$clusterStack[[i]] = clust
  }
}

rm(tmp); gc()

e$row.score.func <- 'orig'
e$meme.scores[[seq.type]]$all.pv <- e$make.pv.ev.matrix( e$meme.scores[[seq.type]] )

save( e, file=sprintf('output/%s_out.RData', organism.dir) )
