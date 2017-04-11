#~ make a set of trees for model fitting
#~ merged crf02 and g
require(ape) 
require(rcolgem)

NTREES <- 20

tre_crf02 <- read.nexus ('trust_neuroaids_worldRef_datesAdded_drmmasked_crf02_f1e5.trees') 
tre_g <- read.nexus( 'trust_neuroaids_worldRef_datesAdded_drmmasked_G_f1e5.trees' )

.calcst <- function(tre)
{
	setNames( sapply( strsplit( tre$tip.label, split='_'), function(splits) as.numeric(splits[length(splits)] ) )
	  , tre$tip.label)
}
	
newtres <- as.list(1:NTREES)
#~ for (i in 1:NTREES)
lapply( 1:NTREES, function(x)
{
	i <- round( runif( 2, min = 500, max = 2000) )
	t0 <- tre_crf02[[i[1]]]
	t1 <- tre_g[[i[2]]]
	
	bdtCRF02 <- binaryDatedTree( t0, .calcst(t0) )
	bdtG <- binaryDatedTree( t1, .calcst(t1) )
	
	bdtG$root.edge <- max(bdtG$maxHeight,bdtCRF02$maxHeight) - bdtG$maxHeight+1
	bdtCRF02$root.edge <- max(bdtG$maxHeight,bdtCRF02$maxHeight) - bdtCRF02$maxHeight+1
	
	newnwk <- paste(sep='', '('
	  , paste(sep=',', sub(';', '', write.tree( bdtG )) , sub(';', '', write.tree(bdtCRF02)) )
	  , ');'
	)
print(date())
#~ 	newtres[[i]] <- read.tree(text = newnwk )
	read.tree(text = newnwk )
}) -> newtres
class(newtres) <- 'multiPhylo'

write.tree( newtres, file = 'mergeTrees0.nwk') 
