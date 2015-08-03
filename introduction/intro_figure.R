#Generating the tree from the intro figure
library(diversitree)


quartz(width=5.8 , height=8.3)
seed=20 ; set.seed(seed)
phy <- tree.musse(pars, 60, x0=1, include.extinct=TRUE)
h <- history.from.sim.discrete(phy, 1:3)
plot(h, phy, cex=.7, main=seed)

#Remove extinct
phy2 <- prune(phy)
h2 <- history.from.sim.discrete(phy2, 1:3)
quartz(width=5.8 , height=5.8)
plot(h2, phy2, cex=.7)

#Remove living (most)

#Select the livings
livings<-phy$tip.label[grep("sp", phy$tip.label)]
#Selecting an arbitrary number of living taxa to keep (11)
keep<-c("sp78","sp41","sp53","sp55","sp69","sp86","sp70","sp102","sp74","sp116","sp68")

livings[-c(match(keep, livings))]

to.drop.names<-livings[-c(match(keep, livings))]
to.drop<-match(to.drop.names, phy$tip.label)

phy3 <- drop.tip(phy, tip=livings[-c(match(keep, livings))])
phy3 <- prune(phy, to.drop=to.drop)
h3 <- history.from.sim.discrete(phy3, 1:2)
quartz(width=5.8 , height=5.8)
plot( phy3, cex=.7)
