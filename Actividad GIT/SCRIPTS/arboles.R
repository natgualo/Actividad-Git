
pdf("GRAFICOS/ARBOLES.pdf")

library(ggtree)
library(ape)
library(ips)
library(ggplot2)

tree <- read.tree("DATA/tree_newick.nwk")
tree
ggtree(tree, branch.length="none", color="blue", size=2, linetype=3)

tree <- read.beast("DATA/flu_tree_beast.tree")

# supply a most recent sampling date so you get the dates
# and add a scale bar
ggtree(tree, mrsd="2013-01-01") + 
  theme_tree2() 

# Finally, add tip labels and adjust axis
ggtree(tree, mrsd="2013-01-01") + 
  theme_tree2() + 
  geom_tiplab(align=TRUE, linesize=.5) + 
  xlim(1990, 2020)

msaplot(p=ggtree(tree), fasta="DATA/flu_aasequence.fasta", window=c(150, 175))

## arbol 2

set.seed(42)
trees <- lapply(rep(c(10, 25, 50, 100), 3), rtree)
class(trees) <- "multiPhylo"
ggtree(trees) + facet_wrap(~.id, scale="free", ncol=4) + ggtitle("Many trees. Such phylogenetics. Wow.")


## Le gusta a nata
nwk <- system.file("extdata", "sample.nwk", package="treeio")

tree <- read.tree(nwk)
circ <- ggtree(tree, layout = "circular")

df <- data.frame(first=c("a", "b", "a", "c", "d", "d", "a", 
                         "b", "e", "e", "f", "c", "f"),
                 second= c("z", "z", "z", "z", "y", "y", 
                           "y", "y", "x", "x", "x", "a", "a"))
rownames(df) <- tree$tip.label

df2 <- as.data.frame(matrix(rnorm(39), ncol=3))
rownames(df2) <- tree$tip.label
colnames(df2) <- LETTERS[1:3]

dev.off()

#emojiiii 2
p <- ggtree(tree, layout = "circular", size=1) +  
  geom_tiplab(parse="emoji", size=10, vjust=.25)
print(p)

## fan layout  
p2 <- open_tree(p, angle=200) 
print(p2)

p2 %>% rotate_tree(-90)


## emojiiiii 3

library(ggtree)
tree_text <- paste0("(((((cow, (whale, dolphin)), (pig2, dog)),",
                    "camel), fish), seedling);")
x <- read.tree(text=tree_text)
library(ggimage)
p <-  ggtree(x, size=2) + geom_tiplab(size=20, parse='emoji') +
  xlim(NA, 7) + ylim(NA, 8.5) 

svglite::svglite("emoji.svg", width = 10, height = 7)
print(p)
