install.packages("devtools")
library(devtools)
install_github("didacs/ggsunburst")
library(ggsunburst)

# newick format
nw <- "(((a, b, c), (d, e, f, g)), (f, i, h));"

# extract data
sb <- sunburst_data(nw)
icicle(sb)
sunburst(sb)

# demonstrate use of data with ggplot()
p <- ggplot(sb$rects)

# icicle diagram
p + geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), size=1, color="white") + scale_y_reverse()

# sunburst diagram
p + geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), size=1, color="white") + coord_polar()