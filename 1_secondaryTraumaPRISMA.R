# https://cran.r-project.org/web/packages/PRISMAstatement/vignettes/PRISMA.html
# install.packages("PRISMAstatement")
library(PRISMAstatement)

prsm <- prisma(found = 882,
       found_other = 0,
       no_dupes = 882,
       screened = 882,
       screen_exclusions = 747,
       full_text = 135,
       full_text_exclusions = 86,
       labels = list(screen_exclusions="Records excluded based \non Title or Abstract\n(n=747)"),
       qualitative = 49,
       quantitative = 49,
       width = 800,
       height = 800)

pdf("STD.pdf")
prsm
dev.off()