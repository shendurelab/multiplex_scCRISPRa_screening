##Loading required libraries

library(tidyverse)


# Create Data on gRNA categories for pie chart

data <- data.frame(
  gRNA_class=as.factor(c("NDD promoter (n=313)", "TSS positive control (n=30)", "K562 candidate enhancer hit (n=50)", "K562 candidate enhancer non-hit (n=50)", "Non-targeting control (n=50)")),
  value=c(313, 30, 50, 50, 50))

data$gRNA_class <- factor(data$gRNA_class, levels = c("NDD promoter (n=313)", "TSS positive control (n=30)", "K562 candidate enhancer hit (n=50)", "K562 candidate enhancer non-hit (n=50)", "Non-targeting control (n=50)"))
data$gRNA_class

# Piechart
ggplot(data, aes(x="", y=value, fill=gRNA_class)) +
  geom_bar(stat="identity", width=1, color="white") +
  scale_fill_manual(values = c("#E69F00", "#00B06B", "#56B4E9", "#56B4E9", "#999999")) +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
