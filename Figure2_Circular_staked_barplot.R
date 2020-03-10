# Here the script to reproduce Figure 2 of the paper "Fast analysis of scATAC-seq data using a predefined set of genomic regions".
# Peaks with p-value<0.05 were selected from DHS500 clusters which had corresponding clusters after CellRanger count matrix analysis:
# 5 times these conditions were respected (cluster A, B, C, D, E).
# These clusters were processed with bedtools closest and the results were summurized into intervals. For each interval the % of peaks at that distance is represented,

# library
library(tidyverse)

# Create dataset
data2 <- data.frame(
  intervals=c('None','d=0','0<d<1000','1000<d<5000','5000<d<10000','10000<d<20000','d>20000',
                   'None','d=0','0<d<1000','1000<d<5000','5000<d<10000','10000<d<20000','d>20000',
                   'None','d=0','0<d<1000','1000<d<5000','5000<d<10000','10000<d<20000','d>20000',
                   'None','d=0','0<d<1000','1000<d<5000','5000<d<10000','10000<d<20000','d>20000',
                   'None','d=0','0<d<1000','1000<d<5000','5000<d<10000','10000<d<20000','d>20000'
                   ),
  group=c(rep('A',7),
          rep('B',7),
          rep('C',7),
          rep('D',7),
          rep('E',7)),
  value=c(0.0,0.8291956305858987,0.003972194637537239,0.011916583912611719,0.006951340615690168,0.007944389275074478,0.1400198609731877,
          0.0,0.7693836978131213,0.0,0.008946322067594433,0.005964214711729622,0.012922465208747515,0.20278330019880716,
          0.0,0.7966101694915254,0.0,0.004985044865403789,0.006979062811565304,0.010967098703888335,0.18045862412761715,
          0.0,0.8019900497512438,0.001990049751243781,0.006965174129353234,0.007960199004975124,0.006965174129353234,0.17412935323383086,
          0.0,0.5754245754245755,0.000999000999000999,0.006993006993006993,0.007992007992007992,0.007992007992007992,0.4005994005994006)
)
head(data2)

data2$value= data2$value+0.0001
data2$value= data2$value*100

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data2$group), ncol(data2)) )
colnames(to_add) <- colnames(data2)
to_add$group <- rep(levels(data2$group), each=empty_bar)
data2 <- rbind(data2, to_add)
data2 <- data2 %>% arrange(group)
data2$id <- seq(1, nrow(data2))

# Get the name and the y position of each label
label_data <- data2
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data2 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data2, aes(x=as.factor(id), y=value, fill=group)) +

  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=80/60/40/20/0 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # Add text showing the value of each 80/60/40/20/0 lines
  annotate("text", x = rep(max(data2$id),5), y = c(0, 20, 40, 60, 80), label = c("0","20", "40", "60", "80") , 
           color="grey", size=5 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,100) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=intervals, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=6, angle= label_data$angle, inherit.aes = FALSE) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
               colour = "black", alpha=0.8, size=0.3 , inherit.aes = FALSE)  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), 
            hjust=c(0.5,1,0.5,0.2,0), 
            colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)


p + ggtitle("CellRanger \nvs \nDHS500") +
  theme(plot.title = element_text(hjust = 0.5,vjust=-130,size=14))

ggsave("Bar_cicle", plot = last_plot(), device = "jpeg",
       scale = 1, width = 20, height = 10, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE)
