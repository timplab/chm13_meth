# library
library(tidyverse)
library(wesanderson)
library(ggsci)
library(scales)

# set theme
theme_set(theme_bw()+ 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background = element_blank(),
                panel.border = element_rect(size=0.5,colour="black"),
                axis.text = element_text(color="black")
                )
          )
pal <- wes_palette("FantasticFox1")
heat_pal <- c(pal[3],pal[2],pal[4],pal[5])
pal <- wes_palette("Rushmore1")

pal_cont_polar <- wes_palette("Zissou1", 21, type = "continuous")

# functions
plotter <- function(print_obj,pre = "temp",out_dir = "~/Dropbox/Data/tmp",h = 4, w = 6){
  # if not interactive, just print it
  if( ! interactive()) {
    print(print_obj)
  } else {
    # if interactive, I want to export this to a pdf file
    outpath <- paste0(out_dir,"/",pre,".pdf")
    message(paste("outputting plot to", outpath))
    pdf(outpath,height = h, width = w, useDingbats = F)
    print(print_obj)
    dev.off()
  }
}
