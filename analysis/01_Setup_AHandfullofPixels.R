# Windows is sensitive to the order of installation
# of some of the packages so follow this order 
# as noted below, on other systems this might be
# less of an issue but follow the instructions
# regardless

# First install the rlang packages which provides
# R functionalities needed by/for other packages
# in particular dplyr 
install.packages("rlang")

# next install dplyr needed to manipulate
# tidy data as discussed in the R crash course
install.packages("dplyr")

# for convenient plotting we also install
# ggplot2 (this is not required but most plots shown
# are generated using ggplot)
install.packages("ggplot2")
install.packages("patchwork")
install.packages("hexbin")
install.packages("leaflet")

# to render markdown documents, install these packages
install.packages("rmarkdown", "knitr") 

# to deal with relative paths in R markdown
# install the here package
install.packages("here")

# next up are the installs of the geospatial
# libraries, these can be installed in one
# go
install.packages(c("terra","tidyterra","sf","MODISTools"))
install.packages("geodata")
install.packages("signal")

install.packages("phenocamr")
install.packages("daymetr")
install.packages("rsample")
install.packages("parsnip")
install.packages("workflows")
install.packages("tune")
install.packages("dials")
install.packages("caret")
install.packages("GenSA")
install.packages("appeears")
install.packages("xgboost")