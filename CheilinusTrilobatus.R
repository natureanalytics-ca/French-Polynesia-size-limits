
library(FishLife)
library(ggplot2)
library(stringr)
library(dplyr)


#-----------------------------
#Papeete market length data
#-----------------------------

#----------------
#species name

nm<-"Cheilinus trilobatus"

#-------------------------------------------------
#Function for getting mode from raw length data

Lc.func=function(x) {
  #  x is length sample.Returns Lc
  z=table(x)
  z1=cumsum(z)
  z1=z1/sum(z)
  a=as.numeric(names(z))
  a1=seq(trunc(min(a))+1,max(a),by=1)
  d=loess(z1~a)
  d1=predict(d,newdata=a1)
  d2=d1[2:length(a1)]-d1[1:(length(a1)-1)]
  d3=a1[1:length(d2)][d2==max(d2)]
  d3
}

sel.func=function(x){
  mode <- Lc.func(x)
  cumulativeProb<-cumsum(sort(x[x<=mode]))/sum(sort(x[x<=mode]))
  L50<-sort(x[x<=mode])[which(abs(cumulativeProb-0.5)==min(abs(cumulativeProb-0.5)))]
  L95<-sort(x[x<=mode])[which(abs(cumulativeProb-0.95)==min(abs(cumulativeProb-0.95)))]
  return(c(L50, L95))
}

#-------------------------------
#Load data - 2012 data

lengthData2012<-readr::read_csv("Raw data fish 2012.csv", col_names = TRUE)

lengthData2012_mod<-lengthData2012 %>%
  mutate("name" = str_c(Genus, Species, sep=" ")) %>%
  select(-`Coeff a`, -`Coeff b`, -`Weight (g)`)
#View(lengthData2012_mod)

#-------------------------------
#Load data - 2021 data

lengthData2021<-readr::read_csv("RAW-DATA_Market_2021.csv", col_names = TRUE)

lengthData2021_mod<-lengthData2021 %>%
  mutate("name" = str_c(Genus, Species, sep=" "))
#View(lengthData2021_mod)

#--------------------------
#Selectivity summaries

#Table 2012
length2012_sub<-lengthData2012_mod %>%
  filter(Site == "Papeete market") %>%
  filter(name %in% nm) %>%
  group_by(name, Family) %>%
  summarize(
    n =  n(), 
    min = min(`Fork length`, na.rm=TRUE), 
    max = max(`Fork length`, na.rm=TRUE),
    mode = Lc.func(`Fork length`),
    L50 = sel.func(`Fork length`)[1],
    L95 = sel.func(`Fork length`)[2]
  ) %>%
  arrange(Family, name)

#Table 2021
length2021_sub<-lengthData2021_mod %>%
  filter(name %in% nm) %>%
  group_by(name, Family) %>%
  summarize(
    n =  n(), 
    min = min(`Fork length`, na.rm=TRUE), 
    max = max(`Fork length`, na.rm=TRUE),
    mode = Lc.func(`Fork length`),
    L50 = sel.func(`Fork length`)[1],
    L95 = sel.func(`Fork length`)[2]
  ) %>%
  arrange(Family, name)

#Combined data sets
combined_sub<-lengthData2012_mod %>%
  filter(Site == "Papeete market") %>%
  bind_rows(lengthData2021_mod) 

combined_sub_table <- combined_sub %>%
  filter(name %in% nm) %>%
  group_by(name, Family) %>%
  summarize(
    n =  n(), 
    min = min(`Fork length`, na.rm=TRUE), 
    max = max(`Fork length`, na.rm=TRUE),
    mode = Lc.func(`Fork length`),
    L50 = sel.func(`Fork length`)[1],
    L95 = sel.func(`Fork length`)[2]
  ) %>%
  arrange(Family, name)

#-------------
#Create histo

#You can uncomment the png and dev.off to save the image to file.
png(file=paste0("selectivity_", nm, ".png"), width=6.5, height=10, units="in", res=300, bg="white", pointsize=12)

dt<-combined_sub %>%
  filter(name %in% nm)
dt<-dt$`Fork length`
L<-seq(min(dt, na.rm=TRUE)-1, max(dt, na.rm=TRUE)+1, 1)
X<-hist(dt, breaks = L, plot=TRUE, main = nm[i], xlab = "Length (cm)")
par(new = TRUE)
L50<-combined_sub_table$L50[which(combined_sub_table$name == nm[i])]
L95<-combined_sub_table$L95[which(combined_sub_table$name == nm[i])]
s<--(L95-L50)/log(1/0.95-1)
prob<-plogis(L, location=L50, scale=s)
plot(L, prob, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col = "blue" , lwd=3)
axis(side=4)
mtext("Selectivity", side=4, line=3, las = 3, cex = 0.8)

#dev.off()



#----------------
#FishLife
#----------------

Predictions <- Plot_taxa( Search_species(Genus="Cheilinus",Species="trilobatus")$match_taxonomy)

mean_L50<-exp(Predictions[[2]]$Mean_pred[7]+diag(Predictions[[2]]$Cov_pred/2)[7])
mean_Linf<-exp(Predictions[[2]]$Mean_pred[1]+diag(Predictions[[2]]$Cov_pred/2)[1])
mean_L50
mean_Linf
mean_L50/mean_Linf


mean_M<-exp(Predictions[[1]]$Mean_pred[6]+diag(Predictions[[1]]$Cov_pred/2)[6])
mean_K<-exp(Predictions[[1]]$Mean_pred[2]+diag(Predictions[[1]]$Cov_pred/2)[2])
mean_M/mean_K
