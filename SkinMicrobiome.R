library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(qiime2R)
library(ggtext)
library(RColorBrewer)
library(ggvenn)
library(broom)
library(ggbreak)
library(phyloseq)
library(vegan)
library(glue)
library(ANCOMBC)
1. Alpha diversity
metadata<-read_q2metadata("sample-metadata.tsv")
shannon <- read_tsv("shannon/alpha-diversity.tsv")
colnames(shannon)[colnames(shannon) == "...1"] <- "sample-id"
chao1 <- read_tsv("chao1/alpha-diversity.tsv")
colnames(chao1)[colnames(chao1) == "...1"] <- "sample-id"
faithpd <- read_tsv("faithpd/alpha-diversity.tsv")
colnames(faithpd)[colnames(faithpd) == "#SampleID"] <- "sample-id"
pielou <- read_tsv("pielou/alpha-diversity.tsv")
colnames(pielou)[colnames(pielou) == "...1"] <- "sample-id"
alphadiversity <- merge(shannon, faithpd, by = "sample-id", all =  TRUE)
alphadiversity <- merge(alphadiversity, chao1, by = "sample-id", all =  TRUE)
alphadiversity <- merge(alphadiversity, pielou, by = "sample-id", all =  TRUE)
alphadiversity <- left_join(alphadiversity,metadata,by=c("sample-id"="SampleID"))
my_comparisons <- list(c("HS Non-Legional Skin","HS Legional Skin"),
                       c("HS Non-Legional Skin","Healthy"),
                       c("HS Legional Skin","Healthy"))
shannon<-alphadiversity%>%
  ggplot(aes(x=Class, y=`shannon_entropy`, fill=Class)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  theme_q2r() +
  labs(y="Shannon Diversity")+
  stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        legend.position = "none",
        axis.title.y = element_text(size=12,face="bold"),
        axis.title.x = element_blank())
chao1<-alphadiversity%>%
  ggplot(aes(x=Class, y=`chao1`, fill=Class)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  theme_q2r() +
  labs(y="Chao1")+
  stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        legend.position = "none",
        axis.title.y = element_text(size=12,face="bold"),
        axis.title.x = element_blank())
pielou<-alphadiversity%>%
  ggplot(aes(x=Class, y=`pielou_evenness`, fill=Class)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  theme_q2r() +
  labs(y="Pielou Evenness")+
  stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        legend.position = "none",
        axis.title.y = element_text(size=12,face="bold"),
        axis.title.x = element_blank())
faithpd<-alphadiversity%>%
  ggplot(aes(x=Class, y=`faith_pd`, fill=Class)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  theme_q2r() +
  labs(y="Faith Phylogenetic Diversity")+
  stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        legend.position = "none",
        axis.title.y = element_text(size=12,face="bold"),
        axis.title.x = element_blank())
(shannon|faithpd)/(chao1|pielou)+
  plot_annotation(tag_levels = "A")&
  theme(plot.tag= element_text(size=16,face="bold"))
2. Beta diversity
metadata<-read_q2metadata("sample-metadata.tsv")
tax<-read_csv("ResultPlot/Taxonomy/level-7.csv") #Download from taxonomy_bar_plot.qzv
tax<-tax%>%
  pivot_longer(!c(index,`Sample_Name`,group,Class),names_to = "Taxonomy",values_to = "count")%>%select(index,Class,Taxonomy,count)%>%
  separate(Taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus","species"),sep=";")%>%
  mutate(genus=str_replace(string=genus,pattern="g__",replacement=""),
         phylum=str_replace(string=phylum,pattern="p__",replacement=""),
         family=str_replace(string=family,pattern="f__",replacement=""),
         species=str_replace(string=species,pattern="s__",replacement=""))
species<-tax%>%
  select(index,Class,species,count)%>%
  group_by(index,Class,species)%>%
  summarize(count=sum(count),.groups="drop")%>%
  group_by(index,Class)%>%
  mutate(numOtus=length(species),
         OTU = paste0("OTU", row_number()))
species<-species%>%
  group_by(OTU)%>%
  mutate(total=sum(count))%>%
  filter(total !=0)%>%ungroup()%>%select(index,count,OTU)%>%
  pivot_wider(names_from=OTU,values_from=count)%>%
  as.data.frame()
rownames(species)<-species$index
species <- species[,-1]
species<-as.matrix(species)
bray_dist<-vegdist(species,method="bray")
z<-as.matrix(bray_dist)
z<-z%>%as_tibble(rownames="sampleid")
meta_distance<- inner_join(z,metadata,by=c("sampleid"="SampleID"))
#Pairwise permanova
HSvsHSnonlegional<-meta_distance%>%
  filter(Class != "Healthy")
HSvsHSnonlegional_dist<-HSvsHSnonlegional%>%
  select(all_of(.[["sampleid"]]))%>%
  as.dist()
adonis2(HSvsHSnonlegional_dist~Class,data=HSvsHSnonlegional,permutation = 1000)
HSvsHealthy<-meta_distance%>%
  filter(Class != "HS Non-Legional Skin")
HSvsHealthy_dist<-HSvsHealthy%>%
  select(all_of(.[["sampleid"]]))%>%
  as.dist()
adonis2(HSvsHealthy_dist~Class,data=HSvsHealthy,permutation = 1000)
HSnonlegionalvsHealthy<-meta_distance%>%
  filter(Class != "HS Legional Skin")
HSnonlegionalvsHealthy_dist<-HSnonlegionalvsHealthy%>%
  select(all_of(.[["sampleid"]]))%>%
  as.dist()
adonis2(HSnonlegionalvsHealthy_dist~Class,data=HSnonlegionalvsHealthy,permutation = 1000)
HSvsHealthy<-meta_distance
HSvsHealthy_dist<-HSvsHealthy%>%
  select(all_of(.[["sampleid"]]))%>%
  as.dist()
adonis2(HSvsHealthy_dist~group,data=HSvsHealthy,permutation = 1000)
#######################
write_tsv(z,"ResultPlot/Alpha Beta/bray_dist.tsv")
pcoa<-cmdscale(bray_dist, k = 2, eig=TRUE,add=TRUE)
coordination<-pcoa$points
colnames(coordination)<-c("pcoa1","pcoa2")
percentage_explained<-round((100 * pcoa$eig / sum(pcoa$eig)),2)
####
labs<-c(glue("PC1 ({percentage_explained[1]}%)"),
        glue("PC2 ({percentage_explained[2]}%)"))
#HSvsHSnonlegionalvsHealthy
coordination %>%
  as_tibble(rownames="sampleid")%>%
  inner_join(.,metadata,by=c("sampleid"="SampleID"))%>%
  ggplot(aes(x=pcoa1,y=pcoa2,color=Class,fill=Class))+
  stat_ellipse(geom="polygon",alpha=0.5)+
  geom_point(shape=21,col="black")+
  scale_color_manual(name=NULL,
                     breaks=c("HS Non-Legional Skin","HS Legional Skin",
                              "Healthy"),
                     values=c("grey75","#ffb09c","#ee2400"))+
  scale_fill_manual(name=NULL,
                    breaks=c("HS Non-Legional Skin","HS Legional Skin",
                             "Healthy"),
                    values=c("grey75","#ffb09c","#ee2400"))+
  labs(x=labs[1],y=labs[2])+
  theme_q2r()+
  theme(axis.text.x = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=12,face="bold"),
        legend.position = "bottom")+
  guides(color = "none",
         fill = guide_legend(override.aes = list(color = c("grey75", "#ffb09c", "#ee2400")),
                             order = 1))
#HSvsHealthy
coordination %>%
  as_tibble(rownames="sampleid")%>%
  inner_join(.,metadata,by=c("sampleid"="SampleID"))%>%
  ggplot(aes(x=pcoa1,y=pcoa2,color=group,fill=group))+
  stat_ellipse(geom="polygon",alpha=0.5)+
  geom_point(shape=21,col="black")+
  scale_color_manual(name=NULL,
                     breaks=c("HS","Healthy"),
                     values=c("grey75","#ee2400"))+
  scale_fill_manual(name=NULL,
                    breaks=c("HS","Healthy"),
                    values=c("grey75","#ee2400"))+
  labs(x=labs[1],y=labs[2])+
  theme_q2r()+
  theme(axis.text.x = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=12,face="bold"),
        legend.position = "bottom")+
  guides(color = "none",
         fill = guide_legend(override.aes = list(color = c("grey75", "#ee2400")),
                             order = 1))
3. Taxonomy barplot
########Species##############
metadata<-read_q2metadata("sample-metadata.tsv")
tax<-read_csv("ResultPlot/Taxonomy/level-7.csv")
tax<-tax%>%
  pivot_longer(!c(index,`Sample_Name`,group,Class),names_to = "Taxonomy",values_to = "count")%>%select(index,Class,Taxonomy,count)%>%
  separate(Taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus","species"),sep=";")%>%
  mutate(genus=str_replace(string=genus,pattern="g__",replacement=""),
         phylum=str_replace(string=phylum,pattern="p__",replacement=""),
         family=str_replace(string=family,pattern="f__",replacement=""),
         species=str_replace(string=species,pattern="s__",replacement=""))
species<-tax%>%
  select(index,Class,genus,species,count)%>%
  mutate(species=ifelse(species=="__",paste0(genus,"_Unclassifier"),species))
values <- species$species
matches <- grep("Unclassifier", values)
unclassified_values <- values[matches]
species <- species %>%
  filter(!species %in% unclassified_values)%>%
  group_by(index)%>%
  mutate(abundance=count/sum(count))%>%
  group_by(index,Class,species)%>%
  summarize(abundance=sum(abundance),.groups="drop")%>%
  group_by(Class,species)%>%
  summarize(mean_abundance=100*mean(abundance),.groups="drop")
a<-species%>%
  group_by(species)%>%
  summarize(max=max(mean_abundance))%>%
  arrange(desc(max))
Taxon_pool<-species%>%
  group_by(species)%>%
  summarize(pool=max(mean_abundance)<1.56,.groups="drop")
species%>%
  inner_join(.,Taxon_pool,by="species")%>%
  mutate(species = if_else(pool,"Others",species)) %>%
  group_by(Class,species)%>%
  summarize(mean_abundance=sum(mean_abundance))%>%
  mutate(species=fct_reorder(species,mean_abundance))%>%
  mutate(species = fct_relevel(species, "Others", after = Inf)) %>%
  ggplot(aes(x=Class,y=mean_abundance,fill=species))+
  geom_col()+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(name = NULL, values = c(
    brewer.pal(8, "Dark2"), # 8 colors from Dark2 palette
    "blue", "gray30", "yellow", "red", "darkmagenta", "green", "pink2", "darkgreen",
    "cyan", "orange", "purple", "lightblue", "lightgreen", "salmon", "gold", "darkred", "navy", "orchid"
  )) +
  xlab("Class")+
  ylab("Relative Abundance (%)")+
  theme_q2r()+
  theme(axis.text.x=element_markdown(size=14,angle=45,hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=14),
        axis.text.y=element_markdown(size=14),
        legend.text = element_markdown(size=14,face="italic"))+
  guides(fill = guide_legend(ncol = 1))  
########Genus##############
metadata<-read_q2metadata("sample-metadata.tsv")
tax<-read_csv("ResultPlot/Taxonomy/level-6.csv")
tax<-tax%>%
  pivot_longer(!c(index,`Sample_Name`,group,Class),names_to = "Taxonomy",values_to = "count")%>%select(index,Class,Taxonomy,count)%>%
  separate(Taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),sep=";")%>%
  mutate(genus=str_replace(string=genus,pattern="g__",replacement=""),
         phylum=str_replace(string=phylum,pattern="p__",replacement=""),
         family=str_replace(string=family,pattern="f__",replacement=""))
genus<-tax%>%
  select(index,Class,genus,count)%>%
  group_by(index)%>%
  mutate(abundance=count/sum(count))%>%
  group_by(index,Class,genus)%>%
  summarize(abundance=sum(abundance),.groups="drop")%>%
  group_by(Class,genus)%>%
  summarize(mean_abundance=100*mean(abundance),.groups="drop")
a<-genus%>%
  group_by(genus)%>%
  summarize(max=max(mean_abundance))%>%
  arrange(desc(max))
Taxon_pool<-genus%>%
  group_by(genus)%>%
  summarize(pool=max(mean_abundance)<1.56,.groups="drop")
genus%>%
  inner_join(.,Taxon_pool,by="genus")%>%
  mutate(genus = if_else(pool,"Others",genus)) %>%
  group_by(Class,genus)%>%vvvvvvvvvvvvvv
  summarize(mean_abundance=sum(mean_abundance))%>%
  mutate(genus=fct_reorder(genus,mean_abundance))%>%
  mutate(genus = fct_relevel(genus, "Others", after = Inf)) %>%
  ggplot(aes(x=Class,y=mean_abundance,fill=genus))+
  geom_col()+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(name = NULL, values = c(
    brewer.pal(8, "Dark2"), # 8 colors from Dark2 palette
    "blue", "gray30", "yellow", "red", "darkmagenta", "green", "pink2", "darkgreen",
    "cyan", "orange", "purple", "lightblue", "lightgreen", "salmon", "gold", "darkred", "navy", "orchid"
  )) +
  xlab("Class")+
  ylab("Relative Abundance (%)")+
  theme_q2r()+
  theme(axis.text.x=element_markdown(size=14,angle=45,hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=14),
        axis.text.y=element_markdown(size=14),
        legend.text = element_markdown(size=14,face="italic"))+
  guides(fill = guide_legend(ncol = 1))
4. Sig Kruskal Wallis
###HealthyvsHSnonleginal
tax<-read_csv("ResultPlot/Taxonomy/level-6.csv")
tax<-tax%>%
  pivot_longer(!c(index,`Sample_Name`,group,Class),names_to = "Taxonomy",values_to = "count")%>%select(index,Class,Taxonomy,count)%>%
  separate(Taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),sep=";")%>%
  mutate(genus=str_replace(string=genus,pattern="g__",replacement=""),
         phylum=str_replace(string=phylum,pattern="p__",replacement=""),
         family=str_replace(string=family,pattern="f__",replacement=""))
genus<-tax%>%
  select(index,Class,genus,count)%>%
  group_by(index)%>%
  mutate(abundance=100*count/sum(count))%>%
  group_by(index,Class,genus)%>%
  summarize(abundance=sum(abundance),.groups="drop")%>%
  mutate(genus = ifelse(genus=="__","Unclassifier",genus))
sig_Healthy_HSnonlegional<-genus%>%
  filter(Class=="Healthy"| Class=="HS Non-Legional Skin")%>%
  nest(data = -genus)%>%
  mutate(test = map(.x=data, ~kruskal.test(abundance~Class,data=.x) %>% tidy)) %>%
  unnest(test) %>%
  filter (p.value < 0.05) %>%
  select (genus,p.value)
genus<-genus%>%filter(Class=="Healthy"| Class=="HS Non-Legional Skin")%>%
  group_by(Class,genus)%>%
  summarize(mean_abundance=mean(abundance),.groups="drop")%>%
  inner_join(.,sig_Healthy_HSnonlegional,by="genus")%>%
  mutate(mean_abundance=if_else(mean_abundance==0.00000,0.00001,mean_abundance))
genus_Healthy<-genus%>%filter(Class=="Healthy")
genus_HSnonlegional<-genus%>%filter(Class=="HS Non-Legional Skin")
LFC<-data.frame(genus=genus_Healthy$genus,
                LFC=log2(genus_Healthy$mean_abundance/genus_HSnonlegional$mean_abundance))
LFC%>%
  mutate(Class=if_else(LFC>0,"Healthy","HS Non-Legional Skin"),genus=fct_reorder(genus,LFC))%>%
  mutate(label_x=if_else(Class=="HS Non-Legional Skin",0.1,-0.1),
         label_hjust=if_else(Class=="HS Non-Legional Skin",0,1))%>%
  ggplot(aes(x=LFC,y=genus,label=genus,fill=Class))+geom_col()+
  geom_vline(xintercept =seq(-18,10,by=2),linetype="dotted",color="darkgray")+
  geom_richtext(aes(x = label_x, hjust = label_hjust, label = paste0("<i>", genus, "</i>")),
                fill = NA, label.color = NA, size = 5)+
  scale_fill_manual(name=NULL,
                    breaks=c("Healthy","HS Non-Legional Skin"),
                    values=c("blue", "grey"))+
  labs(x="Log Fold Change (LFC)",y="genus")+
  scale_x_continuous(breaks = c(-18, seq(-18, 10, by = 2)))+
  theme_classic()+
  theme(
    axis.text.y=element_blank(),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_markdown(size=14),
    legend.text=element_markdown(size=14),
    axis.title=element_markdown(size=14),
    legend.position = "bottom"
  )
###HealthyvsHSLegional
tax<-read_csv("ResultPlot/Taxonomy/level-6.csv")
tax<-tax%>%
  pivot_longer(!c(index,`Sample_Name`,group,Class),names_to = "Taxonomy",values_to = "count")%>%select(index,Class,Taxonomy,count)%>%
  separate(Taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),sep=";")%>%
  mutate(genus=str_replace(string=genus,pattern="g__",replacement=""),
         phylum=str_replace(string=phylum,pattern="p__",replacement=""),
         family=str_replace(string=family,pattern="f__",replacement=""))
genus<-tax%>%
  select(index,Class,genus,count)%>%
  group_by(index)%>%
  mutate(abundance=100*count/sum(count))%>%
  group_by(index,Class,genus)%>%
  summarize(abundance=sum(abundance),.groups="drop")%>%
  mutate(genus = ifelse(genus=="__","Unclassifier",genus))
sig_Healthy_HSlegional<-genus%>%
  filter(Class=="Healthy"| Class=="HS Legional Skin")%>%
  nest(data = -genus)%>%
  mutate(test = map(.x=data, ~kruskal.test(abundance~Class,data=.x) %>% tidy)) %>%
  unnest(test) %>%
  filter (p.value < 0.01) %>%
  select (genus,p.value)
genus<-genus%>%filter(Class=="Healthy"| Class=="HS Legional Skin")%>%
  group_by(Class,genus)%>%
  summarize(mean_abundance=mean(abundance),.groups="drop")%>%
  inner_join(.,sig_Healthy_HSlegional,by="genus")%>%
  mutate(mean_abundance=if_else(mean_abundance==0.00000,0.00001,mean_abundance))
genus_Healthy<-genus%>%filter(Class=="Healthy")
genus_HSlegional<-genus%>%filter(Class=="HS Legional Skin")
LFC<-data.frame(genus=genus_Healthy$genus,
                LFC=log2(genus_Healthy$mean_abundance/genus_HSlegional$mean_abundance))
LFC%>%
  mutate(Class=if_else(LFC>0,"Healthy","HS Legional Skin"),genus=fct_reorder(genus,LFC))%>%
  mutate(label_x=if_else(Class=="HS Legional Skin",0.1,-0.1),
         label_hjust=if_else(Class=="HS Legional Skin",0,1))%>%
  ggplot(aes(x=LFC,y=genus,label=genus,fill=Class))+geom_col()+
  geom_vline(xintercept =seq(-18,4,by=2),linetype="dotted",color="darkgray")+
  geom_richtext(aes(x = label_x, hjust = label_hjust, label = paste0("<i>", genus, "</i>")),
                fill = NA, label.color = NA, size = 5)+
  scale_fill_manual(name=NULL,
                    breaks=c("Healthy","HS Legional Skin"),
                    values=c("blue", "grey"))+
  labs(x="Log Fold Change (LFC)",y="Genus")+
  scale_x_continuous(breaks = c(seq(-18, 4, by = 2),10),
                     labels = function(x) ifelse(x == 10, "**_p.value_**", x))+
  theme_classic()+
  theme(
    axis.text.y=element_blank(),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_markdown(size=14),
    legend.text=element_markdown(size=14),
    axis.title=element_markdown(size=14),
    legend.position = "bottom"
  )+
  geom_text(data=tibble(x=10,y=1),
            aes(x=x,y=y, label="p = 2.04e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=2),
            aes(x=x,y=y, label="p = 1.09e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=3),
            aes(x=x,y=y, label="p = 2.05e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=4),
            aes(x=x,y=y, label="p = 2.05e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=5),
            aes(x=x,y=y, label="p = 2.05e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=6),
            aes(x=x,y=y, label="p = 6.85e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=7),
            aes(x=x,y=y, label="p = 6.85e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=8),
            aes(x=x,y=y, label="p = 1.72e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=9),
            aes(x=x,y=y, label="p = 9.23e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=10),
            aes(x=x,y=y, label="p = 2.34e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=11),
            aes(x=x,y=y, label="p = 2.36e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=12),
            aes(x=x,y=y, label="p = 1.46e-05"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=13),
            aes(x=x,y=y, label="p = 1.26e-04"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=14),
            aes(x=x,y=y, label="p = 4.26e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=15),
            aes(x=x,y=y, label="p = 9.33e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=16),
            aes(x=x,y=y, label="p = 4.86e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=17),
            aes(x=x,y=y, label="p = 9.99e-04"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)
###HSNonLegionalvsHSLegional
tax<-read_csv("ResultPlot/Taxonomy/level-6.csv")
tax<-tax%>%
  pivot_longer(!c(index,`Sample_Name`,group,Class),names_to = "Taxonomy",values_to = "count")%>%select(index,Class,Taxonomy,count)%>%
  separate(Taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),sep=";")%>%
  mutate(genus=str_replace(string=genus,pattern="g__",replacement=""),
         phylum=str_replace(string=phylum,pattern="p__",replacement=""),
         family=str_replace(string=family,pattern="f__",replacement=""))
genus<-tax%>%
  select(index,Class,genus,count)%>%
  group_by(index)%>%
  mutate(abundance=100*count/sum(count))%>%
  group_by(index,Class,genus)%>%
  summarize(abundance=sum(abundance),.groups="drop")%>%
  mutate(genus = ifelse(genus=="__","Unclassifier",genus))
sig_HSnonlegional_HSlegional<-genus%>%
  filter(Class=="HS Non-Legional Skin"| Class=="HS Legional Skin")%>%
  nest(data = -genus)%>%
  mutate(test = map(.x=data, ~kruskal.test(abundance~Class,data=.x) %>% tidy)) %>%
  unnest(test) %>%
  filter (p.value < 0.01) %>%
  select (genus,p.value)
genus<-genus%>%filter(Class=="HS Non-Legional Skin"| Class=="HS Legional Skin")%>%
  group_by(Class,genus)%>%
  summarize(mean_abundance=mean(abundance),.groups="drop")%>%
  inner_join(.,sig_HSnonlegional_HSlegional,by="genus")%>%
  mutate(mean_abundance=if_else(mean_abundance==0.00000,0.00001,mean_abundance))
genus_HSnonlegional<-genus%>%filter(Class=="HS Non-Legional Skin")
genus_HSlegional<-genus%>%filter(Class=="HS Legional Skin")
LFC<-data.frame(genus=genus_HSnonlegional$genus,
                LFC=log2(genus_HSnonlegional$mean_abundance/genus_HSlegional$mean_abundance))
LFC%>%
  mutate(Class=if_else(LFC>0,"HS Non-Legional Skin","HS Legional Skin"),genus=fct_reorder(genus,LFC))%>%
  mutate(label_x=if_else(Class=="HS Legional Skin",0.1,-0.1),
         label_hjust=if_else(Class=="HS Legional Skin",0,1))%>%
  ggplot(aes(x=LFC,y=genus,label=genus,fill=Class))+geom_col()+
  geom_vline(xintercept =seq(-16,4,by=2),linetype="dotted",color="darkgray")+
  geom_richtext(aes(x = label_x, hjust = label_hjust, label = paste0("<i>", genus, "</i>")),
                fill = NA, label.color = NA, size = 5)+
  scale_fill_manual(name=NULL,
                    breaks=c("HS Non-Legional Skin","HS Legional Skin"),
                    values=c("blue", "grey"))+
  labs(x="Log Fold Change (LFC)",y="genus")+
  scale_x_continuous(breaks = c(seq(-16, 4, by = 2),10),
                     labels = function(x) ifelse(x == 10, "**_p.value_**", x))+
  theme_classic()+
  theme(
    axis.text.y=element_blank(),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_markdown(size=14),
    legend.text=element_markdown(size=14),
    axis.title=element_markdown(size=14),
    legend.position = "bottom"
  )+
  geom_text(data=tibble(x=10,y=1),
            aes(x=x,y=y, label="p = 3.12e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=2),
            aes(x=x,y=y, label="p = 8.69e-05"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=3),
            aes(x=x,y=y, label="p = 2.78e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=4),
            aes(x=x,y=y, label="p = 7.91e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=5),
            aes(x=x,y=y, label="p = 7.06e-04"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=6),
            aes(x=x,y=y, label="p = 1.21e-04"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=7),
            aes(x=x,y=y, label="p = 1.09e-04"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=8),
            aes(x=x,y=y, label="p = 1.49e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=9),
            aes(x=x,y=y, label="p = 9.02e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=10,y=10),
            aes(x=x,y=y, label="p = 9.72e-03"), size = 4, fontface="bold.italic",
            inherit.aes=FALSE)
###species###HealthyvsHSLegional
tax<-read_csv("ResultPlot/Taxonomy/level-7.csv")
tax<-tax%>%
  pivot_longer(!c(index,`Sample_Name`,group,Class),names_to = "Taxonomy",values_to = "count")%>%select(index,Class,Taxonomy,count)%>%
  separate(Taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus","species"),sep=";")%>%
  mutate(genus=str_replace(string=genus,pattern="g__",replacement=""),
         phylum=str_replace(string=phylum,pattern="p__",replacement=""),
         family=str_replace(string=family,pattern="f__",replacement=""),
         species=str_replace(string=species,pattern="s__",replacement=""))
species<-tax%>%
  select(index,Class,genus,species,count)%>%
  filter(genus!="__")%>%
  mutate(species = ifelse(species=="__",paste(genus,"_Unclassifier"),species))%>%
  group_by(index)%>%
  mutate(abundance=100*count/sum(count))%>%
  group_by(index,Class,species)%>%
  summarize(abundance=sum(abundance),.groups="drop")
sig_Healthy_HSlegional<-species%>%
  filter(Class=="Healthy"| Class=="HS Legional Skin")%>%
  nest(data = -species)%>%
  mutate(test = map(.x=data, ~kruskal.test(abundance~Class,data=.x) %>% tidy)) %>%
  unnest(test) %>%
  filter (p.value < 0.01) %>%
  select (species,p.value)
species<-species%>%filter(Class=="Healthy"| Class=="HS Legional Skin")%>%
  group_by(Class,species)%>%
  summarize(mean_abundance=mean(abundance),.groups="drop")%>%
  inner_join(.,sig_Healthy_HSlegional,by="species")%>%
  mutate(mean_abundance=if_else(mean_abundance==0.00000,0.00001,mean_abundance))
species_Healthy<-species%>%filter(Class=="Healthy")
species_HSlegional<-species%>%filter(Class=="HS Legional Skin")
LFC<-data.frame(species=species_Healthy$species,
                LFC=log2(species_Healthy$mean_abundance/species_HSlegional$mean_abundance))
LFC%>%
  mutate(Class=if_else(LFC>0,"Healthy","HS Legional Skin"),species=fct_reorder(species,LFC))%>%
  mutate(label_x=if_else(Class=="HS Legional Skin",0.1,-0.1),
         label_hjust=if_else(Class=="HS Legional Skin",0,1))%>%
  ggplot(aes(x=LFC,y=species,label=species,fill=Class))+geom_col()+
  geom_vline(xintercept =seq(-18,18,by=2),linetype="dotted",color="darkgray")+
  geom_richtext(aes(x = label_x, hjust = label_hjust, label = paste0("<i>", species, "</i>")),
                fill = NA, label.color = NA, size = 5)+
  scale_fill_manual(name=NULL,
                    breaks=c("Healthy","HS Legional Skin"),
                    values=c("blue", "grey"))+
  labs(x="Log Fold Change (LFC)",y="Species")+
  scale_x_continuous(breaks = c(seq(-18, 18, by = 2)))+
  theme_classic()+
  theme(
    axis.text.y=element_blank(),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_markdown(size=14),
    legend.text=element_markdown(size=14),
    axis.title=element_markdown(size=14),
    legend.position = "bottom"
  )
###species###HSvsHSnon
tax<-read_csv("ResultPlot/Taxonomy/level-7.csv")
tax<-tax%>%
  pivot_longer(!c(index,`Sample_Name`,group,Class),names_to = "Taxonomy",values_to = "count")%>%select(index,Class,Taxonomy,count)%>%
  separate(Taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus","species"),sep=";")%>%
  mutate(genus=str_replace(string=genus,pattern="g__",replacement=""),
         phylum=str_replace(string=phylum,pattern="p__",replacement=""),
         family=str_replace(string=family,pattern="f__",replacement=""),
         species=str_replace(string=species,pattern="s__",replacement=""))
species<-tax%>%
  select(index,Class,genus,species,count)%>%
  filter(genus!="__")%>%
  mutate(species = ifelse(species=="__",paste(genus,"_Unclassifier"),species))%>%
  group_by(index)%>%
  mutate(abundance=100*count/sum(count))%>%
  group_by(index,Class,species)%>%
  summarize(abundance=sum(abundance),.groups="drop")
sig_HSnon_HSlegional<-species%>%
  filter(Class=="HS Non-Legional Skin"| Class=="HS Legional Skin")%>%
  nest(data = -species)%>%
  mutate(test = map(.x=data, ~kruskal.test(abundance~Class,data=.x) %>% tidy)) %>%
  unnest(test) %>%
  filter (p.value < 0.05) %>%
  select (species,p.value)
species<-species%>%filter(Class=="HS Non-Legional Skin" | Class=="HS Legional Skin")%>%
  group_by(Class,species)%>%
  summarize(mean_abundance=mean(abundance),.groups="drop")%>%
  inner_join(.,sig_HSnon_HSlegional,by="species")%>%
  mutate(mean_abundance=if_else(mean_abundance==0.00000,0.00001,mean_abundance))
species_HSnonlegional<-species%>%filter(Class=="HS Non-Legional Skin")
species_HSlegional<-species%>%filter(Class=="HS Legional Skin")
LFC<-data.frame(species=species_HSnonlegional$species,
                LFC=log2(species_HSnonlegional$mean_abundance/species_HSlegional$mean_abundance))
LFC%>%
  mutate(Class=if_else(LFC>0,"HS Non-Legional Skin","HS Legional Skin"),species=fct_reorder(species,LFC))%>%
  mutate(label_x=if_else(Class=="HS Legional Skin",0.1,-0.1),
         label_hjust=if_else(Class=="HS Legional Skin",0,1))%>%
  ggplot(aes(x=LFC,y=species,label=species,fill=Class))+geom_col()+
  geom_vline(xintercept =seq(-18,18,by=2),linetype="dotted",color="darkgray")+
  geom_richtext(aes(x = label_x, hjust = label_hjust, label = paste0("<i>", species, "</i>")),
                fill = NA, label.color = NA, size = 5)+
  scale_fill_manual(name=NULL,
                    breaks=c("HS Non-Legional Skin","HS Legional Skin"),
                    values=c("blue", "grey"))+
  labs(x="Log Fold Change (LFC)",y="Species")+
  scale_x_continuous(breaks = c(seq(-18, 18, by = 2)))+
  theme_classic()+
  theme(
    axis.text.y=element_blank(),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_markdown(size=14),
    legend.text=element_markdown(size=14),
    axis.title=element_markdown(size=14),
    legend.position = "bottom"
  )
