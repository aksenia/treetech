#### Script documenting the main data tables in the paper
#### using supplementary data 
#### code for the main figures is also included 
#### Author: KL, Ksenia, Lavrichenko

### set the working directory
### the structure needs to be mimicked as in this code
### and the Excel sheets are split and saved as txt files
setwd("~/Documents/GiB/data/cnvr/")
### root/
### --merge/ Folder with files from Supplementart file 2
### --seg_filter/ Folder with list of segments to filter agains (bed file)
### --vapor_out/ Folder with results of the VaPoR runs - Supplementart file 3
### --array_support Folder wih array support table - Supplementart file 4

## package manager
library(pacman)
## work packages
p_load(dplyr, tidyr, tools, scales, stringr, bedr)
## viz packages
p_load(ggplot2, scales, ggpubr, ggsci, ggrepel, gridExtra)

### final sets cnvr directory
## merged 
## files in Supplementary file 2
cnvrfiles <- list.files(file.path("merge"), pattern = "txt$", full.names = T)
names(cnvrfiles) <- basename(file_path_sans_ext(cnvrfiles))

## read and merge into one table
mcnvr <- do.call("rbind", lapply(names(cnvrfiles), function(f)
  read.table(cnvrfiles[[f]], 
             stringsAsFactors = F,
             sep="\t",
             col.names = c("chr", "start", "end", "count", "dh.score", 
                           "i_chr", "i_start", "i_end", "i_id", "i_dh.score", "overlap")) %>%
    mutate(LR_bin=f))) %>%
  separate(LR_bin, c("dummy1", "lr_bin", "dummy2"), sep="\\.") %>%
  dplyr::select(-matches("dummy")) %>%
  mutate(cnvrid=paste(chr, ":", start, "-", end, sep=""))

### FILTER OUT the BLACKLIST REGIONS
### blacklist 
blacklist <- read.table("seg_filter/centro_telo_ECD+blacklist.bed",
                        col.names = c("chr", "start", "end")) %>%
  mutate(coordseg = paste0(chr, ":", start, "-", end))
blist.sort <- bedr.sort.region(blacklist$coordseg)
## use bedr package to subtract blacklisted regions
mcnvr.sort <- bedr.sort.region(unique(mcnvr$cnvrid))
mcnvr.sub <- bedr.subtract.region(mcnvr.sort, blist.sort)

### VaPoR files (Supplementary file 3)
vaporpol <- list.files("vapor_out", full.names = T)
names(vaporpol) <- gsub("\\..*","", gsub("cnvr\\.","", basename(file_path_sans_ext(vaporpol))))
## read in and merge
vpor <- do.call("rbind", lapply(names(vaporpol), function(fn)
  read.table(vaporpol[[fn]],
             stringsAsFactors = F,
             header = T) %>%
    mutate(tech=fn))) %>%
  dplyr::rename(i_chr=CHR,
                i_start=POS,
                i_end=END,
                type=SVTYPE) %>%
  mutate(type=gsub("TAN", "", type)) %>%
  distinct() %>%
  rowwise() %>% 
  mutate(Nread=str_count(VaPoR_Rec, ','),
         VaPoR_mQS=median(as.numeric(strsplit(VaPoR_Rec, split=",")[[1]])))

### final table
mcnvR <- mcnvr %>%
  filter(cnvrid %in% mcnvr.sub) %>%
  mutate(lr_bin=factor(lr_bin, levels = c("all", "above1", "above5", "above1_HQ")),
    len=end-start,
         i_len=i_end-i_start,
         i_lenkb=i_len/1000,
         perc=overlap/len*100,
         len_bin=cut(i_lenkb, breaks = c(0, 1, 5, 20, 100, 500, 1000, 2000), labels = c("<1kb", "1-5kb", "5-20kb", "20-100kb", "100-500kb", "500kb-1Mb", ">1Mb"))) %>%
  separate(i_id, c("qbin", "type", "tech", "supp","coordcnv", "i_count"), sep="\\.", remove=F)  %>%
  left_join(vpor)

#### summarize per CNV locus wih many CNVRs
mcnvRcol <- mcnvR %>%
  dplyr::group_by(cnvrid, chr, len, type, qbin, tech, dh.score, lr_bin) %>%
  dplyr::summarise(totlen=sum(i_len), 
            med_dh.score=median(i_dh.score),
            med_vpor.score=median(VaPoR_GS),
            med_vpor.nread=median(Nread),
            nseg = n_distinct(i_id),
            lGT=paste(sort(unique(VaPoR_GT)), collapse = ",")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(lr_bin, cnvrid, tech) %>%
  dplyr::mutate(tbins=paste(sort(unique(qbin)), collapse ="|"),
         totperc = totlen/len*100) 

### If both HQ and LQ segments are included per technology, 
### select one, prefer the HQ one in hat case
lmcnvRcol <- with(mcnvRcol, split(mcnvRcol, tbins))
## qbins is summarised within each tech, meaning 
## all of them have both LQ and HQ
## can just select one 
mixed_qbins <- lmcnvRcol[["HQ|LQ"]] %>%
  filter(qbin=="HQ")
## first desambiguate within tech
## then summarise across tech
McnvRcol <- do.call("rbind", list(lmcnvRcol[["HQ"]], lmcnvRcol[["LQ"]], mixed_qbins)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(lr_bin, cnvrid) %>%
  dplyr::mutate(ntech=n_distinct(tech),
         ltech=paste(sort(unique(tech)), collapse = "|"),
         qbins=paste(sort(unique(qbin)), collapse ="|"),
         lenkb=len/1000,
         len_bin=cut(lenkb, breaks = c(0, 1, 5, 20, 100, 500, 2000), labels = c("<1", "1-5", "5-20", "20-100", "100-500", ">500")))

##### FIGURES 
#### FIGURE 2 
## filter by CNV locus length - uniformly across all data
pMcnvRcol <- McnvRcol %>%
  dplyr::filter(len >= 500) %>%
  dplyr::group_by(len_bin, tech, ltech, lr_bin) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(len_bin, tech, lr_bin) %>%
  dplyr::mutate(n_tot = sum(n),
                perc = n/n_tot*100,
                ltech = ifelse(ltech=="array"|ltech=="sr"|ltech=="lr", "private", ltech),
                ltech=factor(ltech, levels = c("array|lr|sr", "array|lr", "array|sr", "lr|sr",
                                               "private")))

Fig2b <-ggplot(pMcnvRcol %>%
           filter(lr_bin=="above1") %>%
             mutate(len_bin=gsub("-", "-\n", len_bin),
                    len_bin=gsub("1-\n5", "1-5", len_bin),
                    len_bin=gsub("5-\n20", "5-20", len_bin)) %>%
             mutate(len_bin=factor(len_bin, levels =  c("<1", "1-5", "5-20", "20-\n100", "100-\n500", ">500"))), 
         aes(x=len_bin, y=n, fill=ltech)) +
  geom_bar(stat = "identity", color="white", linetype="dotted", alpha=0.9) +
  geom_text(aes(x=len_bin, y=n, label=ifelse(n>50,n,"")), vjust=1.2, position = "stack", color="white") +
  theme_pubclean() +
  scale_fill_npg() +
  theme(axis.text.x = element_text(size=12),
        legend.position = "none",
        plot.margin = unit(c(0,3,1,3), "lines"),
        text=element_text(size=17)) +
  ylab("CNVR count") + xlab("CNV size bins in Kb") +
  facet_wrap(~tech, nrow=1)

## group two similar bins 1-5 and 5-20 
## they have very similar proportions profiles
psMcnvRcol <- McnvRcol %>%
  dplyr::filter(len >= 500) %>%
  rowwise() %>%
  dplyr::mutate(len_bin=gsub("1-5", "1-20", len_bin),
                len_bin=gsub("5-20", "1-20", len_bin)) %>%
  dplyr::group_by(len_bin, tech, ltech, lr_bin) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(len_bin, tech, lr_bin) %>%
  dplyr::mutate(n_tot = sum(n),
                perc = n/n_tot*100,
                ltech = ifelse(ltech=="array"|ltech=="sr"|ltech=="lr", "private", ltech),
                ltech=gsub("\\|", "+", ltech),
                ltech=factor(ltech, levels = c("array+lr+sr", "array+lr", "array+sr", 
                                               "lr+sr", "private")))

Fig2e <-ggplot(psMcnvRcol %>%
                 filter(lr_bin=="above1") %>%
                 mutate(len_bin=gsub("-", "-\n", len_bin),
                        len_bin=gsub("1-\n20", "1-20", len_bin)) %>%
                 mutate(len_bin=factor(len_bin, levels =  c("<1", "1-20", "20-\n100", "100-\n500", ">500"))), 
               aes(x=len_bin, y=perc, fill=ltech)) +
  geom_bar(stat = "identity", color="white", linetype="dotted", alpha=0.9) +
  coord_polar() +
  scale_fill_npg() +
  theme_pubclean() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,1,0,1), "lines"),
        panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(-0.3, "lines"),
        legend.position = "right",
        text=element_text(size=17)) +
  labs(fill='CNVR supported by') +
  facet_wrap(~tech, nrow=1)

### SUPPL FIG S2 Composition and percentages for CNV loci, per long-read quality bin
FigS2 <-  ggplot(psMcnvRcol %>%
                    mutate(len_bin=gsub("-", "-\n", len_bin),
                           len_bin=gsub("1-\n20", "1-20", len_bin)) %>%
                    mutate(len_bin=factor(len_bin, levels =  c("<1", "1-20", "20-\n100", "100-\n500", ">500"))), 
                  aes(x=len_bin, y=perc, fill=ltech)) +
  geom_bar(stat = "identity", linetype="dotted", alpha=0.8) +
  scale_fill_npg() +
  coord_polar() +
  theme_pubclean() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=6),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,2,0,2), "lines"),
        panel.spacing.x = unit(0.5, "lines"),
        panel.spacing.y = unit(0.4, "lines"),
        legend.position = "right",
        text=element_text(size=16)) +
  labs(fill="CNVR supported by") +
  facet_wrap(~tech + lr_bin, nrow=3)


##########  SCORE densities
#### FIGURE 3 
ab1 <- McnvRcol %>% filter(med_dh.score<5, lr_bin=="above1")
y <- McnvRcol %>% 
  dplyr::filter(med_dh.score<5 || (type=="DUP" && ntech==3), len>=500) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(lbin = ifelse(tech=="lr", as.character(lr_bin), tech)) %>%
  dplyr::select(-tech) %>%
  dplyr::distinct()

b1 <- y %>% 
  dplyr::select(-lr_bin) %>%
  dplyr::mutate(lr_bin=ifelse(lbin %in% c("all","above5","above1_HQ"), "lr_bins", "all_tech"),
         lbin=factor(lbin, levels=c("array", "all", "above1", "above5", "above1_HQ", "lr", "sr"))) %>%
  dplyr::mutate(dh_bin=cut(med_dh.score, breaks = seq(-0.1,5, by=0.1), labels = seq(0,5, by=0.1))) %>%
  dplyr::group_by(lbin, type, lr_bin) %>%
  dplyr::mutate(n_total=n()) %>%
  dplyr::group_by(lbin, dh_bin, type, lr_bin, n_total) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(perc=n/n_total*100,
         dh_bin = as.numeric(as.character(dh_bin)))

b2 <- y %>% 
  dplyr::filter(grepl("above1", lr_bin)) %>%
  dplyr::mutate(lbin=factor(lbin, levels=c("array", "all", "above1", "above5", "above1_HQ", "lr", "sr"))) %>%
  dplyr::mutate(dh_bin=cut(med_dh.score, breaks = seq(-0.1,5, by=0.1), labels = seq(0,5, by=0.1))) %>%
  dplyr::group_by(lbin, type, ntech) %>%
  dplyr::mutate(n_total=n()) %>%
  dplyr::group_by(lbin, dh_bin, type, ntech, n_total) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(perc=n/n_total*100,
         dh_bin = as.numeric(as.character(dh_bin)))

## array score data, Supplementary file 4
siglrr <- read.table(file = "array_support/array_signalfile.txt",
            sep="\t",
            header = T) %>%
  dplyr::inner_join(mcnvR %>% filter(lr_bin=="above1")) %>%
  dplyr::left_join(McnvRcol %>% dplyr::ungroup() %>% 
                     dplyr::filter(lr_bin=="above1") %>% 
                     dplyr::select(cnvrid, ltech, qbins) %>% dplyr::distinct())


c1 <- siglrr %>% 
  dplyr::filter(len >=500) %>%
  dplyr::mutate(sc_bin=cut(hel.dist, breaks = seq(-0.05,5, by=0.05), labels = seq(0,5, by=0.05))) %>%
  dplyr::group_by(qbin, tech) %>%
  dplyr::mutate(n_total=n()) %>%
  dplyr::group_by(qbin, tech, sc_bin, n_total) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(perc=n/n_total*100,
                sc_bin = as.numeric(as.character(sc_bin)),
                tech = gsub("lr", "long reads", tech)) %>%
  dplyr::mutate(tech = gsub("sr", "short reads", tech)) %>%
  dplyr::rename(Technology=tech)



createColorBins <- function(name){
  vals <- pal_simpsons("springfield")(16)[c(8,9,2,3,9,4,16)]
  names(vals) <- c("array", "lr", "sr", "lr.all", "lr.>1", "lr.>5", "lr.>1_all.HQ")
  colscale <- list()
  colscale[["color"]] <- scale_colour_manual(name = name, values = vals)
  colscale[["fill"]] <- scale_fill_manual(name = name, values = vals)
  colscale
}
createColorBinsFull <- function(name){
  vals <- pal_simpsons("springfield")(16)[c(8,9,2,3,9,4,16)]
  names(vals) <- c("array", "long reads", "short reads", "lr.all", "lr.>1", "lr.>5", "lr.>1_all.HQ")
  colscale <- list()
  colscale[["color"]] <- scale_colour_manual(name = name, values = vals)
  colscale[["fill"]] <- scale_fill_manual(name = name, values = vals)
  colscale
}

colBin1 <- createColorBins("lbin")
colBin2 <- createColorBins("tech")
colBin3 <- createColorBinsFull("Technology")

Fig3b <- ggplot(c1, 
       aes(x=sc_bin, y=perc, color=Technology)) +
  geom_vline(xintercept = 0.5, linetype="dotted", size=1.3) + 
  geom_point(size=2) +
  geom_line(aes(group = Technology)) +
  colBin3[["color"]] +
  theme_pubclean() +
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        plot.margin = unit(c(0,1,1,1), "lines"),
        text = element_text(size=18)) +
  labs(x="Array-based score") +
  facet_wrap(~ qbin, scales = "free_y", nrow=2)

Fig3c <- ggviolin(siglrr %>% 
                    filter(len >=500) %>%
           mutate(GT=ifelse(VaPoR_GT=="", "no evidence", VaPoR_GT),
                  GT = factor(GT, levels=c( "1/1", "0/1", "0/0", "no evidence")),
                  tech = gsub("lr", "long reads", tech)) %>%
             mutate(tech = gsub("sr", "short reads", tech)), 
         x="GT", y="hel.dist",
         fill = "GT", add = "boxplot", add.params = list(fill = "white"), order=c( "1/1", "0/1", "0/0", "no evidence")) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  geom_hline(yintercept = 0.25, linetype="dotted") +
  scale_fill_uchicago() +
  theme_classic() +
  theme(legend.position = "none",
        plot.margin = unit(c(0,1,1,1), "lines"),
        text = element_text(size=18)) +
  facet_grid(~tech) +
  labs(y="Array-based score", x="Long read genotype by VaPoR") 


### SUPPL FIG S3 FOR CNV LOCUS PERCENTAGE SPAN
## for different lr bins 
a1 <- McnvRcol %>% filter(lr_bin=="above1")
aa1 <- a1 %>% 
  dplyr::filter(ntech>1, len>=500) %>%
  dplyr::mutate(len_bin=cut(lenkb, breaks = c(0, 5, 100, 500, 2000), labels = c("500bp-5kb", "5-100kb", "100-500kb", ">500kb")),
         perc_bin=cut(totperc, breaks = seq(0,110, by=10), labels = seq(10,110, by=10))) %>%
  dplyr::group_by(tech, perc_bin, len_bin) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::group_by(tech, len_bin) %>%
  dplyr::mutate(n_total = sum(n),
         perc=n/n_total*100)

aa1b <- a1 %>% 
  dplyr::filter(ntech>1, len>=500) %>%
  dplyr::mutate(perc_bin=cut(totperc, breaks = seq(0,110, by=10), labels = seq(10,110, by=10))) %>%
  dplyr::group_by(tech, perc_bin, ltech) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::group_by(tech, ltech) %>%
  dplyr::mutate(n_total = sum(n),
                perc=n/n_total*100,
                ltech=factor(ltech, levels=c("array|lr|sr","array|lr", "array|sr", "lr|sr")))


FigS3a <-  ggplot(aa1, 
                  aes(x=perc_bin, y=perc, color=tech)) +
  geom_point() +
  geom_line(aes(group=tech)) +
  colBin2[["color"]] +
  theme_pubclean() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text=element_text(size=17),
        legend.position = "none") +
  labs(x="Percentage span") +
  facet_wrap(~ len_bin  , ncol=4)

FigS3b <-  ggplot(aa1b %>%
                    mutate(ltech=gsub("\\|", "+", ltech),
                  ltech=factor(ltech, levels = c("array+lr+sr", 
                                                 "lr+sr", "array+lr", "array+sr"))), 
       aes(x=perc_bin, y=perc, color=tech)) +
  geom_point() +
  geom_line(aes(group=tech)) +
  colBin2[["color"]] +
  theme_pubclean() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        text=element_text(size=17),
        legend.position = "none") +
  scale_x_discrete(labels=c("", "20", "", "40", "", "60", "", "80", "", "100")) + 
  labs(x="Percentage span") +
  facet_wrap(~ ltech  , ncol=4)

aa2 <- a1 %>% 
  dplyr::filter(ntech>1, len>=500) %>%
  dplyr::mutate(perc_bin=cut(totperc, breaks = seq(0,110, by=10), labels = seq(10,110, by=10))) %>%
  dplyr::group_by(tech, perc_bin, qbin) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::group_by(tech, qbin) %>%
  dplyr::mutate(n_total = sum(n),
         perc=n/n_total*100)

FigS3c <-  ggplot(aa2, 
                 aes(x=perc_bin, y=perc, color=tech)) +
  geom_point() +
  geom_line(aes(group=tech)) +
  colBin2[["color"]] +
  theme_pubclean() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        text=element_text(size=17)) +
  scale_x_discrete(labels=c("", "20", "", "40", "", "60", "", "80", "", "100")) +
  labs(x="Percentage span") +
  facet_wrap(~ qbin , ncol = 1)
### SUPPL FIG S3 FOR CNV LOCUS PERCENTAGE SPAN
ggarrange(FigS3a, FigS3b, FigS3c, common.legend = F, ncol = 1, heights = c(4,4,5))



### SUPPL FIG S6 FOR SHORT READ SCORE vs LONG READ GENOTYPES
ab2 <- mcnvR %>% 
  filter(lr_bin=="above1") %>%
  dplyr::rename(lGT=VaPoR_GT) %>%
  dplyr::mutate(lGT=ifelse(lGT %in% c("1/1", "0/1"), "1/1_0/1", lGT)) %>% 
  dplyr::filter(lr_bin=="above1", !grepl(",", lGT), len>=500) %>% 
  dplyr::mutate(lGT=ifelse(lGT=="", "NA", lGT)) %>%
  dplyr::mutate(dh_bin=cut(dh.score, breaks = seq(-0.1,5, by=0.1), labels = seq(0,5, by=0.1))) %>%
  dplyr::group_by(tech, lGT, type) %>%
  dplyr::mutate(n_total=n()) %>%
  dplyr::group_by(dh_bin, tech, lGT, type, n_total) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(perc=n/n_total*100,
                dh_bin = as.numeric(as.character(dh_bin)),
                lGT=factor(lGT, levels=c("1/1_0/1", "0/0", "NA")))

FigS6a <- ggplot(ab2, 
                 aes(x=dh_bin, y=perc, color=tech)) +
  geom_vline(xintercept = 0.7, linetype="dotted") + 
  geom_vline(xintercept = 1, linetype="dotted", size=0.8) + 
  geom_vline(xintercept = 1.3, linetype="dotted") + 
  geom_point(size=1) +
  geom_line(aes(group=tech)) +
  coord_cartesian(xlim = c(0,3)) +
  colBin2[["color"]] +
  theme_pubclean() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size=15),
        legend.position = "left",
        plot.margin = unit(c(2,1,1,1), "lines")) +
  labs(x="Depth fold change score") +
  facet_wrap(~ type + lGT , scales = "free_y", nrow=2)


btp <- b1 %>%
  mutate(lbin=gsub("above", "lr.>", lbin)) %>%
  mutate(lbin=gsub("_HQ", "_all.HQ", lbin)) %>%
  mutate(lbin=gsub("^all$", "lr.all", lbin),
         lbin=factor(lbin, levels = c("array", "sr","lr.all", "lr.>1", "lr.>5", "lr.>1_all.HQ")))
FigS6b <- ggplot(btp %>%
                   filter(grepl("lr", lbin)) %>%
                   mutate(lbin=factor(lbin, levels = c("lr.all", "lr.>1", "lr.>5", "lr.>1_all.HQ"))), 
                 aes(x=dh_bin, y=perc, color=lbin)) +
  geom_vline(xintercept = 0.7, linetype="dotted") + 
  geom_vline(xintercept = 1, linetype="dotted", size=0.8) + 
  geom_vline(xintercept = 1.3, linetype="dotted") + 
  geom_point() +
  geom_line(aes(group = lbin)) +
  scale_x_continuous(breaks=c(0,0.5,1,1.5,2), labels = c(0,0.5,1,1.5,2)) +
  coord_cartesian(xlim = c(0,2.2)) +
  colBin1[["color"]] +
  theme_pubclean() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(2,2,1,4), "lines"),
        text = element_text(size=15),
        legend.position = "none") +
  labs(x="Depth fold change score", color = "lr score CNVL bins") +
  facet_wrap(~ type, scales = "free_y", ncol=1)


sb2 <- b2 %>%
  mutate(lbin=gsub("above", "lr.>", lbin)) %>%
  mutate(lbin=gsub("_HQ", "_all.HQ", lbin),
         lbin=factor(lbin, levels = c("array", "sr", "lr.>1", "lr.>1_all.HQ")))
FigS6c <-ggplot(sb2, 
                aes(x=dh_bin, y=perc, color=lbin)) +
  geom_vline(xintercept = 0.7, linetype="dotted") + 
  geom_vline(xintercept = 1, linetype="dotted", size=0.8) + 
  geom_vline(xintercept = 1.3, linetype="dotted") + 
  geom_point() +
  geom_line(aes(group = lbin)) +
  scale_x_continuous(breaks=c(0,0.5,1,1.5,2), labels = c(0,0.5,1,1.5,2)) +
  coord_cartesian(xlim = c(0,2.2)) +
  colBin1[["color"]] +
  theme_pubclean() +
  theme(axis.title.y = element_blank(),
        plot.margin = unit(c(1,2,0.5,1), "lines"),
        text = element_text(size=15),
        legend.position = "none") +
  labs(x="Depth fold change score", color = "long-read score binning") +
  facet_wrap(~ type + ntech, nrow=1)

### SUPPL FIG S6 FOR SHORT READ SCORE vs LONG READ GENOTYPES
ggarrange(ggarrange(FigS6a, FigS6b, widths = c(7, 2.5), common.legend = F), 
          FigS6c, ncol=1, heights = c(2, 1), common.legend = F)


######### Additional analysis: distance between the breakpoints
#### HIGH PERCENTAGES AND ARRAY BREAKPOINT PRECISION
ar1 <- mcnvR %>%
  dplyr::filter(lr_bin=="above1") %>%
  dplyr::left_join(McnvRcol %>% dplyr::select(cnvrid, ltech, lr_bin, tech, totperc)) %>%
  dplyr::filter(grepl("array\\|", ltech)) %>%
  dplyr::group_by(cnvrid, type) %>%
  dplyr::mutate(n_array=n_distinct(coordcnv[tech=="array"])) %>%
  dplyr::filter(n_array==1) %>%
  dplyr::group_by(cnvrid, type) %>%
  dplyr::mutate(array_coord=coordcnv[tech=="array"],
         chr_ar = i_chr[tech=="array"],
         start_ar = i_start[tech=="array"],
         end_ar = i_end[tech=="array"],
         perc_ar=perc[tech=="array"],
         len_ar=i_len[tech=="array"])%>%
  dplyr::filter(tech!="array") %>%
  dplyr::mutate(diff_start=abs(start_ar-i_start),
         diff_end=abs(end_ar-i_end), 
         sign_diff_start = sign(start_ar-i_start),
         sign_diff_end = sign(end_ar-i_end))

#### 'GOOD' AGREEMENT
## eg. array>50% and one of lr and sr >50%
ar1s <-ar1 %>%
  dplyr::filter(perc_ar > 50, perc > 50) %>%
  dplyr::group_by(cnvrid, array_coord, len_ar, perc_ar, type, tech, qbin, sign_diff_start, sign_diff_end) %>%
  dplyr::summarise(break_up = min(diff_start), 
            break_down = min(diff_end)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(break_left=break_up*sign_diff_start,
         break_right=break_down*sign_diff_end) %>%
  dplyr::select(-matches("sign"), -break_up, -break_down) %>%
  gather(-cnvrid, -array_coord, -len_ar, -perc_ar, -qbin, -type, -tech, key="bkpt", value = "distbp") %>%
  dplyr::mutate(distkb = distbp/1000,
         dist_perc=distbp/len_ar*100, 
         dist_bins = cut(distbp, breaks = c(-22000, seq(-1000, 1000, by=200),99000),  
                         labels = c(-2000, seq(-800, 1000, by=200), 2000)))

ar1ss <- ar1s %>%
  dplyr::select(-distkb, -dist_bins, -dist_perc) %>%
  spread(tech, distbp)
  

#### ARRAY BREAKPOINT RESOLUTION
ggplot(ar1s %>% 
         group_by(tech, bkpt, dist_bins) %>%
         dplyr::summarise(n=n()) %>%
         dplyr::mutate(n=ifelse(grepl("left", bkpt), -n,n)), 
       aes(x=dist_bins, y=n, fill=tech)) +
  geom_bar(stat="identity", position = "dodge", color="white") +
   coord_flip() +
  colBin2[["fill"]] +
  theme_pubclean() +
  xlab("Distance to closest array CNV breakpoint") +
  ylab("                    < Upstream breakpoint -- Downstream breakpoint >") 

#### COMPARE lr and sr for the same array CNVR
ggplot(ar1ss %>% 
         filter(!is.na(lr), !is.na(sr)), aes(x=lr, y=sr, color=type)) +
  geom_point(shape=21) +
  geom_abline(linetype="longdash") +
  theme_pubr() +
  facet_wrap(~qbin, scales = "free")


 

