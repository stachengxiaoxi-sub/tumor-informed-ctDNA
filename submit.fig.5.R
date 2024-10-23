#########Figure 6 and Extend Data 8a (JCO_wes_394 匹配 tumor_sample_uuid2)

##########1.看一下C3D1时间点，ctDNA 状态和PFS和OS的关系
##########2.结合baseline和C3D1看一下C3D1时间点，ctDNA 状态和PFS和OS的关系
##########1.看一下C3D1时间点，ctDNA 状态和PFS和OS的关系
##########2.结合baseline和C3D1看一下C3D1时间点，ctDNA 状态和PFS和OS的关系

library(readxl)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
##############读入数据########

################看一下C3D1时间点，tumor informed 方法下ctDNA status和PFS和OS的关系###########

# 三个数据分别为cli.all, mut, C1D3
# cli.all中的样本名为tumor_sample_uuid2和sample_id
# mut和C1D1中的样本名分别对应为Tumor_Sample_UUID和sample_id
# mut和C1D1中的基因符号为Hugo_Symbol和gene


cli.jco<- read_excel("F:/燃石医学/Choice-1-ITH 研究/analysis data/cli_465.xlsx", sheet = 1)
C3D1<-read_excel("F:/燃石医学/Choice-1-ITH 研究/analysis data/C3D1_panel.xlsx", sheet = 1)


#########ctDNA changes with the prognosis#####################


cli.2$ctDNA_inter <- ifelse(cli.2$ctDNA_positive_C1D1 == "yes" & cli.2$ctDNA_positive_C3D1 == "yes", "keep_yes",
                            ifelse(cli.2$ctDNA_positive_C1D1 == "no" & cli.2$ctDNA_positive_C3D1 == "no", "keep_no",
                                   ifelse(cli.2$ctDNA_positive_C1D1 == "yes" & cli.2$ctDNA_positive_C3D1 == "no", "turn_no",
                                          ifelse(cli.2$ctDNA_positive_C1D1 == "no" & cli.2$ctDNA_positive_C3D1 == "yes", "turn_yes", NA))))


#write.csv(cli.2,"F:/燃石医学/Choice-1-ITH 研究/analysis data/cli.C3D1.csv")
cli.2$treatment

# 使用subset函数筛选treatment为ICI+chemo和chemo的数据
ICI <- subset(cli.2, treatment %in% c("ICI"))
chemo <- subset(cli.2, treatment %in% c("chemo"))
# 使用table函数计算ctDNA_inter的比例
ctDNA_inter_ICI <- table(ICI$ctDNA_inter) / nrow(ICI)
ctDNA_inter_chemo <- table(chemo$ctDNA_inter) / nrow(chemo)

###Figure 6 b ctDNA positive 且 ICI-chemo 看一下ctDNA 清除与否OS和PFS

cli.3 <-cli.2 %>%
  filter(treatment=="ICI" & ctDNA_positive_C1D1=="yes") %>%
  filter(ctDNA_inter=="keep_yes"|ctDNA_inter=="turn_no")

#######PFS
cox_fit <- coxph(Surv(PFS_new, PFS_status_new) ~ ctDNA_inter, data = cli.3)
summary_cox<- summary(cox_fit)
summary_cox
survfit(data =cli.3, Surv(PFS_new, PFS_status_new) ~ ctDNA_inter) %>%
  ggsurvplot(risk.table = T,
             ggtheme = theme_prism(),
             pval = TRUE,
             conf.int = F,
             surv.median.line = "hv",
             pval.method = T,
             linetype = "strata",
             lwd = 3,
             risk.table.y.text = FALSE,
             # test.for.trend = T,
             palette = c( "#2c7fb8","#e6550d", "#fdc086","#A0CEE8"),
             legend = "top",
             legend.title = "",
             # legend.labs = c("chemo", "ICI+chemo"),
             #  legend.labs = c("chemo_ctDNA-", "ICI+chemo_ctDNA-", "chemo_ctDNA+", "ICI+chemo_ctDNA+"),
             font.main = c(14, "bold", "black"),
             font.x = c(14, "bold", "black"),
             font.y = c(14, "bold", "black"))

#OS
cox_fit <- coxph(Surv(OS_new, OS_status_new) ~ ctDNA_inter, data = cli.3)
summary_cox<- summary(cox_fit)
summary_cox
survfit(data =cli.3, Surv(OS_new, OS_status_new) ~ ctDNA_inter) %>%
  ggsurvplot(risk.table = T,
             ggtheme = theme_prism(),
             pval = TRUE,
             conf.int = F,
             surv.median.line = "hv",
             pval.method = T,
             linetype = "strata",
             lwd = 3,
             risk.table.y.text = FALSE,
             # test.for.trend = T,
             palette = c( "#2c7fb8","#e6550d", "#fdc086","#A0CEE8"),
             legend = "top",
             legend.title = "",
             # legend.labs = c("chemo", "ICI+chemo"),
             #  legend.labs = c("chemo_ctDNA-", "ICI+chemo_ctDNA-", "chemo_ctDNA+", "ICI+chemo_ctDNA+"),
             font.main = c(14, "bold", "black"),
             font.x = c(14, "bold", "black"),
             font.y = c(14, "bold", "black"))







#########Figure 6 d ICI-chemo OS 生存survival#####################
library(survival)
library(survminer)
library(ggprism)

cli.2$ctDNA_inter_new <- ifelse(cli.2$ctDNA_positive_C1D1 == "yes" & cli.2$ctDNA_positive_C3D1 == "yes", "yes",
                                ifelse(cli.2$ctDNA_positive_C1D1 == "no" & cli.2$ctDNA_positive_C3D1 == "no", "no",
                                       ifelse(cli.2$ctDNA_positive_C1D1 == "yes" & cli.2$ctDNA_positive_C3D1 == "no", "no",
                                              ifelse(cli.2$ctDNA_positive_C1D1 == "no" & cli.2$ctDNA_positive_C3D1 == "yes", "yes", NA))))


names(cli.2)
table(cli.2$treatment)
ICI.1<-cli.2 %>%
  filter(treatment == "ICI")

survfit(data =ICI.1, Surv(OS_new, OS_status_new) ~ ctDNA_inter_new+BOR) %>%
  ggsurvplot(risk.table = T,
             ggtheme = theme_prism(),
             pval = TRUE,
             conf.int = F,
             surv.median.line = "hv",
             pval.method = T,
             linetype = "strata",
             lwd = 3,
             risk.table.y.text = FALSE,
             # test.for.trend = T,
             palette = c("grey", "#fdc086","#2c7fb8","red","orange"),
             legend = "top",
             legend.title = "",
             # legend.labs = c("ctDNA keeping no", "ctDNA keeping yes","ctDNA turn no","ctDNA turn yes"),
             font.main = c(14, "bold", "black"),
             font.x = c(14, "bold", "black"),
             font.y = c(14, "bold", "black"))

###Figure 6 f ICI-chemo PFS 生存survival
survfit(data =ICI.1, Surv(PFS_new, PFS_status_new) ~ ctDNA_inter_new+BOR) %>%
  ggsurvplot(risk.table = T,
             ggtheme = theme_prism(),
             pval = TRUE,
             conf.int = F,
             surv.median.line = "hv",
             pval.method = T,
             linetype = "strata",
             lwd = 3,
             risk.table.y.text = FALSE,
             # test.for.trend = T,
             palette = c("grey", "#fdc086","#2c7fb8","red","orange"),
             legend = "top",
             legend.title = "",
             # legend.labs = c("ctDNA keeping no", "ctDNA keeping yes","ctDNA turn no","ctDNA turn yes"),
             font.main = c(14, "bold", "black"),
             font.x = c(14, "bold", "black"),
             font.y = c(14, "bold", "black"))
###Figure 6 f ICI-chemo PFS 分 BOR 看生存survival

library(ggprism)
names(ICI)
ICI.2<-ICI.1 %>%
  filter(BOR=="SD")
fit <- survfit(Surv(PFS_new, PFS_status_new) ~ ctDNA_inter_new, data = ICI.2)

p <- ggsurvplot(fit, risk.table = "abs_pct", 
                ggtheme = theme_prism(), pval = TRUE, 
                conf.int = FALSE, surv.median.line = "hv", 
                pval.method = TRUE, linetype = "strata", 
                # palette = c("#2c7fb8","#e6550d","#D8433F","#A0CEE8"), 
                risk.table.y.text.col = T, risk.table.y.text = T, 
                legend = "right")

p

###Figure 6 f ICI-chemo PFS 生存survival
survfit(data =ICI.1, Surv(PFS_new, PFS_status_new) ~ ctDNA_inter_new+BOR) %>%
  ggsurvplot(risk.table = T,
             ggtheme = theme_prism(),
             pval = TRUE,
             conf.int = F,
             surv.median.line = "hv",
             pval.method = T,
             linetype = "strata",
             lwd = 3,
             risk.table.y.text = FALSE,
             # test.for.trend = T,
             palette = c("grey", "#fdc086","#2c7fb8","red","orange"),
             legend = "top",
             legend.title = "",
             # legend.labs = c("ctDNA keeping no", "ctDNA keeping yes","ctDNA turn no","ctDNA turn yes"),
             font.main = c(14, "bold", "black"),
             font.x = c(14, "bold", "black"),
             font.y = c(14, "bold", "black"))
###Figure 6 f ICI-chemo PFS 分 BOR 看生存survival

library(ggprism)
names(ICI)
ICI.2<-ICI.1 %>%
  filter(BOR=="SD")
fit <- survfit(Surv(PFS_new, PFS_status_new) ~ ctDNA_inter_new, data = ICI.2)

p <- ggsurvplot(fit, risk.table = "abs_pct", 
                ggtheme = theme_prism(), pval = TRUE, 
                conf.int = FALSE, surv.median.line = "hv", 
                pval.method = TRUE, linetype = "strata", 
                # palette = c("#2c7fb8","#e6550d","#D8433F","#A0CEE8"), 
                risk.table.y.text.col = T, risk.table.y.text = T, 
                legend = "right")

p

###Figure 6 f ICI-chemo PFS 生存survival
survfit(data =ICI.1, Surv(PFS_new, PFS_status_new) ~ ctDNA_inter_new+BOR) %>%
  ggsurvplot(risk.table = T,
             ggtheme = theme_prism(),
             pval = TRUE,
             conf.int = F,
             surv.median.line = "hv",
             pval.method = T,
             linetype = "strata",
             lwd = 3,
             risk.table.y.text = FALSE,
             # test.for.trend = T,
             palette = c("grey", "#fdc086","#2c7fb8","red","orange"),
             legend = "top",
             legend.title = "",
             # legend.labs = c("ctDNA keeping no", "ctDNA keeping yes","ctDNA turn no","ctDNA turn yes"),
             font.main = c(14, "bold", "black"),
             font.x = c(14, "bold", "black"),
             font.y = c(14, "bold", "black"))
###Figure 6 g ICI-chemo OS 分 BOR 看生存survival

library(ggprism)
names(ICI)
ICI.2<-ICI.1 %>%
  filter(BOR=="SD")
fit <- survfit(Surv(OS_new, OS_status_new) ~ ctDNA_inter_new, data = ICI.2)

p <- ggsurvplot(fit, risk.table = "abs_pct", 
                ggtheme = theme_prism(), pval = TRUE, 
                conf.int = FALSE, surv.median.line = "hv", 
                pval.method = TRUE, linetype = "strata", 
                # palette = c("#2c7fb8","#e6550d","#D8433F","#A0CEE8"), 
                risk.table.y.text.col = T, risk.table.y.text = T, 
                legend = "right")

p

