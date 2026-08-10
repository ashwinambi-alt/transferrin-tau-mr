##############################################################################
# R2_PART3_FigS10.R — transferrin-specific controls forest plot (Fig S10)
# Style matches R1_FIGURES.R (Arial 9pt, COL palette, theme_forest, save_fig).
# Shows IVW (filled) and weighted-median (open) beta (95% CI) per control.
##############################################################################

.libPaths(c('C:/Users/ashwi/AppData/Local/R/win-library/4.6', .libPaths()))
suppressPackageStartupMessages({ library(ggplot2) })

PROJECT_DIR <- "C:/Users/ashwi/OneDrive/Documents/EB1A Docs/MR Paper"
REV  <- file.path(PROJECT_DIR,"results","reviewer_response")
SUPP <- file.path(PROJECT_DIR,"Final Figures and data","figures","supplementary")
FONT <- "Arial"
COL_PRIMARY<-"#1B4F8A"; COL_CONTROL<-"#2E7D32"; COL_NULL<-"#AAAAAA"; COL_RISK<-"#C62828"
mm2in <- function(mm) mm/25.4
save_fig <- function(path,plot,w,h,dpi=300){
  ggsave(paste0(path,".png"),plot,width=w,height=h,dpi=dpi,bg="white")
  ggsave(paste0(path,".pdf"),plot,width=w,height=h,device=cairo_pdf,bg="white")
  cat("  Saved:",basename(path),"\n") }

d <- read.csv(file.path(REV,"TableS7_transferrin_controls.csv"), stringsAsFactors=FALSE)
keep <- d[d$method %in% c("Inverse variance weighted","Weighted median"),]
keep$ci_lo <- keep$b - 1.96*keep$se
keep$ci_hi <- keep$b + 1.96*keep$se

# display order (top -> bottom) and colour by role
ord <- c("Transferrin saturation","Serum iron","Educational attainment","Height")
keep$outcome <- factor(keep$outcome, levels=rev(ord))
rolecol <- c("Transferrin saturation"=COL_CONTROL,"Serum iron"=COL_RISK,
             "Educational attainment"=COL_NULL,"Height"=COL_NULL)
keep$col <- rolecol[as.character(keep$outcome)]
# role factor drives the colour legend (colour = control type)
rolelab <- c("Transferrin saturation"="Positive control",
             "Serum iron"="Serum iron (DAG / mediator)",
             "Educational attainment"="Negative control",
             "Height"="Negative control")
keep$role <- factor(rolelab[as.character(keep$outcome)],
                    levels=c("Positive control","Serum iron (DAG / mediator)","Negative control"))
role_colors <- c("Positive control"=COL_CONTROL,
                 "Serum iron (DAG / mediator)"=COL_RISK,
                 "Negative control"=COL_NULL)
keep$is_ivw <- keep$method=="Inverse variance weighted"
# vertical offset so IVW and WM don't overlap
keep$ypos <- as.numeric(keep$outcome) + ifelse(keep$is_ivw, 0.15, -0.15)
lab_ivw <- keep[keep$is_ivw,]
lab_ivw$plab <- ifelse(lab_ivw$pval<1e-3, formatC(lab_ivw$pval,format="e",digits=2),
                       sprintf("p=%.2f",lab_ivw$pval))

theme_forest <- function(base=9) theme_minimal(base_size=base, base_family=FONT) +
  theme(panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
        panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.x=element_line(color="#EBEBEB",linewidth=0.2),
        axis.ticks=element_line(color="black",linewidth=0.3),
        plot.title=element_text(face="bold",size=base,margin=margin(b=6)),
        plot.caption=element_text(size=7,hjust=0,color="grey30",margin=margin(t=6)))

fig <- ggplot(keep, aes(x=b, y=ypos)) +
  geom_vline(xintercept=0, linetype="dashed", color="#888888", linewidth=0.4) +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi, color=role), height=0.12, linewidth=0.6) +
  geom_point(aes(shape=method, color=role), fill="white", size=2.6) +
  scale_shape_manual(values=c("Inverse variance weighted"=15,"Weighted median"=22), name=NULL) +
  scale_color_manual(values=role_colors, name="Control type") +
  # Only ONE legend: colour = control role (green/red/grey), shown as filled squares
  # (drop the CI line from the key). The filled-vs-open shape still distinguishes
  # IVW from weighted median in the plot, but its legend is suppressed because the
  # caption already states "Filled = IVW; open = weighted median" (avoids a
  # redundant, neutral-grey method key that doesn't match the coloured points).
  guides(color=guide_legend(override.aes=list(shape=15, linetype=0, fill=role_colors)),
         shape="none") +
  geom_text(data=lab_ivw, aes(x=ci_hi+0.03, label=plab), hjust=0, size=2.3,
            family=FONT, color=lab_ivw$col) +
  scale_y_continuous(breaks=1:4, labels=rev(ord)) +
  scale_x_continuous(limits=c(-0.85,0.45), breaks=seq(-0.8,0.4,0.2)) +
  labs(title="Transferrin-instrument controls (transferrin → control outcomes)",
       x=expression(beta*" (95% CI) per SD transferrin"), y=NULL,
       caption=paste0("Filled = IVW; open = weighted median. Positive control (transferrin saturation) is strongly negative as predicted\n",
                      "(shares the Benyamin cohort). Negative controls null. Transferrin→serum iron IVW is null (supports iron as a parallel\n",
                      "exposure, not a mediator); the positive weighted-median estimate reflects shared genetics at some loci (e.g. HFE).")) +
  theme_forest()

save_fig(file.path(SUPP,"FigS10_transferrin_controls"), fig, mm2in(180), mm2in(95))
cat("Done.\n")
