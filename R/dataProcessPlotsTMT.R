#' Visualization for explanatory data analysis - TMT experiment
#'
#' To illustrate the quantitative data and quality control of MS runs,
#' dataProcessPlotsTMT takes the quantitative data  and summarized data from function `proteinSummarization` as input 
#' and generate two types of figures in pdf files as output :
#' (1) profile plot (specify "ProfilePlot" in option type), to identify the potential sources of variation for each protein;
#' (2) quality control plot (specify "QCPlot" in option type), to evaluate the systematic bias between MS runs and channels.
#'
#' @export
#' @import ggplot2
#' @importFrom graphics axis image legend mtext par plot.new title plot
#' @importFrom grDevices dev.off hcl pdf
#' @param data the output of \code{\link{proteinSummarization}} function. It is a list with data frames `FeatureLevelData` and `ProteinLevelData`
#' @param type choice of visualization. "ProfilePlot" represents profile plot of log intensities across MS runs.
#' "QCPlot" represents box plots of log intensities across channels and MS runs.
#' @param ylimUp upper limit for y-axis in the log scale.
#' FALSE(Default) for Profile Plot and QC Plot uses the upper limit as rounded off maximum of log2(intensities) after normalization + 3..
#' @param ylimDown lower limit for y-axis in the log scale. FALSE(Default) for Profile Plot and QC Plot uses 0..
#' @param x.axis.size size of x-axis labeling for "Run" and "channel in Profile Plot and QC Plot.
#' @param y.axis.size size of y-axis labels. Default is 10.
#' @param text.size size of labels represented each condition at the top of Profile plot and QC plot. Default is 4.
#' @param text.angle angle of labels represented each condition at the top of Profile plot and QC plot. Default is 0.
#' @param legend.size size of legend above Profile plot. Default is 7.
#' @param dot.size.profile size of dots in Profile plot. Default is 2.
#' @param ncol.guide number of columns for legends at the top of plot. Default is 5.
#' @param width width of the saved pdf file. Default is 10.
#' @param height height of the saved pdf file. Default is 10.
#' @param which.Protein Protein list to draw plots. List can be names of Proteins or order numbers of Proteins.
#' Default is "all", which generates all plots for each protein. For QC plot, "allonly" will generate one QC plot with all proteins.
#' @param originalPlot TRUE(default) draws original profile plots, without normalization.
#' @param summaryPlot TRUE(default) draws profile plots with protein summarization for each channel and MS run.
#' @param address the name of folder that will store the results. Default folder is the current working directory.
#' The other assigned folder has to be existed under the current working directory.
#' An output pdf file is automatically created with the default name of "ProfilePlot.pdf" or "QCplot.pdf".
#' The command address can help to specify where to store the file as well as how to modify the beginning of the file name.
#' If address=FALSE, plot will be not saved as pdf file but showed in window.
#' @return plot or pdf
#' @examples
#' data(input.pd)
#' quant.msstats = proteinSummarization(input.pd,
#'                                       method="msstats",
#'                                       global_norm=TRUE,
#'                                       reference_norm=TRUE)
#'
#' ## Profile plot
#' dataProcessPlotsTMT(data=quant.msstats,
#'                    type='ProfilePlot',
#'                    width = 21,
#'                    height = 7)
#'
#' ## NottoRun: QC plot
#' # dataProcessPlotsTMT(data=quant.msstats,
#'                     # type='QCPlot',
#'                     # width = 21,
#'                     # height = 7)
#'                     
dataProcessPlotsTMT = function(
    data, type, ylimUp = FALSE, ylimDown = FALSE,
    x.axis.size = 10, y.axis.size = 10, text.size = 4, text.angle = 90,
    legend.size = 7, dot.size.profile = 2, ncol.guide = 5, width = 10,
    height = 10, which.Protein = "all", originalPlot = TRUE, summaryPlot = TRUE,
    address = ""
) {
    data.peptide <- data$FeatureLevelData
    data.summarization <- data$ProteinLevelData
    common_groups = intersect(data.peptide$Condition, data.summarization$Condition)
    processed = .prepareDataForPlot(data.peptide, common_groups, "peptides")
    summarized = .prepareDataForPlot(data.summarization, common_groups, "proteins")
    
    checkmate::assertChoice(toupper(type), c("PROFILEPLOT", "QCPLOT"), 
                            .var.name = "type")
    if (address == FALSE) {
        if (which.Protein == "all") {
            stop(paste("** Cannnot generate all plots in a screen.",
                       "Please set one protein at a time."))
        } else if (length(which.Protein) > 1) {
            stop(paste("** Cannnot generate multiple plots in a screen.",
                       "Please set one protein at a time."))
        }
    }
    
    if (toupper(type) == "PROFILEPLOT") {
        .plotProfileTMT(processed, summarized, 
                        ylimUp, ylimDown, x.axis.size, y.axis.size, 
                        text.size, text.angle, legend.size, dot.size.profile, 
                        ncol.guide, width, height, which.Protein, 
                        originalPlot, summaryPlot,
                        address)
    }
    if (toupper(type) == "QCPLOT") {
        .plotQualityTMT(processed, 
                        ylimUp, ylimDown, x.axis.size, y.axis.size, 
                        text.size, text.angle, legend.size, dot.size.profile, 
                        ncol.guide, width, height, which.Protein,
                        address)
    }
}

#' @keywords internal
.prepareDataForPlot = function(input, common_groups, type) {
    
    Condition <- log2Intensity <- abundance <- Protein <- NULL
    
    input = data.table::as.data.table(input)
    input = input[Condition %in% common_groups]
    data.table::setnames(input, "ProteinName", "Protein", skip_absent = TRUE)
    input$Protein = factor(input$Protein)
    input$Condition = factor(input$Condition)
    if (type == "peptides") {
        input[,abundance := log2Intensity]
    }
    input
}


#' @importFrom MSstats theme_msstats savePlot
#' @keywords internal
.plotProfileTMT = function(processed, summarized, 
                           ylimUp, ylimDown, x.axis.size, y.axis.size, 
                           text.size, text.angle, legend.size, dot.size.profile, 
                           ncol.guide, width, height, which.Protein, 
                           originalPlot, summaryPlot,
                           address) {
    
    Protein <- Condition <- xorder <- Run <- NULL
    Channel <- PeptideSequence <- PSM <- cumGroupAxis <- NULL
    abundance <- Abundance <- analysis <- NULL
    
    if (which.Protein != "all") {
        chosen_proteins = getSelectedProteins(which.Protein, unique(processed$Protein))
        processed = processed[Protein %in% chosen_proteins]
        processed$Protein = factor(processed$Protein)
        summarized = summarized[Protein %in% chosen_proteins]
        summarized$Protein = factor(summarized$Protein)
    }
    yaxis.name = "Log2-intensities"
    if (is.numeric(ylimUp)) {
        y.limup = ylimUp
    } else {
        y.limup = ceiling(max(processed$abundance, na.rm = TRUE) + 3)
    }
    if (is.numeric(ylimDown)) {
        y.limdown = ylimDown
    } else {
        y.limdown = 0
    }
    
    all_proteins = unique(processed$Protein)
    processed = .getXAxisOrder(processed)
    tempGroupName = unique(processed[, list(xorder, Condition, Run, Channel)])
    tempGroupName = tempGroupName[order(xorder), ]
    groupline = .getGroupLabel(tempGroupName, y.limup)
    groupline.all = groupline
    unique(groupline.all$Condition)
    ## remove last condition for vertical line between groups
    groupline = groupline[!(Condition %in% levels(Condition)[nlevels(Condition)])]
    
    datamat = data.table::dcast(Protein + Channel ~ Run, 
                                data = summarized, 
                                value.var = "Abundance", keep = TRUE)
    summarized = data.table::melt(datamat, id.vars = c("Protein", "Channel"))
    data.table::setnames(summarized, c("variable", "value"),
                         c("Run", "Abundance"))
    summarized = merge(summarized, tempGroupName, by = c("Run", "Channel"))
    
    
    if (originalPlot) {
        savePlot(address, "ProfilePlot", width, height)
        message(paste0("Drew the Profile plot for ", length(all_proteins), " proteins."))
        pb = txtProgressBar(max = length(all_proteins), style=3)
        for (i in seq_along(all_proteins)) {
            single_protein = processed[Protein == all_proteins[i]]
            single_protein$PeptideSequence = factor(as.character(single_protein$PeptideSequence))
            single_protein$Charge = factor(as.character(single_protein$Charge))
            single_protein$PSM = factor(as.character(single_protein$PSM))
            single_protein$censored = factor(as.character(single_protein$censored))
            if (all(is.na(single_protein$abundance)) | 
                all(single_protein$abundance == 0)) {
                next()
            }
            
            pept_feat = unique(single_protein[, list(PeptideSequence, PSM)])
            pept_feat = pept_feat[order(PeptideSequence, PSM)]
            counts = pept_feat[, .(N = .N), by = "PeptideSequence"]$N
            s = rep(seq_along(counts), times = counts)
            ss = unlist(lapply(counts, function(x) seq(1, x)), FALSE, FALSE)
            
            ## for annotation of condition
            groupline.tmp = data.frame(groupline,
                                       "PSM" = unique(single_protein$PSM)[1],
                                       "PeptideSequence" = unique(single_protein$PeptideSequence)[1])
            groupline.all.tmp = data.frame(groupline.all,
                                           "PSM" = unique(single_protein$PSM)[1],
                                           "PeptideSequence" = unique(single_protein$PeptideSequence)[1])
            cbp = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
            check.length = length(unique(s)) %/% length(cbp)
            if ( check.length > 0 ){
                cbp = rep(cbp, times=check.length + 1)
            }
            
            ptemp = ggplot(aes_string(x = 'xorder', y = 'abundance',
                                      color = 'PSM', linetype = 'PSM'), data = single_protein) +
                facet_grid(~Run) +
                geom_point(data = single_protein, aes(shape=censored), size=dot.size.profile) +
                geom_point(size=dot.size.profile) +
                geom_line(size = 0.5) +
                scale_colour_manual(values=cbp[s]) +
                scale_linetype_manual(values = ss) +
                scale_shape_manual(values = c(16, 1)) +
                labs(title = unique(single_protein$Protein),
                     x = 'MS runs') +
                scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
                scale_x_continuous('MS runs') +
                geom_vline(data = groupline.tmp,
                           aes(xintercept = cumGroupAxis),
                           colour = "grey", linetype = "longdash") +
                geom_text(data = groupline.all.tmp,
                          aes(x = xorder, y = abundance, label = Condition),
                          size = text.size,
                          angle = text.angle,
                          color = "black") +
                theme_msstats("PROFILEPLOT", x.axis.size, y.axis.size, legend.size) +
                theme(axis.ticks.x = element_blank(),
                      axis.text.x = element_blank())+
                guides(color = guide_legend(title = paste("# peptide:", nlevels(single_protein$PeptideSequence)),
                                            title.theme = element_text(size = 13, angle = 0),
                                            keywidth = 0.4,
                                            keyheight = 0.1,
                                            default.unit = 'inch',
                                            ncol = ncol.guide),
                       linetype = guide_legend(title = paste("# peptide:", nlevels(single_protein$PeptideSequence)),
                                               title.theme = element_text(size = 13, angle = 0),
                                               keywidth = 0.4,
                                               keyheight = 0.1,
                                               default.unit = 'inch',
                                               ncol = ncol.guide))
            
            print(ptemp)
            setTxtProgressBar(pb, i)
        }
        close(pb)
        if (address != FALSE) {
            dev.off()
        }
    }
    
    if (summaryPlot) {
        savePlot(address, "ProfilePlot_wSummarization", width, height)
        message(paste0("Drew the Profile plot with summarization for ", length(all_proteins), " proteins."))
        pb = txtProgressBar(max = length(all_proteins), style=3)
        for (i in seq_along(all_proteins)) {
            single_protein = processed[Protein == all_proteins[i]]
            single_protein$PeptideSequence = factor(as.character(single_protein$PeptideSequence))
            single_protein$Charge = factor(as.character(single_protein$Charge))
            single_protein$PSM = factor(as.character(single_protein$PSM))
            if (all(is.na(single_protein$abundance)) | 
                all(single_protein$abundance == 0)) {
                next()
            }
            groupline.tmp = data.frame(groupline,
                                       "PSM" = unique(single_protein$PSM)[1],
                                       "PeptideSequence" = unique(single_protein$PeptideSequence)[1],
                                       "analysis" = 'Run summary')
            groupline.all.tmp = data.frame(groupline.all,
                                           "PSM" = unique(single_protein$PSM)[1],
                                           "PeptideSequence" = unique(single_protein$PeptideSequence)[1],
                                           "analysis" = 'Run summary')
            
            subrun = summarized[Protein == all_proteins[i], ]
            if (nrow(subrun) != 0) {
                quantrun = summarized[
                    Protein == all_proteins[i],
                    list(Protein, PeptideSequence = "Run summary", 
                         Charge = "Run summary", PSM = "Run summary",
                         Channel, Run, abundance = Abundance, xorder,
                         analysis = "Run summary")]
            } else {
                quantrun = data.table::data.table(
                    Protein = all_proteins[i], PeptideSequence = "Run summary",
                    Charge = "Run summary", PSM = "Run summary",
                    abundance = NA, Intensity = NA
                )
            }
            single_protein$analysis = "Processed feature-level data"
            
            final = rbind(single_protein[, colnames(quantrun), with = FALSE], 
                          quantrun)
            final$analysis = factor(final$analysis)
            final$PSM = factor(final$PSM)
            
            ptempall = ggplot(final, 
                              aes_string(x = "xorder", y = "abundance",
                                         color = "analysis", linetype = "PSM",
                                         size = "analysis")) +
                facet_grid(~Run) +
                geom_point(size = dot.size.profile) +
                geom_line(size = 0.5) +
                scale_colour_manual(values = c("lightgray", "darkred")) +
                scale_shape_manual(values = c(16)) +
                scale_size_manual(values = c(1.7, 2), guide = "none") +
                scale_linetype_manual(values = c(rep(1, times = length(unique(final$PSM))-1), 2), guide = "none") +
                labs(title = unique(single_protein$Protein),
                     x = 'MS runs') +
                scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
                geom_vline(data = groupline.tmp,
                           aes(xintercept = cumGroupAxis),
                           colour = "grey", linetype = "longdash") +
                geom_text(data = groupline.all.tmp,
                          aes(x = xorder, y = abundance, label = Condition),
                          size = text.size,
                          angle = text.angle,
                          color = "black") +
                theme_msstats("PROFILEPLOT", x.axis.size, y.axis.size, 
                              legend.size, legend.title = element_blank()) +
                theme(axis.ticks.x = element_blank(),
                      axis.text.x = element_blank())+
                guides(color = guide_legend(order = 1,
                                            title = NULL,
                                            label.theme = element_text(size = 10, angle = 0)))
            
            ## draw point again because some red summary dots could be hiden
            ptempall = ptempall + 
                geom_point(data = final, 
                           aes(x = xorder, y = abundance, 
                               size = analysis, color = analysis))
            print(ptempall)
            setTxtProgressBar(pb, i)
            
        } # end-loop for each protein
        close(pb)
        if (address!=FALSE) {
            dev.off()
        }
    }
}

#' @importFrom MSstats theme_msstats getSelectedProteins savePlot
#' @keywords internal
.plotQualityTMT = function(processed, 
                           ylimUp, ylimDown, x.axis.size, y.axis.size, 
                           text.size, text.angle, legend.size, dot.size.profile, 
                           ncol.guide, width, height, which.Protein,
                           address) {
    
    Condition <- cumGroupAxis <- xorder <- abundance <- Protein <- NULL
    
    yaxis.name = 'Log2-intensities'
    savePlot(address, "QCPlot", width, height)
    
    if (is.numeric(ylimUp)) {
        y.limup = ylimUp
    } else {
        y.limup = ceiling(max(processed$abundance, na.rm = TRUE) + 3)
    }
    if (is.numeric(ylimDown)) {
        y.limdown = ylimDown
    } else {
        y.limdown = 0
    }
    
    processed = .getXAxisOrder(processed)    
    tempGroupName = unique(processed[, c("Condition", "xorder", "Run", "Channel")])
    groupline = .getGroupLabel(tempGroupName, y.limup)
    groupline.all = groupline
    groupline = groupline[!(Condition %in% levels(Condition)[nlevels(Condition)])]
    
    if (which.Protein == "all" | which.Protein == "allonly") {
        message("Drew the Quality Contol plot(boxplot) over all proteins.")
        groupline.tmp = data.frame(groupline,
                                   "PSM" = unique(processed$PSM)[1],
                                   "PeptideSequence" = unique(processed$PeptideSequence)[1])
        groupline.all.tmp = data.frame(groupline.all,
                                       "PSM" = unique(processed$PSM)[1],
                                       "PeptideSequence" = unique(processed$PeptideSequence)[1])
        processed$xorder = factor(processed$xorder) # for boxplot x-axis
        
        ptemp = ggplot(aes_string(x = "xorder", y = "abundance"), data = processed) +
            facet_grid(~Run) +
            geom_boxplot(aes_string(fill = "Condition"), outlier.shape = 1, outlier.size = 1.5) +
            labs(title = "All proteins",
                 x = "MS runs") +
            scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
            geom_vline(data = groupline.tmp,
                       aes(xintercept = cumGroupAxis),
                       colour = "grey", linetype = "longdash") +
            geom_text(data = groupline.all.tmp,
                      aes(x = xorder, y = abundance, label = Condition),
                      size = text.size,
                      angle = text.angle,
                      color = "black") +
            theme_msstats("QCPLOT", x.axis.size, y.axis.size,
                          legend_size = NULL) +
            theme(axis.ticks.x = element_blank(),
                  axis.text.x = element_blank())
        print(ptemp)
    }
    
    if (which.Protein != "allonly") {
        if (which.Protein != "all") {
            chosen_proteins = getSelectedProteins(which.Protein, unique(processed$Protein))
            processed = processed[Protein %in% chosen_proteins]
            processed$Protein = factor(processed$Protein)
        }
        
        all_proteins = unique(processed$Protein)
        message(paste0("Drew the Quality Contol plot(boxplot) for each of ", length(all_proteins), " proteins."))
        pb = txtProgressBar(max = length(all_proteins), style=3)
        for (i in seq_along(all_proteins)) {
            single_protein = processed[Protein == all_proteins[i]]
            #single_protein = single_protein[!is.na(abundance), ]
            if (all(is.na(single_protein$abundance)) | 
                all(single_protein$abundance == 0)) {
                next()
            }
            
            ## for annotation of condition
            groupline.tmp = data.frame(groupline,
                                       "PSM" = unique(single_protein$PSM)[1],
                                       "PeptideSequence" = unique(single_protein$PeptideSequence)[1])
            groupline.all.tmp = data.frame(groupline.all,
                                           "PSM" = unique(single_protein$PSM)[1],
                                           "PeptideSequence" = unique(single_protein$PeptideSequence)[1])
            single_protein$xorder = factor(single_protein$xorder) # for boxplot, x-axis, xorder should be factor
            
            ptemp = ggplot(aes_string(x = 'xorder', y = 'abundance'), data = single_protein) +
                facet_grid(~Run) +
                geom_boxplot(aes_string(fill = 'Condition'), outlier.shape = 1, outlier.size = 1.5) +
                labs(title = unique(single_protein$Protein),
                     x = 'MS runs') +
                scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
                geom_vline(data = groupline.tmp,
                           aes(xintercept = cumGroupAxis),
                           colour = "grey", linetype = "longdash") +
                geom_text(data = groupline.all.tmp,
                          aes(x = xorder, y = abundance, label = Condition),
                          size = text.size,
                          angle = text.angle,
                          color = "black") +
                theme_msstats("QCPLOT", x.axis.size, y.axis.size,
                              legend_size = NULL) +
                theme(axis.ticks.x = element_blank(),
                      axis.text.x = element_blank())
            print(ptemp)
            setTxtProgressBar(pb, i)
        }
        close(pb)
    }
    if (address != FALSE) {
        dev.off()
    }
}

#' @keywords internal
.getXAxisOrder = function(processed) {
    
    Channel <- group.channel <- Run <- Condition <- NULL
    
    processed = processed[order(Run, Condition, Channel)]
    processed$group.channel = paste(processed$Condition, processed$Channel, sep = "_")
    xorder = unique(processed[, list(Run, group.channel)])
    xorder[, xorder := 1:.N, by = "Run"]
    processed = merge(processed, xorder, by = c("Run", "group.channel"))
    processed 
}

#' @keywords internal
.getGroupLabel = function(input, y.limup) {
    cumGroupAxis <- groupAxis <- NULL
    
    input[, Condition := factor(Condition, levels = unique(Condition), ordered = TRUE)]
    groupline = input[, list(groupAxis = .N), by = c("Condition", "Run")]
    groupline[, cumGroupAxis := cumsum(groupAxis) + 0.5, by = "Run"]
    groupline$xorder = groupline$cumGroupAxis - groupline$groupAxis / 2
    groupline$abundance = y.limup - 0.5
    groupline
}
