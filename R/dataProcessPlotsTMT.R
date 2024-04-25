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
#' @importFrom htmltools save_html tagList div
#' @importFrom graphics axis image legend mtext par plot.new title plot
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom plotly ggplotly style add_trace plot_ly subplot
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
#' @param isPlotly Parameter to use Plotly or ggplot2. If set to TRUE, MSstats 
#' will save Plotly plots as HTML files. If set to FALSE MSstats will save ggplot2 plots
#' as PDF files
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
    x.axis.size = 10, y.axis.size = 10, text.size = 2, text.angle = 90,
    legend.size = 7, dot.size.profile = 2, ncol.guide = 5, width = 10,
    height = 10, which.Protein = "all", originalPlot = TRUE, summaryPlot = TRUE,
    address = "", isPlotly = FALSE
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
    warning("Avoid plotting all proteins as it can take a large amount of time 
            to download the files")
    if(isPlotly & address != FALSE) {
      print("Plots will be saved as .HTML file as plotly is selected, set isPlotly = FALSE, if 
            you want to generate PDF using ggplot2")
    }
    
    if (toupper(type) == "PROFILEPLOT") {
        plots <- .plotProfileTMT(processed, summarized, 
                        ylimUp, ylimDown, x.axis.size, y.axis.size, 
                        text.size, text.angle, legend.size, dot.size.profile, 
                        ncol.guide, width, height, which.Protein, 
                        originalPlot, summaryPlot,
                        address, isPlotly)
        plotly_plots = list()
        if(isPlotly) {
          og_plotly_plot = NULL
          summ_plotly_plot = NULL
          if("original_plot" %in% names(plots)) {
            for(i in seq_along(plots[["original_plot"]])) {
              plot_i <- plots[["original_plot"]][[paste("plot",i)]]
              og_plotly_plot <- .convertGgplot2Plotly(plot_i, tips=c("PSM","xorder","abundance","censored"))
              og_plotly_plot = .fixLegendPlotlyPlotsDataprocess(og_plotly_plot, "OriginalPlot")
              plotly_plots = c(plotly_plots, list(og_plotly_plot))
            }
          }
          if("summary_plot" %in% names(plots)) {
            for(i in seq_along(plots[["summary_plot"]])) {
              plot_i <- plots[["summary_plot"]][[paste("plot",i)]]
              summ_plotly_plot <- .convertGgplot2Plotly(plot_i)
              summ_plotly_plot = .fixLegendPlotlyPlotsDataprocess(summ_plotly_plot, "SummaryPlot")
              plotly_plots = c(plotly_plots, list(summ_plotly_plot))
            }
          }
          
          if(address != FALSE) {
            .savePlotlyPlotHTML(plotly_plots,address,"ProfilePlot" ,width, height)
          }
          plotly_plots
        }
    }
    else if (toupper(type) == "QCPLOT") {
        plots <- .plotQualityTMT(processed, 
                        ylimUp, ylimDown, x.axis.size, y.axis.size, 
                        text.size, text.angle, legend.size, dot.size.profile, 
                        ncol.guide, width, height, which.Protein,
                        address, isPlotly)
        plotly_plots = list()
        if(isPlotly) {
          for(i in seq_along(plots)) {
            plot <- plots[[i]]
            plotly_plot <- .convertGgplot2Plotly(plot)
            plotly_plots[[i]] = list(plotly_plot)
          }
          if(address != FALSE) {
            .savePlotlyPlotHTML(plotly_plots,address,"QCPlot" ,width, height)
          }
          plotly_plots <- unlist(plotly_plots, recursive = FALSE)
          plotly_plots
        }
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
                           address, isPlotly) {
    
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
    
    output_plots <- list()
    output_plots[["original_plot"]] = list()
    output_plots[["summary_plot"]] = list()
    if (originalPlot) {
        if(!isPlotly) {
          savePlot(address, "ProfilePlot", width, height)
        }
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
                # geom_point(size=dot.size.profile) +
                geom_line(size = 0.5) +
                scale_colour_manual(values=cbp[s]) +
                scale_linetype_manual(values = ss) +
                scale_shape_manual(values = c(16, 1),labels = c("Detected data", "Censored missing data")) +
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
              theme(axis.ticks.x = element_blank(),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
                    axis.text.x = element_blank(),strip.text.x = element_text(
                      size = 5, color = "black",angle = 15)) +
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
            output_plots[["original_plot"]][[paste("plot",i)]] <- ptemp
            setTxtProgressBar(pb, i)
        }
        close(pb)
        if (address != FALSE & !isPlotly) {
          dev.off()
        } 
    }
    
    if (summaryPlot) {
        if(!isPlotly) {
          savePlot(address, "ProfilePlot_wSummarization", width, height)
        }
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
                      axis.text.x = element_blank(),strip.text.x = element_text(
                        size = 5, color = "black",angle = 15),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14))+
                guides(color = guide_legend(order = 1,
                                            title = NULL,
                                            label.theme = element_text(size = 10, angle = 0)))
            
            ## draw point again because some red summary dots could be hiden
            ptempall = ptempall + 
                geom_point(data = final, 
                           aes(x = xorder, y = abundance, 
                               size = analysis, color = analysis))
            print(ptempall)
            output_plots[["summary_plot"]][[paste("plot",i)]] <- ptempall
            setTxtProgressBar(pb, i)
            
        } # end-loop for each protein
        close(pb)
        if (address != FALSE & !isPlotly) {
          dev.off()
        } 
    }
    if(isPlotly) {
      output_plots
    }
}

#' @importFrom MSstats theme_msstats getSelectedProteins savePlot
#' @keywords internal
.plotQualityTMT = function(processed, 
                           ylimUp, ylimDown, x.axis.size, y.axis.size, 
                           text.size, text.angle, legend.size, dot.size.profile, 
                           ncol.guide, width, height, which.Protein,
                           address, isPlotly) {
    
    Condition <- cumGroupAxis <- xorder <- abundance <- Protein <- NULL
    
    yaxis.name = 'Log2-intensities'
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
    
    if (!isPlotly) {
      savePlot(address, "QCPlot", width, height)
    }
    plots <- vector("list", length(unique(processed$Protein)) + 1) # +1 for all/allonly plot
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
                axis.text.x = element_blank(),strip.text.x = element_text(
                  size = 5, color = "black",angle = 15))
        
        print(ptemp)
        plots[[1]] = ptemp
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
                      axis.text.x = element_blank()
                      ,strip.text.x = element_text(
                        size = 5, color = "black",angle = 15)
                      ) 
            print(ptemp)
            plots[[i+1]] = ptemp # to accomodate all proteins
            setTxtProgressBar(pb, i)
        }
        close(pb)
    }
    if (address != FALSE) {
        dev.off()
    }
    if (isPlotly) {
      plots <- Filter(function(x) !is.null(x), plots) # remove if protein was not "all"
      plots
    }
}

facet_strip_bigger <- function(gp){
  if (!is.null(gp$x$layout$annotations)) {
    for (i in seq_along(gp$x$layout$annotations)) {
      if(gp$x$layout$annotations[[i]]$text != "Log2-intensities" && gp$x$layout$annotations[[i]]$text != "MS runs") {
        gp$x$layout$annotations[[i]]$font$size <- 7
        # gp$x$layout$annotations[[i]]$xanchor <- "center"
        gp$x$layout$annotations[[i]]$xshift <- 50
      }
    }
  }
  return(gp)
}

#' converter for plots from ggplot to plotly
#' @noRd
.convertGgplot2Plotly = function(plot, tips = "all") {
  converted_plot <- ggplotly(plot,tooltip = tips)
  converted_plot <- plotly::layout(
    converted_plot,
    width = 1800,   # Set the width of the chart in pixels
    height = 600,  # Set the height of the chart in pixels
    title = list(
      font = list(
        size = 18
      )
    ),
    legend = list(
      x = 0,     # Set the x position of the legend
      y = -0.25,    # Set the y position of the legend (negative value to move below the plot)
      orientation = "h",  # Horizontal orientation
      font = list(
        size = 12  # Set the font size for legend item labels
      ),
      title = list(
        font = list(
          size = 12  # Set the font size for the legend title
        )
      )
    )
  ) 
  converted_plot <- facet_strip_bigger(converted_plot)
  converted_plot
}

.savePlotlyPlotHTML = function(plots, address, file_name, width, height) {
  print("Saving plots as HTML")
  pb <- txtProgressBar(min = 0, max = 4, style = 3)
  
  setTxtProgressBar(pb, 1)
  file_name = getFileName(address, file_name, width, height)
  file_name = paste0(file_name,".html")
  
  setTxtProgressBar(pb, 2)
  doc <- .getPlotlyPlotHTML(plots, width, height)
  
  setTxtProgressBar(pb, 3)
  save_html(html = doc, file = file_name) # works but lib same folder
  
  setTxtProgressBar(pb, 4)
  zip(paste0(gsub("\\.html$", "", file_name),".zip"), c(file_name, "lib"))
  unlink(file_name)
  unlink("lib",recursive = T)
  
  close(pb)
}

.getPlotlyPlotHTML = function(plots, width, height) {
  doc <- tagList(lapply(plots,function(x) div(x, style = "float:left;width:100%;")))
  # Set a specific width for each plot
  plot_width <- 2000
  plot_height <- 600
  
  # Create a div for each plot with style settings
  divs <- lapply(plots, function(x) {
    div(x, style = paste0("width:", plot_width, "px; height:", plot_height, "px; margin: 10px;"))
  })
  
  # Combine the divs into a tagList
  doc <- tagList(divs)
  doc
}

.fixLegendPlotlyPlotsDataprocess = function(plot, type) {
  df <- data.frame(id = seq_along(plot$x$data), legend_entries = unlist(lapply(plot$x$data, `[[`, "name")))
  df$legend_group <- gsub("^\\((.*?),.*", "\\1", df$legend_entries)
  df$is_first <- !duplicated(df$legend_group)
  df$is_bool <- ifelse(grepl("TRUE|FALSE", df$legend_group), TRUE, FALSE)
  # df[nrow(df), "is_first"] <- FALSE 
  df$is_valid_column <- ifelse(grepl("Processed feature-level data|Run summary", df$legend_entries), TRUE, FALSE)
  plot$x$data[[nrow(df)]]$showlegend <- FALSE # remove text legend
  
  for (i in df$id) {
    is_first <- df$is_first[[i]]
    is_bool <- df$is_bool[[i]]
    plot$x$data[[i]]$name <- df$legend_group[[i]]
    plot$x$data[[i]]$legendgroup <- plot$x$data[[i]]$name
    if (!is_first) plot$x$data[[i]]$showlegend <- FALSE
      if(type == "SummaryPlot") {
        is_valid_column <- df$is_valid_column[[i]]
        if (!is_valid_column) plot$x$data[[i]]$showlegend <- FALSE
      }
    if(is_bool) plot$x$data[[i]]$showlegend <- FALSE
  }
  plot
}

getFileName = function(name_base, file_name, width, height) {
  all_files = list.files(".")
  if(file_name == 'ProfilePlot'){
    num_same_name = sum(grepl(paste0("^", name_base, file_name, "_[0-9]?"), all_files))
  } else {
    num_same_name = sum(grepl(paste0("^", name_base, file_name, "[0-9]?"), all_files))
  }
  if (num_same_name > 0) {
    file_name = paste(file_name, num_same_name + 1, sep = "_")
  }
  file_path = paste0(name_base, file_name)
  return(file_path)
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
