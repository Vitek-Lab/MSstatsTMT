#' Visualization for explanatory data analysis - TMT experiment
#'
#' To illustrate the quantitative data and quality control of MS runs,
#' dataProcessPlotsTMT takes the quantitative data from converter functions (\code{\link{PDtoMSstatsTMTFormat}}, \code{\link{MaxQtoMSstatsTMTFormat}}, \code{\link{SpectroMinetoMSstatsTMTFormat}}) as input
#' and generate two types of figures in pdf files as output :
#' (1) profile plot (specify "ProfilePlot" in option type), to identify the potential sources of variation for each protein;
#' (2) quality control plot (specify "QCPlot" in option type), to evaluate the systematic bias between MS runs.
#'
#' @export
#' @import ggplot2
#' @importFrom graphics axis image legend mtext par plot.new title plot
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom dplyr mutate
#' @importFrom reshape2 dcast
#' @param data.psm name of the data with PSM-level, which can be the output of converter functions(\code{\link{PDtoMSstatsTMTFormat}}, \code{\link{MaxQtoMSstatsTMTFormat}}, \code{\link{SpectroMinetoMSstatsTMTFormat}}).
#' @param data.summarization name of the data with protein-level, which can be the output of \code{\link{proteinSummarization}} function.
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
#' quant.msstats <- proteinSummarization(input.pd,
#'                                          method="msstats",
#'                                          normalization=TRUE)
#'
#' ## Profile plot
#' dataProcessPlotsTMT(data.psm=input.pd,
#'                      data.summarization=quant.msstats,
#'                      type='ProfilePlot',
#'                      width = 21,
#'                      height = 7)
#'
#' ## NottoRun: QC plot
#' # dataProcessPlotsTMT(data.psm=input.pd,
#'                     # data.summarization=quant.msstats,
#'                     # type='QCPlot',
#'                     # width = 21,
#'                     # height = 7)

dataProcessPlotsTMT <- function(data.psm = data.psm,
                                 data.summarization = data.summarization,
                                 type = type,
                                 ylimUp = FALSE,
                                 ylimDown = FALSE,
                                 x.axis.size = 10,
                                 y.axis.size = 10,
                                 text.size = 4,
                                 text.angle = 90,
                                 legend.size = 7,
                                 dot.size.profile = 2,
                                 ncol.guide = 5,
                                 width = 10,
                                 height = 10,
                                 which.Protein = "all",
                                 originalPlot = TRUE,
                                 summaryPlot = TRUE,
                                 address = "") {

    Condition = Run = xorder = Channel = NULL
    groupAxis = cumGroupAxis = abundance = analysis = NULL
    
    datafeature <- data.psm
    datarun <- data.summarization

    # conditions in feature data
    fea.conds <- as.character(unique(datafeature$Condition))
    # conditions in protein data
    run.conds <- as.character(unique(datarun$Condition))
    
    # only keep the overlapped conditions between feature data and protein data
    shared.conds <- intersect(fea.conds, run.conds)
    datafeature <- datafeature[datafeature$Condition %in% shared.conds,]
    datarun <- datarun[datarun$Condition %in% shared.conds,]
    
    # make sure condition is factor
    datafeature$Condition <- factor(datafeature$Condition)
    datarun$Condition <- factor(datarun$Condition)
    
    colnames(datafeature)[colnames(datafeature) == 'ProteinName'] <- 'Protein'
    datafeature$Protein <- factor(datafeature$Protein)
    datarun$Protein <- factor(datarun$Protein)

    ## feature level data : log2 transform
    datafeature$abundance <- log2(datafeature$Intensity)
    datafeature[!is.na(datafeature$Intensity) &
                    datafeature$Intensity < 1, 'abundance'] <- 0

    if (length(setdiff(toupper(type), c(toupper("ProfilePlot"), toupper("QCPlot")))) != 0) {
        stop(paste0("Input for type=", type,
                    ". However,'type' should be one of \"ProfilePlot\", \"QCPlot\"."))
    }

    if (address == FALSE){
        ## here I used == FALSE, instead of !address. Because address can be logical or characters.
        if (which.Protein == 'all') {
            stop('** Cannnot generate all plots in a screen. Please set one protein at a time.')
        } else if (length(which.Protein) > 1) {
            stop('** Cannnot generate multiple plots in a screen. Please set one protein at a time.')
        }
    }

    ## Profile plot ##
    ## ---------------
    if (toupper(type) == "PROFILEPLOT") {

        ## choose Proteins or not
        if (which.Protein != "all") {
            ## check which.Protein is name of Protein
            if (is.character(which.Protein)) {
                temp.name <- which.Protein

                ## message if name of Protein is wrong.
                if (length(setdiff(temp.name,unique(datafeature$Protein))) > 0) {
                    stop(paste0("Please check protein name. Data set does not have this protein. - ",
                                toString(temp.name)))
                }
            }

            ## check which.Protein is order number of Protein
            if (is.numeric(which.Protein)) {
                temp.name <- levels(datafeature$Protein)[which.Protein]

                ## message if name of Protein is wrong.
                if (length(levels(datafeature$Protein)) < max(which.Protein)) {
                    stop(paste0("Please check your ion of proteins. There are ",
                                length(levels(datafeature$Protein))," proteins in this dataset."))
                }
            }

            ## use only assigned proteins
            datafeature <- datafeature[which(datafeature$Protein %in% temp.name), ]
            datafeature$Protein <- factor(datafeature$Protein)

            datarun <- datarun[which(datarun$Protein %in% temp.name), ]
            datarun$Protein <- factor(datarun$Protein)
        }

        ## assign upper or lower limit
        y.limup <- ceiling(max(datafeature$abundance, na.rm = TRUE) + 3)

        if (is.numeric(ylimUp)) {
            y.limup <- ylimUp
        }

        y.limdown <- 0
        if (is.numeric(ylimDown)) {
            y.limdown <- ylimDown
        }

        datafeature <- datafeature[with(datafeature, order(Run, Condition, Channel)), ]
        datafeature$Run <- factor(datafeature$Run)
        datarun$Run <- factor(datarun$Run)

        ## !! important: order of x-axis
        ## can be reorder by group and then channel, WITHIN Run
        ## first make new column for x-axis
        datafeature$group.channel <- paste(datafeature$Condition, datafeature$Channel, sep = "_")

        ## not sure better way for coding
        ## potentially change it.
        datafeature$xorder <- NA

        for (k in 1:length(unique(datafeature$Run))) {

            runid <- unique(datafeature$Run)[k]
            datafeature[datafeature$Run == runid, ]$xorder <- factor(datafeature[datafeature$Run == runid, ]$group.channel,
                                                                     levels <- unique(datafeature[datafeature$Run == runid, ]$group.channel),
                                                                     labels <- seq(1, length(unique(datafeature[datafeature$Run == runid, ]$group.channel))))
        }

        ## check
        ## unique(datafeature[datafeature$Run == '5', c('Channel', 'Condition', 'Run', 'xorder','group.channel')])

        ## need to make data.frame with same variables for condition name
        datafeature$xorder <- as.numeric(datafeature$xorder)
        ## keep unique information for x-axis labeling. will be used in plotting
        tempGroupName <- unique(datafeature[, c("Condition", "xorder", "Run", "Channel")])

        ## count # per condition per Run
        #groupline <- unique(datafeature[, c('Condition', 'Run')])
        #groupline$groupAxis <- as.numeric(xtabs(~Condition+Run, tempGroupName))
        groupline <- tempGroupName %>% dplyr::group_by(Condition, Run) %>% dplyr::mutate(groupAxis = n())
        groupline <- groupline %>% dplyr::select(-xorder, -Channel)
        groupline <- groupline[!duplicated(groupline), ]

        ## make accumurated # as condition increase
        groupline <- groupline %>% dplyr::group_by(Run) %>% dplyr::mutate(cumGroupAxis = cumsum(groupAxis))

        groupline$cumGroupAxis <- groupline$cumGroupAxis + 0.5

        ## add coordinate for group id
        groupline$xorder <- groupline$cumGroupAxis - groupline$groupAxis / 2
        groupline$abundance <- y.limup - 0.5

        ## save all information, for labeling group in plot
        groupline.all <- groupline

        ## remove last condition for vertical line between groups
        groupline <- groupline[-which(groupline$Condition %in% levels(groupline$Condition)[nlevels(groupline$Condition)]), ]

        ## need to fill in incomplete rows for Runlevel data
        haverun <- FALSE

        if (sum(is.element(colnames(datarun), "Run")) != 0) {
            datamat <- dcast(Protein + Channel ~ Run, data = datarun, value.var = 'Abundance', keep = TRUE)

            datarun <- melt(datamat, id.vars=c('Protein', 'Channel'))
            colnames(datarun)[colnames(datarun) %in% c("variable", "value")] <- c('Run', 'Abundance')

            ## match x axis order
            datarun <- merge(datarun, tempGroupName, by = c('Run', 'Channel'))

            haverun <- TRUE
        }

        ## save the plots as pdf or not
        ## If there are the file with the same name, add next numbering at the end of file name

        ## y-axis labeling
        yaxis.name <- 'Log2-intensities'

        if (originalPlot) {
            if (address != FALSE) {
                allfiles <- list.files()

                num <- 0
                filenaming <- paste0(address, "ProfilePlot")
                finalfile <- paste0(address, "ProfilePlot.pdf")

                while (is.element(finalfile, allfiles)) {
                    num <- num + 1
                    finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
                }

                pdf(finalfile, width = width, height = height)
            }

            ## factoring for run, channel, condition should be done before loop

            for (i in 1:nlevels(datafeature$Protein)) {

                sub <- datafeature[datafeature$Protein == levels(datafeature$Protein)[i], ]
                sub$PeptideSequence <- factor(as.character(sub$PeptideSequence))
                sub$Charge <- factor(as.character(sub$Charge))
                sub$PSM <- factor(as.character(sub$PSM))

                ## if all measurements are NA,
                if (nrow(sub) == sum(is.na(sub$abundance))) {
                    message(paste0("Can't the Profile plot for ", unique(sub$Protein),
                                   "(", i, " of ", length(unique(datafeature$Protein)),
                                   ") because all measurements are NAs."))
                    next()
                }

                if (nrow(sub) == sum(!is.na(sub$abundance) & sub$abundance == 0)) {
                    message(paste0("Can't the Profile plot for ", unique(sub$Protein),
                                   "(", i, " of ", length(unique(datafeature$Protein)),
                                   ") because all measurements are zeros."))
                    next()
                }

                ## seq for peptide and charge
                ## for seting up color and linetype
                b <- unique(sub[, c("PeptideSequence", "PSM")])
                b <- b[with(b, order(PeptideSequence, PSM)), ] ## add because if there are missing value, orders are different.

                temp1 <- xtabs(~PeptideSequence, b)
                ## unique charge id within peptide sequence, for line type
                ss <- NULL
                ## unique peptide sequence id, for color
                s <- NULL

                for (j in 1:length(temp1)) {
                    temp3 <- rep(j, temp1[j])
                    s <- c(s, temp3)
                    temp2 <- seq(1, temp1[j])
                    ss <- c(ss, temp2)
                }

                ## for annotation of condition
                groupline.tmp <- data.frame(groupline,
                                            "PSM" = unique(sub$PSM)[1],
                                            "PeptideSequence" = unique(sub$PeptideSequence)[1])

                groupline.all.tmp <- data.frame(groupline.all,
                                                "PSM" = unique(sub$PSM)[1],
                                                "PeptideSequence" = unique(sub$PeptideSequence)[1])

                ## 1st plot for original plot
                ptemp <- ggplot(aes_string(x = 'xorder', y = 'abundance',
                                           color = 'PSM', linetype = 'PSM'), data = sub) +
                    facet_grid(~Run) +
                    geom_point(size=dot.size.profile) +
                    geom_line(size = 0.5) +
                    scale_colour_manual(values = s) +
                    scale_linetype_manual(values = ss) +
                    scale_shape_manual(values = c(16)) +
                    labs(title = unique(sub$Protein),
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
                    theme(
                        panel.background = element_rect(fill = 'white', colour = "black"),
                        legend.key = element_rect(fill = 'white', colour = 'white'),
                        panel.grid.minor = element_blank(),
                        strip.background = element_rect(fill = 'gray95'),
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = y.axis.size, colour = "black"),
                        axis.ticks = element_line(colour = "black"),
                        axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
                        axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
                        title = element_text(size = x.axis.size + 8, vjust = 1.5),
                        legend.position = "top",
                        legend.text = element_text(size = legend.size)) +
                    guides(color = guide_legend(title = paste("# peptide:", nlevels(sub$PeptideSequence)),
                                                title.theme = element_text(size = 13, angle = 0),
                                                keywidth = 0.4,
                                                keyheight = 0.1,
                                                default.unit = 'inch',
                                                ncol = ncol.guide),
                           linetype = guide_legend(title = paste("# peptide:", nlevels(sub$PeptideSequence)),
                                                   title.theme = element_text(size = 13, angle = 0),
                                                   keywidth = 0.4,
                                                   keyheight = 0.1,
                                                   default.unit = 'inch',
                                                   ncol = ncol.guide))

                print(ptemp)

                message(paste("Drew the Profile plot for ", unique(sub$Protein),
                              "(", i, " of ", length(unique(datafeature$Protein)), ")"))
            }
            # end-loop for each protein

            if (address != FALSE) {
                dev.off()
            }

        } # end original plot

        ############################################
        ## 2st plot for original plot : summary
        ############################################

        if (summaryPlot) {
            if (address != FALSE) {
                allfiles <- list.files()

                num <- 0
                filenaming <- paste0(address, "ProfilePlot_wSummarization")
                finalfile <- paste0(address, "ProfilePlot_wSummarization.pdf")

                while (is.element(finalfile, allfiles)) {
                    num <- num + 1
                    finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
                }

                pdf(finalfile, width = width, height = height)
            }

            for (i in 1:nlevels(datafeature$Protein)) {

                sub <- datafeature[datafeature$Protein == levels(datafeature$Protein)[i], ]
                sub$PeptideSequence <- factor(as.character(sub$PeptideSequence))
                sub$Charge <- factor(as.character(sub$Charge))
                sub$PSM <- factor(as.character(sub$PSM))

                ## if all measurements are NA,
                if (nrow(sub) == sum(is.na(sub$abundance))) {
                    message(paste0("Can't the Profile plot for ", unique(sub$Protein),
                                   "(", i, " of ", length(unique(datafeature$Protein)),
                                   ") because all measurements are NAs."))
                    next()
                }

                if (nrow(sub) == sum(!is.na(sub$abundance) & sub$abundance == 0)) {
                    message(paste0("Can't the Profile plot for ", unique(sub$Protein),
                                   "(", i, " of ", length(unique(datafeature$Protein)),
                                   ") because all measurements are zeros."))
                    next()
                }

                ## for annotation of condition
                groupline.tmp <- data.frame(groupline,
                                            "PSM" = unique(sub$PSM)[1],
                                            "PeptideSequence" = unique(sub$PeptideSequence)[1],
                                            "analysis" = 'Run summary')

                groupline.all.tmp <- data.frame(groupline.all,
                                                "PSM" = unique(sub$PSM)[1],
                                                "PeptideSequence" = unique(sub$PeptideSequence)[1],
                                                "analysis" = 'Run summary')

                if (haverun) {
                    subrun <- datarun[datarun$Protein == levels(datafeature$Protein)[i], ]

                    if (nrow(subrun) != 0) {

                        quantrun <- sub[1, ]
                        quantrun[, 2:ncol(quantrun)] <- NA
                        quantrun <- quantrun[rep(seq_len(nrow(subrun))), ]

                        quantrun$Protein <- subrun$Protein
                        quantrun$PeptideSequence <- "Run summary"
                        quantrun$Charge <- "Run summary"
                        quantrun$PSM <- "Run summary"
                        quantrun$Channel <- subrun$Channel
                        quantrun$Run <- subrun$Run
                        quantrun$abundance <- subrun$Abundance
                        quantrun$xorder <- subrun$xorder

                    } else {
                        # if there is only one Run measured across all runs, no Run information for linear with censored
                        quantrun <- datafeature[1, ]
                        quantrun[, 2:ncol(quantrun)] <- NA

                        quantrun$Protein <- levels(datafeature$Protein)[i]
                        quantrun$PeptideSequence <- "Run summary"
                        quantrun$Charge <- "Run summary"
                        quantrun$PSM <- "Run summary"
                        quantrun$abundance <- NA
                        quantrun$Intensity <- NA
                    }

                    quantrun$analysis <- "Run summary"
                    sub$analysis <- "Processed feature-level data"

                    final <- rbind(sub, quantrun)
                    final$analysis <- factor(final$analysis)
                    final$PSM <- factor(final$PSM)

                    ptempall <- ggplot(aes_string(x = 'xorder', y = 'abundance',
                                                  color = 'analysis', linetype = 'PSM', size = 'analysis'), data = final) +
                        facet_grid(~Run) +
                        geom_point(size = dot.size.profile) +
                        geom_line(size = 0.5) +
                        scale_colour_manual(values = c("lightgray", "darkred")) +
                        scale_shape_manual(values = c(16)) +
                        scale_size_manual(values = c(1.7, 2), guide = "none") +
                        scale_linetype_manual(values = c(rep(1, times = length(unique(final$PSM))-1), 2), guide = "none") +
                        labs(title = unique(sub$Protein),
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
                        theme(
                            panel.background = element_rect(fill = 'white', colour = "black"),
                            legend.key = element_rect(fill = 'white', colour = 'white'),
                            panel.grid.minor = element_blank(),
                            strip.background = element_rect(fill = 'gray95'),
                            axis.ticks.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.text.y = element_text(size = y.axis.size, colour = "black"),
                            axis.ticks = element_line(colour = "black"),
                            axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
                            axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
                            title = element_text(size = x.axis.size + 8, vjust = 1.5),
                            legend.position = "top",
                            legend.text = element_text(size = legend.size),
                            legend.title = element_blank()) +
                        guides(color = guide_legend(order = 1,
                                                    title = NULL,
                                                    label.theme = element_text(size = 10, angle = 0)))

                    ## draw point again because some red summary dots could be hiden
                    ptempall <- ptempall + geom_point(data = final, aes(x = xorder, y = abundance, size = analysis, color = analysis))

                    print(ptempall)

                    message(paste("Drew the Profile plot with summarization for ", unique(sub$Protein),
                                  "(", i, " of ", length(unique(datafeature$Protein)), ")"))

                }

            } # end-loop for each protein

            if (address!=FALSE) {
                dev.off()
            }
        }
    } # end Profile plot


    ## QC plot (Quality control plot) ##
    ## ---------------------------------
    if (toupper(type) == "QCPLOT") {

        ## y-axis labeling
        yaxis.name <- 'Log2-intensities'

        ## save the plots as pdf or not
        ## If there are the file with the same name, add next numbering at the end of file name
        if (address != FALSE) {
            allfiles <- list.files()

            num <- 0
            filenaming <- paste0(address,"QCPlot")
            finalfile <- paste0(address,"QCPlot.pdf")

            while (is.element(finalfile, allfiles)) {
                num <- num + 1
                finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
            }

            pdf(finalfile, width = width, height = height)
        }

        ## assign upper or lower limit
        y.limup <- ceiling(max(datafeature$abundance, na.rm = TRUE) + 3)

        if (is.numeric(ylimUp)) {
            y.limup <- ylimUp
        }

        y.limdown <- 0
        if (is.numeric(ylimDown)) {
            y.limdown <- ylimDown
        }

        datafeature <- datafeature[with(datafeature, order(Run, Condition, Channel)), ]

        ## !! important: order of x-axis
        ## can be reorder by group and then channel, WITHIN Run
        ## first make new column for x-axis
        datafeature$group.channel <- paste(datafeature$Condition, datafeature$Channel, sep = "_")

        ## not sure better way for coding
        ## potentially change it.
        datafeature$xorder <- NA

        for (k in 1:length(unique(datafeature$Run))) {

            runid <- unique(datafeature$Run)[k]
            datafeature[datafeature$Run == runid, ]$xorder <- factor(datafeature[datafeature$Run == runid, ]$group.channel,
                                                                     levels <- unique(datafeature[datafeature$Run == runid, ]$group.channel),
                                                                     labels <- seq(1, length(unique(datafeature[datafeature$Run == runid, ]$group.channel))))
        }

        ## check
        ## unique(datafeature[datafeature$Run == 'PAMI-176_Mouse_K-T', c('Channel', 'Condition', 'Run', 'xorder','group.channel')])

        ## need to make data.frame with same variables for condition name
        datafeature$xorder <- as.numeric(datafeature$xorder)
        ## keep unique information for x-axis labeling. will be used in plotting
        tempGroupName <- unique(datafeature[, c("Condition", "xorder", "Run", "Channel")])

        ## count # per condition per Run
        #groupline <- unique(datafeature[, c('Condition', 'Run')])
        #groupline$groupAxis <- as.numeric(xtabs(~Condition+Run, tempGroupName))
        groupline <- tempGroupName %>% dplyr::group_by(Condition, Run) %>% dplyr::mutate(groupAxis = n())
        groupline <- groupline %>% dplyr::select(-xorder, -Channel)
        groupline <- groupline[!duplicated(groupline), ]

        ## make accumurated # as condition increase
        groupline <- groupline %>% dplyr::group_by(Run) %>% dplyr::mutate(cumGroupAxis = cumsum(groupAxis))

        groupline$cumGroupAxis <- groupline$cumGroupAxis + 0.5

        ## add coordinate for group id
        groupline$xorder <- groupline$cumGroupAxis - groupline$groupAxis / 2
        groupline$abundance <- y.limup - 0.5

        ## save all information, for labeling group in plot
        groupline.all <- groupline

        ## remove last condition for vertical line between groups
        groupline <- groupline[-which(groupline$Condition %in% levels(groupline$Condition)[nlevels(groupline$Condition)]), ]

        ## all protein
        if (which.Protein == 'all' | which.Protein == 'allonly') {

            ## for annotation of condition
            groupline.tmp <- data.frame(groupline,
                                        "PSM" = unique(datafeature$PSM)[1],
                                        "PeptideSequence" = unique(datafeature$PeptideSequence)[1])

            groupline.all.tmp <- data.frame(groupline.all,
                                            "PSM" = unique(datafeature$PSM)[1],
                                            "PeptideSequence" = unique(datafeature$PeptideSequence)[1])

            ## 1st plot for original plot
            ## for boxplot, x-axis, xorder should be factor
            datafeature$xorder <- factor(datafeature$xorder)

            ptemp <- ggplot(aes_string(x = 'xorder', y = 'abundance'), data = datafeature) +
                facet_grid(~Run) +
                geom_boxplot(aes_string(fill = 'Condition'), outlier.shape = 1, outlier.size = 1.5) +
                labs(title = 'All proteins',
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
                theme(
                    panel.background = element_rect(fill = 'white', colour = "black"),
                    legend.key = element_rect(fill = 'white', colour = 'white'),
                    panel.grid.minor = element_blank(),
                    strip.background = element_rect(fill = 'gray95'),
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.text.y = element_text(size = y.axis.size, colour = "black"),
                    axis.ticks = element_line(colour = "black"),
                    axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
                    axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
                    title = element_text(size = x.axis.size + 8, vjust = 1.5),
                    legend.position = "none")

            print(ptemp)

            message("Drew the Quality Contol plot(boxplot) for all proteins.")
        }

        ## each protein
        ## choose Proteins or not
        if (which.Protein != 'allonly') {
            if (which.Protein != "all") {
                ## check which.Protein is name of Protein
                if (is.character(which.Protein)) {

                    temp.name <- which.Protein

                    ## message if name of Protein is wrong.
                    if (length(setdiff(temp.name, unique(datafeature$Protein))) > 0) {
                        dev.off()
                        stop(paste0("Please check protein name. Data set does not have this protein. - ",
                                    toString(temp.name)))
                    }
                }

                ## check which.Protein is order number of Protein
                if (is.numeric(which.Protein)) {
                    temp.name <- levels(datafeature$Protein)[which.Protein]

                    ## message if name of Protein is wrong.
                    if (length(levels(datafeature$Protein)) < max(which.Protein)) {
                        dev.off()
                        stop(paste0("Please check your ion of proteins. There are ",
                                    length(levels(datafeature$Protein)), " proteins in this dataset."))
                    }
                }

                ## use only assigned proteins
                datafeature <- datafeature[which(datafeature$Protein %in% temp.name), ]
                datafeature$Protein <- factor(datafeature$Protein)
            }

            for (i in 1:nlevels(datafeature$Protein)) {
                sub <- datafeature[datafeature$Protein == levels(datafeature$Protein)[i], ]
                sub <- sub[!is.na(sub$abundance), ]

                ## if all measurements are NA,
                if (nrow(sub) == sum(is.na(sub$abundance))) {
                    message(paste("Can't the Quality Control plot for ", unique(sub$Protein),
                                  "(", i, " of ", length(unique(datafeature$Protein)),
                                  ") because all measurements are NAs."))
                    next()
                }

                ## for annotation of condition
                groupline.tmp <- data.frame(groupline,
                                            "PSM" = unique(sub$PSM)[1],
                                            "PeptideSequence" = unique(sub$PeptideSequence)[1])

                groupline.all.tmp <- data.frame(groupline.all,
                                                "PSM" = unique(sub$PSM)[1],
                                                "PeptideSequence" = unique(sub$PeptideSequence)[1])

                ## 1st plot for original plot
                ## for boxplot, x-axis, xorder should be factor
                sub$xorder <- factor(sub$xorder)

                ptemp <- ggplot(aes_string(x = 'xorder', y = 'abundance'), data = sub) +
                    facet_grid(~Run) +
                    geom_boxplot(aes_string(fill = 'Condition'), outlier.shape = 1, outlier.size = 1.5) +
                    labs(title = unique(sub$Protein),
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
                    theme(
                        panel.background = element_rect(fill = 'white', colour = "black"),
                        legend.key = element_rect(fill = 'white', colour = 'white'),
                        panel.grid.minor = element_blank(),
                        strip.background = element_rect(fill = 'gray95'),
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = y.axis.size, colour = "black"),
                        axis.ticks = element_line(colour = "black"),
                        axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
                        axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
                        title = element_text(size = x.axis.size + 8, vjust = 1.5),
                        legend.position = "none")

                print(ptemp)

                message(paste("Drew the Quality Contol plot(boxplot) for ", unique(sub$Protein),
                              "(", i, " of ", length(unique(datafeature$Protein)), ")"))

            } # end-loop
        }

        if (address != FALSE) {
            dev.off()
        }
    } # end QC plot
}
