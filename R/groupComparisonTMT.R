#' Finding differentially abundant proteins across conditions in TMT experiment
#'
#' Tests for significant changes in protein abundance across conditions based on a family of linear mixed-effects models in TMT experiment.
#' Experimental design of case-control study (patients are not repeatedly measured) is automatically determined based on proper statistical model.
#'
#' @param data the output of \code{\link{proteinSummarization}} function. It is a list with data frames `FeatureLevelData` and `ProteinLevelData`
#' @param contrast.matrix Comparison between conditions of interests. 1) default is "pairwise", which compare all possible pairs between two conditions. 2) Otherwise, users can specify the comparisons of interest. Based on the levels of conditions, specify 1 or -1 to the conditions of interests and 0 otherwise. The levels of conditions are sorted alphabetically.
#' @param moderated TRUE will moderate t statistic; FALSE (default) uses ordinary t statistic.
#' @param adj.method adjusted method for multiple comparison. "BH" is default.
#' @param remove_norm_channel TRUE(default) removes "Norm" channels from protein level data.
#' @param remove_empty_channel TRUE(default) removes "Empty" channels from protein level data.
#' @param save_fitted_models logical, if TRUE, fitted models will be added to
#' @inheritParams .documentFunction
#' 
#' @return a list that consists of the following elements: 
#' (1) ComparisonResult: statistical testing results; 
#' (2) FittedModel: the fitted linear models
#' 
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats anova lm median model.matrix na.omit p.adjust pt
#' @importFrom MSstats MSstatsContrastMatrix
#' 
#' @examples
#' data(input.pd)
#' # use protein.summarization() to get protein abundance data
#' quant.pd.msstats = proteinSummarization(input.pd,
#'                                        method="msstats",
#'                                        global_norm=TRUE,
#'                                        reference_norm=TRUE)
#'
#' test.pairwise = groupComparisonTMT(quant.pd.msstats, moderated = TRUE)
#' head(test.pairwise$ComparisonResult)
#'
#' # Only compare condition 0.125 and 1
#' levels(quant.pd.msstats$ProteinLevelData$Condition)
#'
#' # Compare condition 1 and 0.125
#' comparison=matrix(c(-1,0,0,1),nrow=1)
#'
#' # Set the nafmes of each row
#' row.names(comparison)="1-0.125"
#'
#' # Set the column names
#' colnames(comparison)= c("0.125", "0.5", "0.667", "1")
#' test.contrast = groupComparisonTMT(data = quant.pd.msstats,
#' contrast.matrix = comparison,
#' moderated = TRUE)
#' head(test.contrast$ComparisonResult)
#' 
groupComparisonTMT = function(
    data, contrast.matrix = "pairwise", 
    moderated = FALSE, adj.method = "BH",
    remove_norm_channel = TRUE, remove_empty_channel = TRUE,
    save_fitted_models = FALSE,
    use_log_file = TRUE, append = FALSE, 
    verbose = TRUE, log_file_path = NULL
){
    MSstatsConvert::MSstatsLogsSettings(
        use_log_file, append, verbose, log_file_path, 
        base = "MSstatsTMT_log_groupComparison_", pkg_name = "MSstatsTMT"
    )
    getOption("MSstatsTMTLog")("INFO", "MSstatsTMT - groupComparisonTMT function")
    getOption("MSstatsTMTLog")("INFO", paste("Moderated t-stat :", moderated))
    getOption("MSstatsTMTLog")("INFO", paste("Adjust p-value :", adj.method))
    getOption("MSstatsTMTLog")("INFO", paste("Remove empty channels :",
                                             remove_empty_channel))
    getOption("MSstatsTMTLog")("INFO", paste("Remove normalization channels :",
                               remove_norm_channel))
 
    summarized = MSstatsPrepareForGroupComparisonTMT(data$ProteinLevelData, 
                                                     remove_norm_channel,
                                                     remove_empty_channel)
    contrast_matrix = MSstats::MSstatsContrastMatrix(contrast.matrix,
                                                        unique(summarized$Group))
    fitted_models = MSstatsFitComparisonModelsTMT(summarized)
    FittedModel <- fitted_models$fitted_model
    names(FittedModel) <- fitted_models$protein
    
    fitted_models = MSstatsModerateTTest(summarized, fitted_models, moderated)
    testing_results = MSstatsGroupComparisonTMT(fitted_models, contrast_matrix)
    testing_results = MSstatsGroupComparisonOutputTMT(testing_results, 
                                                      adj.method)
    
    if(save_fitted_models){
        list(ComparisonResult = testing_results,
             ModelQC = NULL,
             FittedModel = FittedModel)   
    } else{
        list(ComparisonResult = testing_results,
             ModelQC = NULL,
             FittedModel = NULL)   
    }
}


#' Prepare output of proteinSummarization for group comparison
#' 
#' @param input output of proteinSummarization
#' @param remove_norm_channel if TRUE, "Norm" channel will be removed
#' @param remove_empty_channel, if TRUE, empty channel will be removed
#' 
#' @return data.table
#' 
#' @keywords internal
#' 
MSstatsPrepareForGroupComparisonTMT = function(input, remove_norm_channel,
                                               remove_empty_channel) {
    
    Condition <- Abundance <- NULL
    
    input = data.table::as.data.table(input)
    input = .checkGroupComparisonInput(input)
    
    if (remove_empty_channel & is.element("Empty", unique(input$Condition))) {
        input = input[Condition != "Empty",]
        input$Condition = factor(input$Condition)
    }
    if (remove_norm_channel & is.element("Norm", unique(input$Condition))) {
        input = input[Condition != "Norm",]
        input$Condition = factor(input$Condition)
    }
    input = input[!is.na(Abundance),]
    data.table::setnames(input, c("BioReplicate", "Condition"), 
                         c("Subject", "Group"))
    ## Ting: need to change later for time course design
    input$Group = factor(input$Group)
    input$Subject = factor(input$Subject)
    input$Run = factor(input$Run)
    input$Channel = factor(input$Channel)
    input$TechRepMixture = factor(input$TechRepMixture)
    input$Mixture = factor(input$Mixture)
    input$Protein = as.character(input$Protein)
    input    
}


#' Fit linear models for group comparison
#' 
#' @param input output of the MSstatsPrepareForGroupComparisonTMT function
#' 
#' @return list 
#' 
#' @keywords internal
#' 
MSstatsFitComparisonModelsTMT = function(input) {
    Abundance = Group = Protein = NULL
    Channel = Mixture = Run = Subject = TechRepMixture = NULL
    
    all_proteins = as.character(unique(input$Protein))
    num_proteins = length(all_proteins)
    linear_models = vector("list", num_proteins)
    
    annotation = unique(input[, list(Run, Channel, Subject,
                                     Group, Mixture, TechRepMixture)])
    has_single_subject = .checkSingleSubject(annotation)
    has_techreps = .checkTechReplicate(annotation)
    has_biomixtures = .checkMulBioMixture(annotation)
    has_single_run = .checkSingleRun(annotation)
    has_Repeated_Measures <- .checkRepeatedMeasures(annotation)
    
    if(has_biomixtures){
        msg = paste0("Design: ", length(unique(annotation$Mixture)), " mixtures.")
    } else{
        msg = paste0("Design: 1 mixture.")
    }
    getOption("MSstatsTMTLog")("INFO", msg)
    getOption("MSstatsTMTMsg")("INFO", msg)
    
    if(has_techreps){
        msg = paste0("Design: ", length(unique(annotation$TechRepMixture)), 
                     " technical replicated MS runs per mixture.")
    } else{
        msg = paste0("Design: 1 MS run per mixture.")
    }
    getOption("MSstatsTMTLog")("INFO", msg)
    getOption("MSstatsTMTMsg")("INFO", msg)
    
    if(has_single_subject){
        msg = paste0("Design: 1 subject per condition (No biological variation).")
        getOption("MSstatsTMTLog")("INFO", msg)
        getOption("MSstatsTMTMsg")("INFO", msg)
    } else{
        if(has_Repeated_Measures){
            msg = paste0("Design: repeated measures design (A biological subject is measured in multiple conditions).")
        } else{
            msg = paste0("Design: group comparison design (Different conditions contains different biological subjects).")
        }
        getOption("MSstatsTMTLog")("INFO", msg)
        getOption("MSstatsTMTMsg")("INFO", msg)
    }
    
    msg = paste0("Model fitting for ", num_proteins , " proteins.")
    getOption("MSstatsTMTLog")("INFO", msg)
    getOption("MSstatsTMTMsg")("INFO", msg)
    pb = txtProgressBar(max=num_proteins, style=3)
    for (i in seq_along(all_proteins)) {
        protein = all_proteins[i]
        single_protein = input[Protein == protein]
        model_fit_result = MSstatsComparisonModelSingleTMT(single_protein, 
                                                           protein)
        linear_models[[i]] = model_fit_result
        setTxtProgressBar(pb, i)
    }
    close(pb)
    rbindlist(linear_models)
}


#' Fit a linear model for group comparison for a single protein
#' 
#' @param single_protein protein-level data for a single protein (single element
#' of list created by the MSstatsPrepareForGroupComparisonTMT function)
#' @param protein_name name of a protein from the single_protein data.table
#' 
#' @return list
#'  
#' @keywords internal
#'  
MSstatsComparisonModelSingleTMT = function(single_protein, protein_name) {
    
    Run <- Channel <- Subject <- Group <- Mixture <- TechRepMixture <- NULL
    
    annotation = unique(single_protein[, list(Run, Channel, Subject,
                                              Group, Mixture, TechRepMixture)])
    has_single_subject = .checkSingleSubject(annotation)
    has_techreps = .checkTechReplicate(annotation)
    has_biomixtures = .checkMulBioMixture(annotation)
    has_single_run = .checkSingleRun(annotation)
    has_Repeated_Measures <- .checkRepeatedMeasures(annotation)
    
    fitted_model = .fitModelTMT(single_protein, has_single_subject, 
                                has_techreps, has_biomixtures, has_single_run, 
                                has_Repeated_Measures)
    result = .addVarianceInformation(fitted_model, protein_name)
    result
}


#' Moderate T statistic for group comparison
#' 
#' @param summarized protein-level data produced by the proteinSummarization function
#' @param fitted_models output of the MSstatsFitComparisonModelsTMT function
#' @param moderated if TRUE, moderation will be performed
#' 
#' @return list
#' 
#' @keywords internal
#' 
MSstatsModerateTTest = function(summarized, fitted_models, moderated) {
    
    variance_df <- variance <- Protein <- NULL
    
    eb_input_s2 = fitted_models[variance_df != 0 & !is.na(variance_df),
                                variance]
    eb_input_df = fitted_models[variance_df != 0 & !is.na(variance_df),
                                variance_df]
    if (moderated) {
        eb_fit = limma::squeezeVar(eb_input_s2, eb_input_df)
        df_prior = eb_fit$df.prior
        variance_prior = eb_fit$var.prior
    } else { ## ordinary t statistic
        variance_prior = 0
        df_prior = 0
    }
    fitted_models$df_prior = df_prior
    fitted_models$variance_prior = variance_prior
    fitted_models$total_df <- sum(eb_input_df, na.rm=TRUE)
    result = lapply(split(fitted_models, fitted_models$protein), as.list)
    result = lapply(result, function(single_fitted_model) {
        protein = single_fitted_model$protein
        single_fitted_model$data = summarized[Protein == protein]
        single_fitted_model
    })
}


#' Group comparison for TMT data
#' 
#' @param fitted_models output of the MSstatsModerateTTest function
#' @param contrast_matrix contrast matrix
#' 
#' @return data.table
#' 
#' @keywords internal
#'  
MSstatsGroupComparisonTMT = function(fitted_models, contrast_matrix) {
    msg = paste0("Testing for ", length(fitted_models) , " proteins:")
    getOption("MSstatsTMTLog")("INFO", msg)
    getOption("MSstatsTMTMsg")("INFO", msg)
    testing_results = vector("list", length(fitted_models))
    pb = txtProgressBar(max = length(fitted_models), style = 3)
    for (i in seq_along(fitted_models)) {
        testing_result = MSstatsTestSingleProteinTMT(fitted_models[[i]], 
                                                     contrast_matrix)
        testing_results[[i]] = testing_result
        setTxtProgressBar(pb, i)
    }
    close(pb)
    testing_results
}


#' Hypothesis tests for a single protein in TMT data
#' 
#' @param fitted_model single element of the MSstatsModerateTTest output
#' @param contrast_matrix contrast matrix
#' 
#' @return list
#' 
#' @keywords internal
#' 
MSstatsTestSingleProteinTMT = function(fitted_model, contrast_matrix) {
    single_protein = fitted_model[["data"]]
    groups = as.character(unique(single_protein$Group))
    groups = sort(groups)
    coefs = fitted_model[["coefs"]][[1]]
    s2 = fitted_model[["variance"]]
    s2_df = fitted_model[["variance_df"]]
    fit = fitted_model[["fitted_model"]][[1]]
    s2_prior = fitted_model[["variance_prior"]]
    df_prior = fitted_model[["df_prior"]]
    total_df = fitted_model[["total_df"]]
    protein = fitted_model[["protein"]]
    
    results = vector("list", nrow(contrast_matrix))
    if (!inherits(fit, "try-error")) {
        
        if (is.infinite(df_prior)) {
            s2_posterior <- s2_prior
        } else{
            s2_posterior = (s2_prior * df_prior + s2 * s2_df) / (df_prior + s2_df)
        }
        
        for (row_id in seq_len(nrow(contrast_matrix))) {
            contrast = contrast_matrix[row_id, , drop = FALSE]
            positive.groups = colnames(contrast)[contrast > 0]
            negative.groups = colnames(contrast)[contrast < 0]
            
            if (any(positive.groups %in% groups) & 
                any(negative.groups %in% groups)) {
                cm = .makeContrastSingleTMT(fit, contrast, single_protein, coefs)
                FC = (cm %*% coefs)[1, 1]
                if (inherits(fit, "lm")) {
                    se2.post = diag(t(cm) %*% summary(fit)$cov.unscaled %*% cm) * s2_posterior
                    df.post = s2_df + df_prior
                } else {
                    # Acknowlege: Tyler Bradshawthis contributed to this part of implementation
                    vcov = fit@vcov_beta
                    se2 = as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)
                    ## calculate posterior variance
                    vcov.post = fit@pp$unsc() * s2_posterior
                    se2.post = as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm)
                    ## calculate posterior df
                    g <- sapply(fit@Jac_list, function(gm) cm %*% gm %*% cm)
                    denom <- try(t(g) %*% fit@vcov_varpar %*% g, silent=TRUE)
                    if (inherits(denom, "try-error")) {
                        df.post = s2_df + df_prior
                    } else{
                        df.post = 2*(se2)^2/denom + df_prior
                    }
                }
                df.post <- pmin(df.post, total_df)
                
                t = FC / sqrt(se2.post)
                p = 2*pt(-abs(t), df = df.post)
                se = sqrt(se2.post)
                DF = df.post
                if (s2_df == 0) {
                    issue = "SingleMeasurePerCondition"
                } else{
                    issue = NA
                }
            } else {
                result = .handleMissingConditionTMT(single_protein, contrast[1, ])
                FC = result[["logFC"]]
                p = NA
                se = NA
                DF = NA
                issue = result[["issue"]]
            }
            
            results[[row_id]] = list(Protein = protein, Comparison = row.names(contrast),
                                     log2FC = FC, pvalue = p, SE = se, DF = DF, issue = issue)
            # This code cannot be extracted to a new function yet due to 
            # performance issues with environments
        }
    } else {
        # very few measurements so that the model is unfittable
        for (row_id in 1:nrow(contrast_matrix)) {
            contrast = contrast_matrix[row_id, , drop = FALSE]
            result = .handleMissingConditionTMT(single_protein, contrast[1, ])
            results[[row_id]] = list(Protein = protein, 
                                     Comparison = row.names(contrast),
                                     log2FC = result$logFC, pvalue = NA, SE = NA,
                                     DF = NA, issue = result$issue)
        }  
    } 
    data.table::rbindlist(results)
}


#' Combine testing results for individual proteins
#' 
#' @param testing_results output of the MSstatsGroupComparisonTMT function
#' @param adj_method method that will be used to adjust p-values for multiple comparisons
#' 
#' @return data.table
#' 
#' @keywords internal
#' 
MSstatsGroupComparisonOutputTMT = function(testing_results, adj_method) {
    
    adj.pvalue <- pvalue <- Protein <- Comparison <- log2FC <- SE <- DF <- issue <- NULL 
    
    testing_results = data.table::rbindlist(testing_results)
    testing_results$Protein = as.factor(testing_results$Protein)
    testing_results[, adj.pvalue := p.adjust(pvalue, adj_method), 
                    by = "Comparison"]
    testing_results[, list(Protein, Label = Comparison, log2FC, SE, DF,
                           pvalue, adj.pvalue, issue)]
}
