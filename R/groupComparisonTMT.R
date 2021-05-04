#" Finding differentially abundant proteins across conditions in TMT experiment
#"
#" Tests for significant changes in protein abundance across conditions based on a family of linear mixed-effects models in TMT experiment.
#" Experimental design of case-control study (patients are not repeatedly measured) is automatically determined based on proper statistical model.
#"
#" @export
#" @param data Name of the output of \code{\link{proteinSummarization}} function. It should have columns named Protein, Mixture, TechRepMixture, Run, Channel, Condition, BioReplicate, Abundance.
#" @param contrast.matrix Comparison between conditions of interests. 1) default is "pairwise", which compare all possible pairs between two conditions. 2) Otherwise, users can specify the comparisons of interest. Based on the levels of conditions, specify 1 or -1 to the conditions of interests and 0 otherwise. The levels of conditions are sorted alphabetically.
#" @param moderated TRUE will moderate t statistic; FALSE (default) uses ordinary t statistic.
#" @param adj.method adjusted method for multiple comparison. "BH" is default.
#" @param remove_norm_channel TRUE(default) removes "Norm" channels from protein level data.
#" @param remove_empty_channel TRUE(default) removes "Empty" channels from protein level data.
#" @return data.frame with result of inference
#" @examples
#" data(input.pd)
#" # use protein.summarization() to get protein abundance data
#" quant.pd.msstats = proteinSummarization(input.pd,
#"                                        method="msstats",
#"                                        global_norm=TRUE,
#"                                        reference_norm=TRUE)
#"
#" test.pairwise = groupComparisonTMT(quant.pd.msstats, moderated = TRUE)
#"
#" # Only compare condition 0.125 and 1
#" levels(quant.pd.msstats$Condition)
#"
#" # Compare condition 1 and 0.125
#" comparison=matrix(c(-1,0,0,1),nrow=1)
#"
#" # Set the nafmes of each row
#" row.names(comparison)="1-0.125"
#"
#" # Set the column names
#" colnames(comparison)= c("0.125", "0.5", "0.667", "1")
#" test.contrast = groupComparisonTMT(data = quant.pd.msstats,
#" contrast.matrix = comparison,
#" moderated = TRUE)
#" 
groupComparisonTMT = function(
    data, contrast.matrix = "pairwise", moderated = FALSE, adj.method = "BH",
    remove_norm_channel = TRUE, remove_empty_channel = TRUE
){
    # LOG: "MSstatsTMT - groupComparisonTMT function"
    # LOG: paste("Moderated t-stat :", moderated)
    # LOG: paste("Adjust p-value :", adj.method)
    # LOG: "Remove empty channels :", remove_empty_channel)
    # LOG: "Remove normalization channels :", remove_norm_channel)
    # TODO: switch to new input (list with two elements)
    summarized = MSstatsPrepareForGroupComparisonTMT(data, remove_norm_channel,
                                                     remove_empty_channel)
    fitted_models = MSstatsFitComparisonModelsTMT(summarized)
    fitted_models = MSstatsModerateTTest(summarized, fitted_models, moderated)
    testing_results = MSstatsGroupComparisonTMT(fitted_models, contrast_matrix)
    testing_results = MSstatsGroupComparisonOutputTMT(testing_results, 
                                                      adj.method)
    testing_results
}


MSstatsPrepareForGroupComparisonTMT = function(input, remove_norm_channel,
                                               remove_empty_channel) {
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

MSstatsFitComparisonModelsTMT = function(input) {
    Abundance = Group = Protein = NULL
    
    all_proteins = as.character(unique(input$Protein))
    num_proteins = length(all_proteins)
    linear_models = vector("list", num_proteins)
    
    # LOG: paste0("Model fitting for ", num_proteins , " proteins:")
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

MSstatsComparisonModelSingleTMT = function(single_protein, protein_name) {
    annotation = unique(single_protein[, list(Run, Channel, Subject,
                                              Group, Mixture, TechRepMixture)])
    has_single_subject = .checkSingleSubject(annotation)
    has_techreps = .checkTechReplicate(annotation)
    has_biomixtures = .checkMulBioMixture(annotation)
    has_single_run = .checkSingleRun(annotation)
    
    fitted_model = .fitModelTMT(single_protein, has_single_subject, 
                                has_techreps, has_biomixtures, has_single_run)
    result = .addVarianceInformation(fitted_model, protein_name)
    result
}

MSstatsModerateTTest = function(summarized, fitted_models, moderated) {
    if (moderated) {
        eb_input_s2 = fitted_models[variance_df != 0 & !is.na(variance_df),
                                    variance]
        eb_input_df = fitted_models[variance_df != 0 & !is.na(variance_df),
                                    variance_df]
        eb_fit = limma::squeezeVar(eb_input_s2, eb_input_df)
        if (is.infinite(eb_fit$df.prior)) {
            df_prior = 0
            variance_prior = 0
        } else{
            df_prior = eb_fit$df.prior
            variance_prior = eb_fit$var.prior
        }
    } else { ## ordinary t statistic
        variance_prior = 0
        df_prior = 0
    }
    fitted_models$df_prior = df_prior
    fitted_models$variance_prior = variance_prior
    result = lapply(split(fitted_models, fitted_models$protein), as.list)
    result = lapply(result, function(single_fitted_model) {
        protein = single_fitted_model$protein
        single_fitted_model$data = summarized[Protein == protein]
        single_fitted_model
    })
}


MSstatsGroupComparisonTMT = function(fitted_models, contrast_matrix) {
    # LOG: paste0("Testing for ", length(all_proteins) , " proteins:")
    testing_results = vector("list", length(fitted_models))
    pb = txtProgressBar(max = length(fitted_models), style = 3)
    for (i in seq_along(fitted_models)) {
        testing_result = MSstatsTestSingleProteinTMT(fitted_models[[i]], 
                                                     contrast_matrix)
        testing_results[[i]] = testing_result
        setTxtProgressBar(pb, i)
    }
    close(pb)
    data.table::rbindlist(testing_results)
}


MSstatsTestSingleProteinTMT = function(fitted_model, contrast_matrix) {
    single_protein = fitted_model[["data"]]
    groups = as.character(unique(single_protein$Group))
    groups = sort(groups) # sort the groups based on alphabetic order
    coefs = fitted_model[["coefs"]][[1]]
    s2 = fitted_model[["variance"]]
    s2_df = fitted_model[["variance_df"]]
    fit = fitted_model[["fitted_model"]][[1]]
    s2_prior = fitted_model[["variance_prior"]]
    df_prior = fitted_model[["df_prior"]]
    protein = fitted_model[["protein"]]
    
    results = vector("list", nrow(contrast_matrix))
    if (!inherits(fit, "try-error")) {
        s2_posterior = (s2_prior * df_prior + s2 * s2_df) / (df_prior + s2_df)
        if (!inherits(fit, "lm")) { # prepare for df calculation
            rho = list() # environment containing info about model
            rho = suppressMessages(.rhoInit(rho, fit, TRUE)) # save lmer outcome in rho envir variable
            vss = .vcovLThetaL(fit)
        } else {
            rho = NULL
            vss = NULL
        }
         
        for (row_id in seq_len(nrow(contrast_matrix))) {
            contrast = contrast_matrix[row_id, , drop = FALSE]
            positive.groups = colnames(contrast)[contrast > 0]
            negative.groups = colnames(contrast)[contrast < 0]
            
            if (any(positive.groups %in% groups) & 
                any(negative.groups %in% groups)) {
                cm = .makeContrastSingleTMT(fit, contrast, single_protein, coefs)
                FC = (cm %*% coefs)[, 1]
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
                    g = .mygrad(function(x)  vss(t(cm), x)$varcor, c(rho$thopt, rho$sigma))
                    denom = try(t(g) %*% fit@vcov_varpar %*% g, silent=TRUE)
                    if (inherits(denom, "try-error")) {
                        df.post = s2_df + df_prior
                    } else{
                        df.post = 2*(se2)^2/denom + df_prior
                    }
                }
                
                t = FC / sqrt(se2.post)[1, 1]
                p = 2*pt(-abs(t), df = df.post)
                se = sqrt(se2.post)[1, 1]
                DF = df.post[1, 1]
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
            
            results[[row_id]] = list(Protein = protein, comparison = row.names(contrast),
                 log2FC = FC, pvalue = p, SE = se, DF = DF, issue = issue)
            # results[[row_id]] = .handleSingleContrastTMT(contrast, fit, single_protein, 
            #                                              coefs, protein, groups, s2_posterior, rho, vss,
            #                                              df_prior, s2_df)
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


MSstatsGroupComparisonOutputTMT = function(testing_results, adj_method) {
    testing_results$Protein = as.factor(testing_results$Protein)
    testing_results[, adj.pvalue := p.adjust(pvalue, adj_method), 
                    by = "comparison"]
    testing_results[, list(Protein, Label = comparison, log2FC, SE, DF,
                           pvalue, adj.pvalue, issue)]
}
