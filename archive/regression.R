regress_kt <- function (dataframe, assay = FALSE, round = FALSE) 
{
    dataframe <- easysorter:::ensure_long(dataframe)
    dataframe <- dplyr::filter(dataframe, is.finite(phenotype), 
                               !is.na(strain))
    dataframe <- dataframe %>% 
        dplyr::filter(trait != "n.sorted")
    if (assay) {
        # check if round is true
        if(round) {
            # do assay regression for each round
            dataframe <- dataframe %>% 
                dplyr::group_by(condition, round) %>% 
                dplyr::filter(length(unique(assay)) > 1)
            regressed <- dataframe %>% 
                dplyr::group_by(condition, trait, round) %>% 
                do(fit = lm(phenotype ~ assay - 1, .))
            withresids <- regressed %>% 
                broom::augment(fit) %>% 
                dplyr::ungroup() %>% 
                dplyr::left_join(dataframe, ., by = c("condition", "trait", "phenotype", "assay", "round")) %>% 
                dplyr::distinct(condition, trait, phenotype, strain, row, col, plate, .keep_all = T) %>% 
                dplyr::rename(resid = .resid)
            regressedframe <- withresids %>% 
                dplyr::mutate(phenotype = resid) %>% 
                dplyr::select(-resid, -.fitted, -.se.fit, -.hat, -.sigma, -.cooksd, -.std.resid)
            
            # do round regression
            assayregressed <- regressedframe %>% 
                dplyr::group_by(condition, trait) %>% 
                do(fit = lm(phenotype ~ round - 1, .))
            assayresids <- assayregressed %>% 
                broom::augment(fit) %>% 
                dplyr::ungroup() %>% 
                dplyr::left_join(regressedframe, ., by = c("condition", "trait", "phenotype", "round")) %>% 
                dplyr::distinct(condition, trait, phenotype, strain, row, col, plate, .keep_all = T) %>% 
                dplyr::rename(resid = .resid)
            roundregressed <- assayresids %>% 
                dplyr::mutate(phenotype = resid) %>% 
                dplyr::select(-resid, -.fitted, -.se.fit, -.hat, -.sigma, -.cooksd, -.std.resid)
        } else {
            dataframe <- dataframe %>% 
                dplyr::group_by(condition) %>% 
                dplyr::filter(length(unique(assay)) > 1)
            regressed <- dataframe %>% 
                dplyr::group_by(condition, trait) %>% 
                do(fit = lm(phenotype ~ assay - 1, .))
            withresids <- regressed %>% 
                broom::augment(fit) %>% 
                dplyr::ungroup() %>% 
                dplyr::left_join(dataframe, ., by = c("condition", "trait", "phenotype", "assay")) %>% 
                dplyr::distinct(condition, trait, phenotype, strain, row, col, plate, .keep_all = T) %>% 
                dplyr::rename(resid = .resid)
            regressedframe <- withresids %>% 
                dplyr::mutate(phenotype = resid) %>% 
                dplyr::select(-resid, -.fitted, -.se.fit, -.hat, -.sigma, -.cooksd, -.std.resid)
        }
    }
    else {
        data <- dataframe %>% 
            dplyr::filter(!is.na(control))
        controls <- dataframe %>% 
            dplyr::filter(is.na(control) | control == "None")
        controls$control <- controls$condition
        moltendata <- data
        moltencontrols <- controls %>% 
            dplyr::group_by(strain, control, trait, assay) %>% 
            dplyr::summarize(controlphenotype = mean(phenotype,  na.rm = TRUE))
        fusedmoltendata <- dplyr::left_join(moltendata, moltencontrols, by = c("strain", "control", "trait", "assay")) %>% 
            dplyr::filter(!is.na(phenotype), !is.na(controlphenotype))
        regressed <- fusedmoltendata %>% 
            dplyr::group_by(condition, trait) %>% 
            dplyr::do(fit = lm(phenotype ~ controlphenotype -  1, .))
        withresids <- regressed %>% broom::augment(fit) %>% 
            dplyr::ungroup() %>% 
            dplyr::left_join(fusedmoltendata, ., by = c("condition", "trait", "phenotype", "controlphenotype")) %>% 
            dplyr::distinct(condition, trait, phenotype, controlphenotype, strain, row, col, plate, .keep_all = T) %>% 
            dplyr::rename(resid = .resid)
        regressedframe <- withresids %>% 
            dplyr::mutate(phenotype = resid) %>% 
            dplyr::select(-resid, -controlphenotype)
    }
    return(regressedframe)
}
