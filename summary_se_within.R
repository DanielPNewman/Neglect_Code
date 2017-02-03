# summary_se_within.R

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)

summary_se = function(data, measurevar, idvar, betweenvars) {
    require('dplyr')
    require('lazyeval')
    
    summarized = data %>%
        group_by_(.dots=c(idvar, betweenvars)) %>%
        dplyr::summarise_(.dots = list(
            MeasureVar = interp( ~mean( MEASUREVAR, na.rm=TRUE), MEASUREVAR = as.name(measurevar)  )
        )) %>%
        ungroup() %>%
        group_by_(.dots=betweenvars) %>%
        dplyr::summarise(N    = n(),
                         Mean = mean(MeasureVar, na.rm=TRUE),
                         SD   = sd(MeasureVar, na.rm=TRUE),
                         SE   = SD / sqrt(N) 
        ) %>%
        ungroup()
    
    colnames(summarized)[colnames(summarized)=='Mean'] = measurevar
    
    return( summarized )
    
}

summary_se_within = function(data, measurevar, idvar, withinvars= NULL, betweenvars= NULL) {
    require('dplyr')
    require('lazyeval')
    
    # Really just between subjects?
    if (is.null(withinvars)) {
        return( summary_se(data, measurevar, idvar, betweenvars) )
    }
    
    # Warn about impossibility of errorbars in mixed designs:
    if (!is.null(betweenvars)) {
        warning('Error bars cannot be accurately represented in mixed designs. ',
                'Treating each level of any between-subjects factors as separate experiment.')
    }
    
    # Collapse Multiple Observations in the Same Participant x Design Cell, Get Normed MeasureVar:
    normed = collapse_and_norm(data, measurevar, idvar, withinvars, betweenvars)
    
    # Get Correction factor:
    num_within_groups = prod( apply(data[,withinvars,drop=FALSE], MARGIN = 2, FUN = function(x) length(unique(x))) )
    correction_factor = sqrt( num_within_groups / (num_within_groups-1) )
    
    # Get Means, SDs, Etc:
    summarized = normed %>%
        group_by_(.dots= c(betweenvars, withinvars) ) %>%
        dplyr::summarise(N    = n(),
                         Mean = mean(MeasureVar, na.rm= TRUE),
                         SD   = sd(MeasureVarNormed, na.rm= TRUE),
                         SE   = SD / sqrt( N ),
                         CI   = SE * qt(.975, df = N-1) 
        ) %>%
        dplyr::mutate(SD = SD*correction_factor,
                      SE = SE*correction_factor,
                      CI = CI*correction_factor)
    
    colnames(summarized)[colnames(summarized)=='Mean'] = measurevar
    
    summarized
    
}

collapse_and_norm = function(data, measurevar, idvar, withinvars, betweenvars= NULL) {
    require('dplyr')
    require('lazyeval')
    
    # Collapse Multiple Observations in the Same Participant x Design Cell : 
    collapsed = data %>%
        group_by_(.dots= c(idvar, betweenvars, withinvars)) %>%
        dplyr::summarise_(.dots = list(MeasureVar = interp( ~mean( MEASUREVAR,na.rm=TRUE), MEASUREVAR = as.name(measurevar) )
        )) %>%
        ungroup() %>%
        group_by_(.dots= c(idvar, betweenvars) ) %>%
        dplyr::mutate(SubMean = mean(MeasureVar, na.rm=TRUE),
                      Diff    = MeasureVar - SubMean
        ) %>%
        ungroup() 
    
    # Get Normed Data:
    normed = collapsed %>%
        dplyr::mutate(MeasureVarNormed = (MeasureVar-SubMean) + mean(MeasureVar, na.rm=TRUE) ) 
}