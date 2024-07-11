


statistical_analysis <- function(read_pos_file, caller) {
    
    
    df = fread(read_pos_file)
    
    colnames(df)[6] = 'sequences'
    
    df$ref_alt = paste(df$REF, df$ALT)
    
    
    df$outcome = ifelse(df$Found == 'Yes', 1, 0)
    
    
    tmp = df %>%
        
        select(REF, ALT, bin, Found) %>%
        
        gather(., key = 'RiskFactor', value = 'RF_Value', -Found) %>%
        
        freq_table(RiskFactor, RF_Value, Found) %>%
        
        mutate(
            
            result = paste0(n, ' (', prop, '%', ')')
            
        ) %>%
        
        select(-n, -prop) %>%
        
        spread(., key = 'Found', value = 'result')
    
    tmp$p_value = NA
    
    tmp$p_value[which(tmp$RiskFactor == 'REF')] = chisq.test(df$Found, df$REF)$p.value
    
    tmp$p_value[which(tmp$RiskFactor == 'ALT')] = chisq.test(df$Found, df$ALT)$p.value
    
    a = glm(ifelse(Found == 'Yes', 1, 0) ~ bin, data = df, family = 'binomial') %>%
        
        summary() %>%
        
        coefficients()
    
    tmp$p_value[which(tmp$RiskFactor == 'bin')] = a[2,4] 
    
    tmp$p_value = ifelse(
        
        tmp$p_value < 0.001,
        '<0.001',
        as.character(round(tmp$p_value, 3))
    )
    
    tmp$p_value[which(duplicated(tmp$RiskFactor))] = ''
    
    
    
}