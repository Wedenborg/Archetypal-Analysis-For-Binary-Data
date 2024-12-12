# create_onehotencoding.R

create_onehotencoding <- function(cdm){
  
  cohorts_tbl <- tbl(cdmCon(cdm), 'cohorts')
  cohorts <- get_cohorts(cdm)
  
  
  cohort_tab<- data.frame(name=sapply(cohorts,function(x) {ifelse(is.null(x$name),NA,x$name)}))
  cohort_tab <- cohort_tab %>% mutate(type=gsub("(_.*)", "", name)) %>% mutate(cohortName=str_replace(name,"^coho_|^trtm_|^como_",""))
  
  
  event_list <- list("Dx"=cohort_tab %>% filter(type=="como") %>% pull(cohortName),
                     "Tr"=cohort_tab %>% filter(type=="trtm") %>% pull(cohortName))
  
  target_cohorts <- c("cluster_headache","migraine")
  dat_list <- setNames(vector("list", length(target_cohorts)), target_cohorts)
  
  for (target_coho in c("cluster_headache","migraine")){

    one_hot_tbls <- setNames(vector("list",3),c("DxB_","DxA_","TrA_"))
    
    for (typew in c("Dx","Tr")){ #typew="Dx"
      
      event <- cohorts_tbl %>% filter(cohort_name %in% !!event_list[[typew]]) 
      
      for (prefixw in c("B_","A_")) { # prefixw="B_"
        
        if (typew == "Tr" & prefixw == "B_"){next}
        
        if (prefixw == "B_"){
          target <- cohorts_tbl %>% filter(cohort_name==paste0(target_coho,"_before_index") )
          
          comb <- left_join(target %>% select("subject_id","cohort_start_date","cohort_end_date"),
                            event %>% select("subject_id","event_name"="cohort_name","event_start_date"="cohort_start_date","event_end_date"="cohort_end_date"),
                            by="subject_id") %>% 
            mutate(event_present=ifelse(event_start_date>=cohort_start_date & event_end_date<cohort_end_date,TRUE,NA)) 
        }
        if (prefixw == "A_"){
          target <- cohorts_tbl %>% filter(cohort_name==paste0(target_coho,"_after_index"))
          
          comb <- left_join(target %>% select("subject_id","cohort_start_date","cohort_end_date"),
                            event %>% select("subject_id","event_name"="cohort_name","event_start_date"="cohort_start_date","event_end_date"="cohort_end_date"),
                            by="subject_id") %>% 
            mutate(event_present=ifelse(event_start_date>=cohort_start_date & event_end_date<=cohort_end_date,TRUE,NA))
        }  
        
        event_pres_long <- comb %>% 
          group_by(subject_id,event_name) %>% 
          summarise(event_present = any(event_present),.groups="drop") %>% 
          mutate(event_present=ifelse(is.na(event_present),FALSE,event_present)) %>% 
          mutate(event_name=ifelse(is.na(event_name),"None",event_name)) %>% 
          distinct() 
        
        tmp <- pivot_wider(event_pres_long,values_from=event_present,names_from=event_name,values_fill = list(event_present = FALSE),names_prefix=paste0(typew,prefixw)) 
        
        one_hot_tbls[[paste0(typew,prefixw)]] <- tmp %>% select(colnames(tmp)[!(colnames(tmp) %in% paste0(typew,prefixw,"None"))])
        
      }
      
    }
    
    
    # one_hot_encoded_CH <- full_join(one_hot_tbls[[1]],one_hot_tbls[[2]],by="subject_id")
    # 
    # for (le in setdiff(seq_along(one_hot_tbls),c(1,2))){
    #   one_hot_encoded_CH <- full_join(one_hot_encoded_CH,one_hot_tbls[[le]],by="subject_id")
    # }
    
    one_hot_encoded_CH <- purrr::reduce(one_hot_tbls, full_join, by = "subject_id")
    
    # add index date from target cohort
    one_hot_encoded_CH_wAG0 <- one_hot_encoded_CH %>% left_join(target %>% select(subject_id, cohort_start_date) ,by="subject_id")
    
    
    # add day_of_birth and gender (gender_source_concept_id) from person table (cdm$person)
    one_hot_encoded_CH_wAG1 <- one_hot_encoded_CH_wAG0 %>% 
      left_join(cdm$person %>% select(person_id, gender_concept_id, year_of_birth),by=c("subject_id"="person_id"))
    
    # get meaning female/male from concept table
    one_hot_encoded_CH_wAG2 <- one_hot_encoded_CH_wAG1 %>% 
      left_join(cdm$concept %>% 
                  filter(concept_name %in% c("MALE","FEMALE")) %>% 
                  mutate(concept_name=ifelse(concept_name=="MALE","Male","Female")) %>% 
                  select(concept_id,"GENDER"="concept_name"), 
                by=c("gender_concept_id"="concept_id"))
    
    # calculate AGE, 
    # and then remove: gender_concept_id, year_of_birth, cohort_start_date (, index_year)
    one_hot_encoded_CH_wAG3 <- one_hot_encoded_CH_wAG2 %>% 
      mutate(index_year=!!CDMConnector::datepart("cohort_start_date", interval = "year"),
             AGE= index_year-year_of_birth) %>% 
      select(-c(gender_concept_id,index_year,year_of_birth)) # select(-c(gender_concept_id,cohort_start_date,index_year,year_of_birth))
    
    # # put one-hot-encoding in the db write scheme 
    # one_hot_encoded_CH_wAG3 %>% compute(name=paste0("one_hot_encoded_",target_coho),
    #                                     temporary=FALSE,
    #                                     overwrite=TRUE)
    
    dat_list[[target_coho]] <- one_hot_encoded_CH_wAG3 %>% mutate(TARGET_COHORT = target_coho)
    
  }
  # combine the two cohorts to have clustering done on both (in case a subject is in two cohorts, the latest diagnose is used)
  dat <- purrr::reduce(dat_list, union_all)
  
  # remove subjects occurring in multiple target cohorts except for the occurrence with the latest index date (that is most likely their correct diagnose (migraine/CH))
  dat_unique <- dat %>%
    mutate(cohort_start_date = !!CDMConnector::asDate(cohort_start_date)) %>% # make sure the date is indeed interpreted as a date
    group_by(subject_id) %>%
    filter(cohort_start_date == max(cohort_start_date)) %>% 
    ungroup()

  # put one-hot-encoding in the db write scheme
  dat_unique %>% compute(name=paste0("one_hot_encoded_data"),
                                      temporary=FALSE,
                                      overwrite=TRUE)
  
}