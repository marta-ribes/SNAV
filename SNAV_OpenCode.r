version #4.0.3

#Ignore warnings
options(warn=-1)

library(tidyverse)
library(readxl)
library(haven)
library(ggpubr)
library(psych)
library(MASS)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(stringr)
library(scales) #below the function that allows to put in log the scales:
ks <- function (x) { number_format(accuracy = 1,
                                   scale = 1/1000,
                                   suffix = "k",
                                   big.mark = ",")(x) }
library(epitools)
library(vtable)
library(janitor) #to exec Gemma's function
library(broom) #to exec Gemma's function
library(sjPlot) #show tables for regressions with ORs
library(corrplot)
library(GGally)
dodge = position_dodge(width=0.1)
library(Hmisc)

#Adjust function by Gemma Ruiz-Olalla i Miquel Vázquez
# '@export

chmi.stat.p_adjust <- function(
  phen,
  method = c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY',  'fdr', 'none'),
  arrange_by = c('raw_pval', 'adjust_pval'),
  trend = FALSE,
  stat = c('Hmisc'))
{
## nc == no converge
## np == no perform test
## ns == no signifincant

### argument 'stat'
  stat <- match.arg(stat, c('Hmisc'))

### update dataset with 'raw_pval' and 'adjust_pval'
  if('raw_pval' %in% names(phen)) {
    phen <- phen %>%
      mutate(
        raw_pval_text = factor(
          ifelse(is.nan(raw_pval), 'nc',
          ifelse(is.na(raw_pval), 'np',
          ifelse(raw_pval <= 0.0001, 'P < 0.0001',
          ifelse(raw_pval <= 0.001, 'P < 0.001',
          ifelse(raw_pval <= 0.05, paste0('P = ', round(raw_pval, 3)), 'ns')))))),
        raw_pval_signif = factor(
          ifelse(is.nan(raw_pval), 'nc',
          ifelse(is.na(raw_pval), 'np',
          ifelse(raw_pval <= 0.0001, '***',
          ifelse(raw_pval <= 0.001, '**',
          ifelse(raw_pval <= 0.01, '*',
          ifelse(raw_pval <= 0.05, '.', 'ns'))))))),
        adj_pval = as.numeric(
          ifelse(is.na(raw_pval), 'NA', p.adjust(raw_pval, method = method))),
        adj_pval_text = factor(
          ifelse(is.nan(adj_pval), 'nc',
          ifelse(is.na(adj_pval), 'np',
          ifelse(adj_pval <= 0.0001, 'P < 0.0001',
          ifelse(adj_pval <= 0.001, 'P < 0.001',
          ifelse(adj_pval <= 0.05, paste0('P = ', round(adj_pval, 3)), 'ns')))))),
        adj_pval_signif = factor(
          ifelse(is.nan(adj_pval), 'nc',
          ifelse(is.na(adj_pval), 'np',
          ifelse(adj_pval <= 0.0001, '***',
          ifelse(adj_pval <= 0.001, '**',
          ifelse(adj_pval <= 0.01, '*',
          ifelse(adj_pval <= 0.05, '.', 'ns'))))))))

        ### argument 'trend'
          if (trend == TRUE) {
            phen <- phen %>%
             mutate(
              raw_pval_text = factor(
                ifelse(is.nan(raw_pval), 'nc',
                ifelse(is.na(raw_pval), 'np',
                ifelse(raw_pval <= 0.0001, 'P < 0.0001',
                ifelse(raw_pval <= 0.001, 'P < 0.001',
                ifelse(raw_pval <= 0.05, paste0('P = ', round(raw_pval, 3)), 'ns')))))),
              raw_pval_signif = factor(
                ifelse(is.nan(raw_pval), 'nc',
                ifelse(is.na(raw_pval), 'np',
                ifelse(raw_pval <= 0.0001, '***',
                ifelse(raw_pval <= 0.001, '**',
                ifelse(raw_pval <= 0.01, '*',
                ifelse(raw_pval <= 0.05, '.', 'ns'))))))),
              adj_pval = as.numeric(
                ifelse(is.na(raw_pval), 'NA', p.adjust(raw_pval, method = method))),
              adj_pval_text = factor(
                ifelse(is.nan(adj_pval), 'nc',
                ifelse(is.na(adj_pval), 'np',
                ifelse(adj_pval <= 0.01, 'P < 0.01',
                ifelse(adj_pval <= 0.05, 'P < 0.05',
                ifelse(adj_pval <= 0.1, paste0('P = ', round(adj_pval, 3)), 'ns')))))),
              adj_pval_signif = factor(
                ifelse(is.nan(adj_pval), 'nc',
                ifelse(is.na(adj_pval), 'np',
                ifelse(adj_pval <= 0.01, '***',
                ifelse(adj_pval <= 0.05, '**',
                ifelse(adj_pval <= 0.1, '*', 'ns')))))))
        }
    } else {
     phen <- phen %>%
      mutate(
        adj_pval_text = factor(
          ifelse(is.nan(adj_pval), 'nc',
          ifelse(is.na(adj_pval), 'np',
          ifelse(adj_pval <= 0.0001, 'P < 0.0001',
          ifelse(adj_pval <= 0.001, 'P < 0.001',
          ifelse(adj_pval <= 0.05, paste0('P = ', round(adj_pval, 3)), 'ns')))))),
        adj_pval_signif = factor(
          ifelse(is.nan(adj_pval), 'nc',
          ifelse(is.na(adj_pval), 'np',
          ifelse(adj_pval <= 0.0001, '***',
          ifelse(adj_pval <= 0.001, '**',
          ifelse(adj_pval <= 0.01, '*',
          ifelse(adj_pval <= 0.05, '.', 'ns'))))))))
    }


### argument 'arrange_by'
  if('raw_pval' %in% names(phen) & arrange_by == 'raw_pval') {
    phen <- phen %>% arrange(raw_pval)

  } else if (arrange_by == 'adjust_pval') {
    phen <- phen %>% arrange(adj_pval)
  }


### return
  return(phen)
}



# Load function by Gemma Ruiz-Olalla i Miquel Vázquez to have coeficients from regression models as %

chmi.stat.lm_beta <- function(
phen,
l_group,
l_varsy, l_varsx,
    custom_formls = FALSE, #if we want to customize the model formulas. By default we don't.
    formls = NULL, #list of customized models (one per each 5 isotypes). By default, NULL.
    arg_multivar = c('univariate', 'multivariate'),
    arg_info_model = c('r_squared', 'adj_r_squared', 'sigma', 't_value', 'p_value',
      'df', 'log_lik', 'aic', 'bic', 'deviance', 'df_residual', 'NULL'),
    arg_info_beta = c('std_error' , 't_value_beta', 'NULL'),
    arg_filter_vars,
    gr_adj, #variables to group by for p.adjusting purposes (optional argument)
format_pval = TRUE,
    arg_pval_arrange = c('raw_pval', 'adjust_pval'),
    padj_method = c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY',  'fdr', 'none'),
    stat = c('janitor', 'broom', 'tidymodels',  'tidyverse'))
{
### stat
  stat <- match.arg(stat, c('janitor', 'broom', 'tidymodels', 'tidyverse'))

  # list of formulas
  l_formulas <- list()

  if(custom_formls == FALSE){

  ### loop for `l_varsy` & `l_varsx`
    if(arg_multivar == 'univariate'){
      # univariate
      for(j in l_varsy) {
      for(i in l_varsx) {
         forms <- paste0(j, ' ~ ', i)

        l_formulas[[length(l_formulas) + 1]] = forms
      }

      }

    } else if (arg_multivar == 'multivariate'){
      # multivariate
        for(j in l_varsy) {
        for(i in l_varsx) {
          forms <- paste0(j, ' ~ ', paste0(unlist(l_varsx, i), collapse = ' + '))
        }
       l_formulas[[length(l_formulas) + 1]] = forms
      }
    }

phen_1 <- phen %>%
     dplyr::select(l_group, l_varsy, l_varsx) %>%
     group_by_at(l_group) %>%
     group_nest() %>%
     slice(rep(1:n(), each = length(l_formulas))) %>%
     arrange_at(all_of(l_group)) %>%
     mutate(l_forms = unlist(rep(l_formulas, times = dim(unique(phen[, l_group]))[1])))


  }

    else if(custom_formls == TRUE){

        l_formulas <- formls

        phen_1 <- phen %>%
             dplyr::select(l_group, l_varsy, l_varsx) %>%
         arrange_at(l_group) %>%
                  group_by(ig) %>% #He canviat 'isotype' per 'ig'
                  group_nest() %>%
                  mutate(l_forms = rep(l_formulas)) %>% #first we put the corresponding formula per isotype
                  unnest(cols = data) %>%
                  group_by_at(l_group) %>% #and finally we arrange the data so that it looks like when we
                  #don't have customized formulas
                  group_nest() %>%
                  mutate(l_form = map(data, ~list(.x$l_forms)),
                         l_forms = map(l_form, ~ .x[[1]][[1]])) %>%
                  dplyr::select(-l_form)

    }

### extract whole summary of 'lm' models

      phen_1 <- phen_1 %>%
                mutate(l_mods = map2(data, l_forms,
                        ~ try(lm(formula = as.formula(.y), data = .x, model = TRUE))),
                       smry_model = map(l_mods, ~ janitor::clean_names(broom::glance(.x), case = 'snake')),
                       smry_beta = map(l_mods, ~ janitor::clean_names(broom::tidy(.x), case = 'snake')),
                       ci_lower = map(l_mods, ~ tibble::as_tibble(confint(.x, level = .95))[, 1]),
                       ci_upper = map(l_mods, ~ tibble::as_tibble(confint(.x, level = .95))[, 2])) %>%
                dplyr::select(-data, -l_mods) %>%
                unnest(cols = c('l_forms', 'smry_model', 'smry_beta', 'ci_lower', 'ci_upper'),
                       names_repair = ~ tidyr_legacy(., sep = '_')) %>% ## '_1' fix 'data_beta' info for 'statistic' & 'p_value'
                filter(!str_detect(term, 'Intercept')) %>%
                dplyr::select(everything(),
                       t_value = statistic,
                       beta = estimate, vars_x = term, t_value_beta = statistic_1, raw_pval = p_value_1,
                       ci_lower = starts_with('2.5') , ci_upper = starts_with('97.5')) %>%
                separate(l_forms, c('vars_y', 'rest'), sep = ' ~ ', remove = TRUE) %>%
                mutate_at(vars(c('vars_y', 'vars_x')), as.factor) %>%
                dplyr::select(-rest)


### argument 'arg_filter_vars'
  if(!missing(arg_filter_vars)) {
  # update 'data'
    phen_1 <- phen_1 %>%
      filter(str_detect(vars_x, pattern = paste(arg_filter_vars, collapse = '|')))
  }

# list of variables to choose
l_vars <- c(l_group, 'vars_y', 'vars_x', 'beta', 'ci_lower', 'ci_upper', 'raw_pval')
l_vars_arg <- c(l_vars, arg_info_model, arg_info_beta)

### update 'data'
  phen_1 <- phen_1 %>%
    dplyr::select(one_of(l_vars_arg))


### argument 'format_pval'
if(format_pval == TRUE & !missing(padj_method)) {

        #if `gr_adj` is not missing, the adjust will be made grouping by such variables
        if(!missing(gr_adj)){
        phen_1 <- phen_1 %>%
                  group_by_at(gr_adj) %>%
                  group_nest() %>%
                  mutate(adj_data = map(data,
                      ~ chmi.stat.p_adjust(phen = .x, method = padj_method, arrange_by = arg_pval_arrange))) %>%
                  dplyr::select(-data) %>%
                  unnest(cols = adj_data) %>%
                  arrange(raw_pval)

        }else{
        phen_1 <- chmi.stat.p_adjust(phen = phen_1, method = padj_method, arrange_by = arg_pval_arrange) %>%
             mutate_at(vars(c('adj_pval_text', 'adj_pval_signif')), as.factor)
                  }
  }

  # return
  return(phen_1)
}



#' @export
new_covidcat.beta_transform <- function(mod_results_phen)
{

mod_results_phen <- mod_results_phen %>%
                    #this below is all for log-linear interpretation (both numeric and categorical variables)
                    mutate(beta_transform_unit = 10^beta,
                           beta_transform_percent = (beta_transform_unit-1)*100,
                           ci_lower_transform_unit = 10^ci_lower,
                           ci_lower_transform_percent = (ci_lower_transform_unit-1)*100,
                           ci_upper_transform_unit = 10^ci_upper,
                           ci_upper_transform_percent = (ci_upper_transform_unit-1)*100)

### return
  return(mod_results_phen)
}

# Set your own working directory:
setwd("")

#comes from the redcap database cleaned by Julia Montanya:
dta<-read.csv("Cizur_data_20200626_clean_final0804.csv")

#Ig results from the lab (original archive, contains the NC):
all_lab<-read_xlsx("Serologias Encuesta Navarra FINAL.xlsx") 

#Comes from the lab database put into wide in excel and then used in STATA:
lab<-read_dta("20200810_lab_results.dta")

    lab$Global_posneg<- replace(lab$Global_posneg,which(lab$Global_posneg=="indet"), "pos") #Cause in the end they were not indets
    #table(lab$Global_posneg)
    # neg pos 
    # 673  56 

    lab<-lab%>%
    rename(RBD_IgG_posneg=RBG_IgG_posneg,
          RBD_IgM_posneg=RBG_IgM_posneg,
          RBD_IgA_posneg=RBG_IgA_posneg,
          Spike_IgG=SpikeIgG)

    lab$S_IgM_posneg<-as.factor(lab$S_IgM_posneg)
    lab$RBD_IgM_posneg<-as.factor(lab$RBD_IgM_posneg)
    lab$IgM_posneg<-as.factor(lab$IgM_posneg)
    lab$S_IgG_posneg<-as.factor(lab$S_IgG_posneg)
    lab$RBD_IgG_posneg<-as.factor(lab$RBD_IgG_posneg)
    lab$IgG_posneg<-as.factor(lab$IgG_posneg)
    lab$S_IgA_posneg<-as.factor(lab$S_IgA_posneg)
    lab$RBD_IgA_posneg<-as.factor(lab$RBD_IgA_posneg)
    lab$IgA_posneg<-as.factor(lab$IgA_posneg)
    lab$Global_posneg<-as.factor(lab$Global_posneg)

    #Change to numeric (first change decimal commas to points)
    lab$Spike_IgM<- sub(",",".",lab$Spike_IgM)
    lab$Spike_IgM<-as.numeric((lab$Spike_IgM),digits=10)
    lab$RBD_IgM<- sub(",",".",lab$RBD_IgM)
    lab$RBD_IgM<-as.numeric((lab$RBD_IgM),digits=10)
    lab$Spike_IgG<- sub(",",".",lab$Spike_IgG)
    lab$Spike_IgG<-as.numeric((lab$Spike_IgG),digits=10)
    lab$RBD_IgG<- sub(",",".",lab$RBD_IgG)
    lab$RBD_IgG<-as.numeric((lab$RBD_IgG),digits=10)
    lab$Spike_IgA<- sub(",",".",lab$Spike_IgA)
    lab$Spike_IgA<-as.numeric((lab$Spike_IgA),digits=10)
    lab$RBD_IgA<- sub(",",".",lab$RBD_IgA)
    lab$RBD_IgA<-as.numeric((lab$RBD_IgA),digits=10)
   
#Controls (original archive):
controls<-read_xls("Controls negatius bombers_total.xls")

#Merge redcap data and lab data for general use
ddbb<-merge(dta,lab, by = "lab_sample_grouplab_sample_id", all=TRUE)

#Many cells with empty spaces, we transform them into NAs
ddbb <- ddbb %>%
  mutate(across(where(is.character), ~ na_if(.,""))) #Only in character vars

#For the variable (socioecon_groupeducation) some options were lacking, so I create the PhD variable and "School" and re-atribute some that were manually written on "other" but only for those older than 20 y/o
ddbb$socioecon_groupeducation[ddbb$studyn==142 | ddbb$studyn==194 | ddbb$studyn==277 | ddbb$studyn==305 | ddbb$studyn==326 | ddbb$studyn==333 | ddbb$studyn==344 | ddbb$studyn==375 | ddbb$studyn==381 | ddbb$studyn==436 | ddbb$studyn==623 | ddbb$studyn==684 | ddbb$studyn==755 | ddbb$studyn==759 | ddbb$studyn==833 | ddbb$studyn==853 | ddbb$studyn==856 |ddbb$studyn==899 | ddbb$studyn==983 | ddbb$studyn==1070]<-"PhD"
ddbb$socioecon_groupeducation[ddbb$studyn==49 | ddbb$studyn==289 | ddbb$studyn==335 | ddbb$studyn==376]<-"university_degree"
ddbb$socioecon_groupeducation[ddbb$studyn==33 | ddbb$studyn==418 | ddbb$studyn==431 | ddbb$studyn==533 | ddbb$studyn==573 | ddbb$studyn==574  | ddbb$studyn==590 | ddbb$studyn==593 | ddbb$studyn==724 | ddbb$studyn==748 | ddbb$studyn==766 | ddbb$studyn==847 | ddbb$studyn==1035 | ddbb$studyn==1140]<-"secondary_education"
ddbb$socioecon_groupeducation[ddbb$studyn==614]<-"technical_formacion_profesional"
ddbb$socioecon_groupeducation[ddbb$studyn==205 | ddbb$studyn==419 | ddbb$studyn==994 | ddbb$studyn==1198]<-"school"

#Add variables to Put into factor:
cols<-c('participantsex', 'ageg','puebloconcejo','symptom_groupcovid_yn','symptom_groupfever_37','symptom_groupfever_38','symptom_groupchills','symptom_groupfatigue','symptom_groupmuscle_ache','symptom_groupsore_throat','symptom_grouplost_smell','symptom_grouplost_taste','symptom_groupcough','symptom_grouprunny_nose','symptom_groupshort_breath','symptom_groupdyspnea','symptom_groupwheeze',
        'symptom_groupchest_pain','symptom_groupother_resp','symptom_groupheadache','symptom_groupnausea','symptom_groupab_pain','symptom_groupdairrhea','symptom_groupexpectorations','symptom_questions_groupsymptom_m','symptom_questions_groupsymptom_w','symptom_questions_groupsymptom_h','symptom_questions_groupsymptom_v','other_hh_covid_confirmed_groupot',
        'other_hh_covid_suspected_groupot','risk_factors_groupsuffering_dise','risk_factors_grouptaking_medicat','risk_factors_groupblood_group','risk_factors_groupbcg_vac','risk_factors_groupflu_vac','risk_factors_groupflu_yn', 'risk_factors_groupcold_yn','risk_factors_groupsmoker_yn','socioecon_groupeducation','socioecon_groupeducation_other','socioecon_groupoccupation',
        'socioecon_groupemployed_yn','socioecon_groupemployed_before_y','socioecon_groupjob_loss_covid','socioecon_groupsalary_mode','socioecon_groupcovid_affect_fina','socioecon_groupstruggle_expenses','exposures_groupsanitary_centre','exposures_groupsanitary_centre_r','exposures_groupsanitary_accompan','exposures_groupsanitary_symptoms','exposures_groupsanitary_incompat',
        'exposures_groupdiagnosed_covid','exposures_groupplaces_visited','exposures_groupsupermarket_pre','exposures_groupgrocery_pre','exposures_groupbank_pre','exposures_groupchurch_pre','exposures_grouphairdresser_pre','exposures_grouppamplona_pre','exposures_groupmadrid_pre','exposures_groupother_spain_pre','exposures_groupoutside_spain_pre','exposures_groupplaces_visited_su',
        'exposures_groupcontact_covid','exposures_groupcontact_covid_rel','exposures_groupcontact_covid_pro','perceptions_groupinfo_source','perceptions_groupcovid_transmit','perceptions_groupcovid_prevent','perceptions_groupcovid_info_rece')
ddbb[cols]<-lapply(ddbb[cols],factor)

#Add variables to Put into date format:
cols<-c('participantdob','symptom_groupcovid_date_diagnose','symptom_groupfever_37__start','symptom_groupfever_37_end','symptom_groupfever_38__start','symptom_groupfever_38_end', 'symptom_groupchills_start','symptom_groupchills_end','symptom_groupfatigue_start','symptom_groupfatigue_end','symptom_groupmuscle_ache_start','symptom_groupmuscle_ache_end',
       'symptom_groupsore_throat_start','symptom_groupsore_throat_end','symptom_grouplost_smell_start','symptom_grouplost_smell_end','symptom_grouplost_taste_start','symptom_grouplost_taste_end','symptom_groupcough_start','symptom_groupcough_end','symptom_grouprunny_nose_start','symptom_grouprunny_nose_end','symptom_groupshort_breath_start','symptom_groupshort_breath_end',
       'symptom_groupdyspnea_start','symptom_groupdyspnea_end','symptom_groupwheeze_start','symptom_groupwheeze_end','symptom_groupchest_pain_start','symptom_groupchest_pain_end','symptom_groupother_resp_start','symptom_groupother_resp_end','symptom_groupheadache_start','symptom_groupheadache_end','symptom_groupnausea_start','symptom_groupnausea_end',
       'symptom_groupexpectorations_star','symptom_groupexpectorations_end','symptom_groupdairrhea_start','symptom_groupdairrhea_end','risk_factors_groupflu_end_date','risk_factors_groupcold_end_date','socioecon_groupwork_last_day','exposures_groupsanitary_work_dat','exposures_groupsanitary_other_da','exposures_groupsupermarket_date','exposures_groupgrocery_date','exposures_groupbank_date',
       'exposures_groupchurch_date','exposures_grouphairdresser_date','exposures_grouppamplona_date','exposures_groupmadrid_date','exposures_groupother_spain_date','exposures_groupoutside_spain_dat','exposures_groupcontact_covid_dat','symptom_groupab_pain_start')
ddbb[cols]<-(lapply(ddbb[cols],as.Date,"%Y-%m-%d"))

#symptom_questions_groupsymptom_m = Did any of these symptoms require you to seek medical attention?
#symptom_questions_groupsymptom_w = Did any of these symptoms require you to miss work or school?
#symptom_questions_groupsymptom_h = Did any of these symptoms require you to be hospitalized?
#symptom_questions_groupsymptom_v = Did any of these symptoms require you to get ventilation?

#I create the variable S_pos:
ddbb<-ddbb%>%
mutate(S_pos=as.factor(ifelse((S_IgM_posneg=="pos"|S_IgG_posneg=="pos"|S_IgA_posneg=="pos"),1,0)))

#I create the variable RBD_pos:
ddbb<-ddbb%>%
mutate(RBD_pos=as.factor(ifelse((RBD_IgM_posneg=="pos"|RBD_IgG_posneg=="pos"|RBD_IgA_posneg=="pos"),1,0)))

# I create the variable body mass index:
ddbb<-ddbb%>%
 mutate(bmi=(risk_factors_groupweight/((risk_factors_groupheight)*(risk_factors_groupheight))))

# I create a variable for allergy y/n:
ddbb<-ddbb%>%
 mutate(allergy=as.factor(ifelse(grepl("allergy",risk_factors_groupsuffering_dise), "yes", "no")))

# I create a var  for autoinmune disease, manually I detect two among the seropositive: 22 has celiac disease and 968 has ulcerative colitis.
ddbb<-ddbb%>%
 mutate(autoimm_dis=as.factor(ifelse(studyn==22 | studyn==968, "yes", "no")))

# I create the variable "ratio_rbd_s" by dividing the levels of IgG-RBD by IgG-S
ddbb<-ddbb%>%
 mutate(IgG_ratio_rbd_s=(RBD_IgG/Spike_IgG),
       IgM_ratio_rbd_s=(RBD_IgM/Spike_IgM),
       IgA_ratio_rbd_s=(RBD_IgA/Spike_IgA))

# I create the variable "IgA_RBD_S" by adding up the levels of IgA-RBD by IgA-S
ddbb<-ddbb%>%
 mutate(IgA_RBD_S=(RBD_IgA+Spike_IgA))

# I create the variable "IgM_RBD_S" by adding up the levels of IgA-RBD by IgA-S
ddbb<-ddbb%>%
 mutate(IgM_RBD_S=(RBD_IgM+Spike_IgM))

# I create the variable "IgG_RBD_S" by adding up the levels of IgA-RBD by IgA-S
ddbb<-ddbb%>%
 mutate(IgG_RBD_S=(RBD_IgG+Spike_IgG))

#I Create a var time since diagnosed
ddbb<-ddbb%>%
mutate(days_since_diag= difftime(lab_sample_grouplab_sample_date ,symptom_groupcovid_date_diagnose, units="days"))

#I create a var with start of symptoms
start<-ddbb%>%
dplyr::select("studyn" | (starts_with("symptom_") & (ends_with("start") | ends_with("star"))))%>%
mutate(symp_start=pmin( symptom_groupfever_37__start,symptom_groupfever_38__start,symptom_groupchills_start,symptom_groupfatigue_start,symptom_groupmuscle_ache_start,symptom_groupsore_throat_start,symptom_grouplost_smell_start,symptom_grouplost_taste_start,symptom_groupcough_start,symptom_grouprunny_nose_start,symptom_groupshort_breath_start,symptom_groupdyspnea_start,symptom_groupwheeze_start,symptom_groupchest_pain_start,symptom_groupother_resp_start,symptom_groupheadache_start,symptom_groupnausea_start,symptom_groupab_pain_start,symptom_groupdairrhea_start,symptom_groupexpectorations_star,na.rm=TRUE))%>%
dplyr::select(studyn,symp_start)
start$symp_start[start$symp_start < as.Date("2020-02-01")] <- NA #We will consider symptoms earlier than 1st of feb as not COVID-19
ddbb<-merge(ddbb,start, by = "studyn", all=TRUE)

# Days pso
ddbb<-ddbb%>%
mutate(days_pso= difftime(lab_sample_grouplab_sample_date ,symp_start, units="days"))

# Group symptoms to have less variables on the models
ddbb<-ddbb%>%
mutate(fever=as.factor(ifelse((symptom_groupfever_37=="yes" | symptom_groupfever_38=="yes" | symptom_groupchills=="yes"), "yes", "no")))%>%
mutate(loss_taste_smell=as.factor(ifelse((symptom_grouplost_smell=="yes" | symptom_grouplost_taste=="yes"), "yes", "no")))%>%
mutate(lower_resp=as.factor(ifelse((symptom_groupcough=="yes" | symptom_groupshort_breath=="yes"| symptom_groupdyspnea=="yes"| symptom_groupexpectorations=="yes"| symptom_groupwheeze=="yes"| symptom_groupchest_pain=="yes"), "yes", "no")))%>%
mutate(upper_resp=as.factor(ifelse((symptom_groupsore_throat=="yes" | symptom_grouprunny_nose=="yes"), "yes", "no")))%>%
mutate(dige=as.factor(ifelse((symptom_groupnausea=="yes" | symptom_groupab_pain=="yes"| symptom_groupdairrhea=="yes"), "yes", "no")))%>%
mutate(ache_fatigue=as.factor(ifelse((symptom_groupfatigue=="yes" | symptom_groupmuscle_ache=="yes"| symptom_groupheadache=="yes"), "yes", "no")))

#symp yes/no
ddbb<-ddbb%>%
mutate(symp_yn=as.factor(ifelse((symptom_groupfever_37=="yes" | symptom_groupfever_38=="yes" | symptom_groupchills=="yes" |symptom_grouplost_smell=="yes" | symptom_grouplost_taste=="yes" |symptom_groupcough=="yes" | symptom_groupshort_breath=="yes"| symptom_groupdyspnea=="yes"| symptom_groupexpectorations=="yes"| symptom_groupwheeze=="yes"| symptom_groupchest_pain=="yes"|symptom_groupsore_throat=="yes" | symptom_grouprunny_nose=="yes"|symptom_groupnausea=="yes" | symptom_groupab_pain=="yes"| symptom_groupdairrhea=="yes"|symptom_groupfatigue=="yes" | symptom_groupmuscle_ache=="yes"| symptom_groupheadache=="yes"), "yes", "no")))
ddbb$symp_yn[is.na(ddbb$symp_yn)]<-"no" #change NAs to no


#Comorb
#symp yes/no
ddbb<-ddbb%>%
mutate(comorb_yn=as.factor(ifelse(!is.na(risk_factors_groupsuffering_dise), "yes", "no")))
    
#exposures y/n
ddbb<-ddbb%>%
mutate(exp_yn=as.factor(ifelse(!is.na(exposures_groupplaces_visited), "yes", "no")))

#schooling level
ddbb<-ddbb%>%
mutate(schooling=socioecon_groupeducation)
ddbb$schooling[ddbb$schooling=="masters" | ddbb$schooling=="PhD" ]<-"university_degree" #remove two categories
levels(ddbb$schooling)<-list("University degree / Master / PhD"='university_degree',"Formación profesional"="technical_formacion_profesional","Secondary education"='secondary_education',"Primary education"="school")

str(ddbb, list.len=ncol(ddbb))
options(repr.matrix.max.cols=300, repr.matrix.max.rows=734) # I set it like this so that it will display my whole data set

#The ones I have to gather into long (the ones that are different for each isotype):
    #'Spike_IgM'
    #'RBD_IgM'
    #'Spike_IgG'
    #'RBD_IgG'
    #'Spike_IgA'
    #'RBD_IgA'
    #'IgA_RBD_S'
    #'IgM_RBD_S'
    #'IgG_RBD_S'
    
    #'S_IgM_posneg'
    #'RBD_IgM_posneg'
    #'S_IgG_posneg'
    #'RBD_IgG_posneg'
    #'S_IgA_posneg'
    #'RBD_IgA_posneg'

    #IgM_posneg
    #IgG_posneg
    #IgA_posneg
    #S_pos
    #RBD_pos
    
    #'ratio_rbd_s' --> IgG 

clean<-ddbb%>%
 filter(lab_sample_grouplab_sample_id!=952)%>% #we delete those ID's which are not complete
 filter(lab_sample_grouplab_sample_id!=999)%>%
 filter(studyn!=319)%>%
 filter(studyn!=1077)%>%
 filter(studyn!=325)%>%
 filter(studyn!=688)%>%
 filter(studyn!=981)
#length(unique(clean$studyn)) #727

    #First columns that are Ig-Ag-studyn specific
ddbb_long0<-clean%>%
    dplyr::select('studyn','Spike_IgM','RBD_IgM',
              'Spike_IgG','RBD_IgG',
              'Spike_IgA','RBD_IgA',
              'IgA_RBD_S','IgM_RBD_S','IgG_RBD_S')%>%
    rename(S_IgA=Spike_IgA)%>%
    rename(S_IgM=Spike_IgM)%>%
    rename(S_IgG=Spike_IgG)%>%
    rename(RBDS_IgA=IgA_RBD_S)%>%
    rename(RBDS_IgM=IgM_RBD_S)%>%
    rename(RBDS_IgG=IgG_RBD_S)%>%
    gather(ig_ag, mfi, S_IgM:RBDS_IgG)%>%
    separate(ig_ag, c("ag", "ig"), sep="_")%>%unique()
#summary(ddbb_long0)
#head(arrange(ddbb_long0, studyn))
#length(ddbb_long0$studyn)
#length(unique(ddbb_long0$studyn))

ddbb_long1<-clean%>%
    dplyr::select('studyn','S_IgM_posneg','RBD_IgM_posneg',
              'S_IgG_posneg','RBD_IgG_posneg',
              'S_IgA_posneg','RBD_IgA_posneg')%>%
    gather(ig_ag, ig_ag_posneg, S_IgM_posneg:RBD_IgA_posneg)%>%
    separate(ig_ag, c("ag", "ig"), sep="_")%>%unique()
#summary(ddbb_long1)
#head(arrange(ddbb_long1, studyn))
#length(ddbb_long1$studyn)
#length(unique(ddbb_long1$studyn))

ddbb_long<-merge(ddbb_long0, ddbb_long1, by=c('studyn','ag','ig'), all=TRUE)%>%unique
#length(ddbb_long$studyn) 
#head(arrange(ddbb_long, studyn))
#length(unique(ddbb_long$studyn))

    #Add vars specific to studyn-Ag
ddbb_long2<-clean%>%
    dplyr::select('studyn','S_pos', 'RBD_pos')%>%
    gather(ig_ag, ag_posneg, S_pos:RBD_pos)%>%
    separate(ig_ag, c("ag", "trash"), sep="_")%>%
    dplyr::select(-trash)%>%unique()
#table(ddbb_long2$ag)
#table(ddbb_long2$ag_posneg)

ddbb_long<-merge(ddbb_long, ddbb_long2, by=c('studyn','ag'), all=TRUE)%>%unique()


    #Add vars specific to studyn-Ig
ddbb_long3<-clean%>%
    dplyr::select('studyn','IgM_posneg', 'IgG_posneg', 'IgA_posneg')%>%
    gather(ig_ag, ig_posneg, IgM_posneg:IgA_posneg)%>%
    separate(ig_ag, c("ig", "trash"), sep="_")%>%
    dplyr::select(-trash)%>%unique()

ddbb_long<-merge(ddbb_long, ddbb_long3, by=c('studyn','ig'), all=TRUE)%>%unique()

ddbb_long4<-clean%>%
    dplyr::select('studyn','IgM_ratio_rbd_s', 'IgA_ratio_rbd_s', 'IgG_ratio_rbd_s')%>%
    gather(ig_ag, ratio_mfi, IgM_ratio_rbd_s:IgG_ratio_rbd_s)%>%
    separate(ig_ag, c("ig", "trash"), sep="_")%>%
    dplyr::select(-trash)%>%unique()

ddbb_long<-merge(ddbb_long, ddbb_long4, by=c('studyn','ig'), all=TRUE)%>%unique()

#summary(ddbb_long)
#head(ddbb_long)
#table(ddbb_long$ig)
#table(ddbb_long$ig_posneg)

    #Finally we bind all columns that are only studyn-specific
rest<-ddbb%>%
dplyr::select(-'Spike_IgM',-'RBD_IgM',
              -'Spike_IgG',-'RBD_IgG',
             - 'Spike_IgA',-'RBD_IgA',
              -'IgA_RBD_S',-'IgM_RBD_S',-'IgG_RBD_S',
             -'S_IgM_posneg',-'RBD_IgM_posneg',
              -'S_IgG_posneg',-'RBD_IgG_posneg',
              -'S_IgA_posneg',-'RBD_IgA_posneg',
             -'S_pos', -'RBD_pos',
             -'IgM_posneg', -'IgG_posneg', -'IgA_posneg',
             -'IgM_ratio_rbd_s',- 'IgA_ratio_rbd_s',- 'IgG_ratio_rbd_s')

ddbb_long<-merge(rest,ddbb_long, by=c('studyn'), all=TRUE)%>%unique()
length(unique(ddbb_long$studyn)) #733

#Add variables to Put into factor:
cols<-c('ig','ag','ig_ag_posneg','ag_posneg','ig_posneg')
ddbb_long[cols]<-lapply(ddbb_long[cols],factor)

#summary(ddbb_long)
#head(ddbb_long)

#Tune a bit the dataframe so that labels are nicer in the plot
fig<-ddbb%>% rename (Seropositivity=Global_posneg)

figure<-ggplot(subset(fig, !is.na(Seropositivity)), aes(x=symp_start, fill=Seropositivity)) + 
  geom_histogram(binwidth=7, colour="white")+
  theme_bw()+
  scale_x_date(labels = date_format("%d / %b"),
                 date_breaks = "1 week",
                 limits = c(as.Date("2020-01-15"), as.Date("2020-05-30"))) +
  ylab("Frequency of reported COVID-19 compatible symptoms") + xlab("Date (2020)") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  theme(legend.position="top")
figure
#ggsave(figure, filename = "Dates of symp start.png", width = 7, height= 5)

#summary(ddbb%>%dplyr::select("studyn",'Spike_IgM','S_IgM_posneg','RBD_IgM','RBD_IgM_posneg','IgM_posneg','Spike_IgG','S_IgG_posneg','RBD_IgG','RBD_IgG_posneg','IgG_posneg','Spike_IgA','S_IgA_posneg','RBD_IgA','RBD_IgA_posneg','IgA_posneg','Global_posneg'))
venn<-ddbb%>%dplyr::select("lab_sample_grouplab_sample_id","studyn",'Spike_IgM','S_IgM_posneg','RBD_IgM','RBD_IgM_posneg','IgM_posneg','Spike_IgG','S_IgG_posneg','RBD_IgG','RBD_IgG_posneg','IgG_posneg','Spike_IgA','S_IgA_posneg','RBD_IgA','RBD_IgA_posneg','IgA_posneg','Global_posneg')

length(unique(venn$studyn))

venn.table<-ddbb%>%
 filter(studyn!=319)%>%
 filter(studyn!=1077)%>%
 filter(studyn!=325)%>%
 filter(studyn!=688)%>%
 filter(studyn!=981)%>%
 dplyr::select(S_IgM_posneg, RBD_IgM_posneg,S_IgG_posneg, RBD_IgG_posneg,S_IgA_posneg, RBD_IgA_posneg)

summary(venn.table)

st(venn.table, out="kable") #I paste the code in a new Marckdown cell
#st(regr, out="kable") #to save it in a csv

a<-ddbb%>%filter(IgM_posneg=="pos")
length(unique(a$studyn))
a<-ddbb%>%filter(IgG_posneg=="pos")
length(unique(a$studyn))
a<-ddbb%>%filter(IgA_posneg=="pos")
length(unique(a$studyn))

a<-ddbb%>%filter(IgM_posneg=="pos" & IgG_posneg=="pos")
length(unique(a$studyn))
a<-ddbb%>%filter(IgM_posneg=="pos" & IgA_posneg=="pos")
length(unique(a$studyn))
a<-ddbb%>%filter(IgG_posneg=="pos" & IgA_posneg=="pos")
length(unique(a$studyn))

a<-ddbb%>%filter(IgM_posneg=="pos" & IgG_posneg=="pos" & IgA_posneg=="pos")
length(unique(a$studyn))

venn.table<-venn.table%>%
 mutate(S_pos=ifelse((S_IgM_posneg=="pos"|S_IgG_posneg=="pos"|S_IgA_posneg=="pos"),"pos","neg"))%>%
 mutate(RBD_pos=ifelse((RBD_IgM_posneg=="pos"|RBD_IgG_posneg=="pos"|RBD_IgA_posneg=="pos"),"pos","neg"))

venn.table<-venn.table%>%
 mutate(ag_posneg=ifelse(S_pos=="pos"& RBD_pos=="neg", "Only S",
                          ifelse(S_pos=="neg" & RBD_pos=="pos","Only RBD",
                                ifelse (S_pos=="pos" & RBD_pos=="pos", "Both pos", "neg"))))
venn.table$S_pos<-as.factor(venn.table$S_pos)
venn.table$RBD_pos<-as.factor(venn.table$RBD_pos)
venn.table$ag_posneg<-as.factor(venn.table$ag_posneg)
summary(venn.table)

#Now I get a dataframe from the big one with the variables I want to put in the model
regr<-ddbb%>% # we need one column for ag and one for ig
 filter(lab_sample_grouplab_sample_id!=952)%>% #we delete those ID's which are not complete
 filter(lab_sample_grouplab_sample_id!=999)%>%
 filter(studyn!=319)%>%
 filter(studyn!=1077)%>%
 filter(studyn!=325)%>%
 filter(studyn!=688)%>%
 filter(studyn!=981)%>%
 dplyr::select( "studyn","Global_posneg","participantsex", "ageg","puebloconcejo","bmi", "risk_factors_groupblood_group","risk_factors_groupbcg_vac","risk_factors_groupflu_vac","symptom_groupcovid_yn",'symp_yn','comorb_yn',"other_hh_covid_suspected_groupot","socioecon_grouphh_size",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","allergy", "schooling", "socioecon_groupoccupation", "exp_yn")%>% unique() %>% droplevels()

regr[regr=="unknown"]<-NA #we don't want unknown to be considered a category for the models
regr<-droplevels(regr) #to drop the unknown levels

regr$Global_posneg<- as.character(regr$Global_posneg) #we need 1 and 0, so I convert to character first
regr$Global_posneg<-replace(regr$Global_posneg, which(regr$Global_posneg=="pos"),1)
regr$Global_posneg[regr$Global_posneg=="neg"]<-"0"
regr$Global_posneg<- as.numeric(regr$Global_posneg) #convert to numeric (seen in: https://www.youtube.com/watch?v=WnmwuD8OwMw)

summary(regr)
length(unique(regr$studyn))

univar_outcome <- lapply(c("participantsex", "ageg","puebloconcejo","bmi", "risk_factors_groupblood_group","risk_factors_groupbcg_vac","risk_factors_groupflu_vac","symptom_groupcovid_yn",'symp_yn','comorb_yn',"socioecon_grouphh_size",
               "risk_factors_groupcold_yn","allergy", "schooling"), function(var) {
                           
                           formula    <- as.formula(paste("Global_posneg~", var))
                           res.logist <- glm(formula, data = regr, family = binomial) #glm/lm depenent de l'outcome
                           
                           results_univar_outcome <-summary(res.logist)$coefficients
                           #tab_model(res.logist) #If we want the ORs instead of the coefficients in an html format
                         })
univar_outcome
#Save the univariate in txt: change name, paste in a csv and separate by space, then re-name the columns:
#capture.output(univar_outcome, file = "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/univar_seropos.csv")  

#Another way to transform the coefficients (betas) into 
#Define the vars to be compaired against GLobal_posneg:
regr.risk<-regr%>%dplyr::select("participantsex", "ageg","puebloconcejo", "risk_factors_groupblood_group","risk_factors_groupbcg_vac","risk_factors_groupflu_vac","symptom_groupcovid_yn",'symp_yn','comorb_yn',
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","allergy", "schooling", "exp_yn")
#Define the function for comparation, in this case the Wald RR:
rr<-function(vars){
    riskratio(table(vars,regr$Global_posneg))
}
#Apply it:
lapply(regr.risk, FUN=rr)

#regr<-regr%>%filter(!is.na(lower_resp)) #have to delete 
length(unique(regr$studyn)) #727
model.outcome = glm(Global_posneg ~ ageg+risk_factors_groupcold_yn+symp_yn+symptom_groupcovid_yn+socioecon_grouphh_size, #glm/lm depenent de l'outcome
                    data=regr,
                    family = binomial)
summary(model.outcome)

# Now stepwise regression
step.model.outcome <- stepAIC(model.outcome)
coef(step.model.outcome)

model.outcome.final = glm(Global_posneg ~ ageg+risk_factors_groupcold_yn+symp_yn+symptom_groupcovid_yn+socioecon_grouphh_size,
                          data=regr,
                          family = binomial)
summary(model.outcome.final)
tab_model(model.outcome.final)

#Now I get a dataframe from the big one with the variables I want to put in the model
regr<-ddbb_long%>% # we need one column for ag and one for ig
 filter(lab_sample_grouplab_sample_id!=952)%>% #we delete those ID's which are not complete
 filter(lab_sample_grouplab_sample_id!=999)%>%
 filter(studyn!=319)%>%
 filter(studyn!=1077)%>%
 filter(studyn!=325)%>%
 filter(studyn!=688)%>%
 filter(studyn!=981)%>%
 filter(Global_posneg=="pos" & !is.na(days_pso))%>% # because we want to see what determines being  positive for S once you're infected and only those with pso info
 dplyr::select( "studyn","ig","ag", "mfi", "days_pso","participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy")%>% unique() %>% droplevels()
#Get rid of uncompensated vars: symptom_questions_groupsymptom_v,symptom_questions_groupsymptom_h,risk_factors_groupflu_yn because all or too many answeres were "no", autoinmune disease
regr$days_pso<-as.numeric(regr$days_pso)

#Put in log10 the variables that are continuous and do not follow a normal distribution
shapiro.test(regr$mfi)#--> test variable per variable

#Not normal, then we convert to log10:
col_logs<-c("mfi")
regr[col_logs] <- log10(regr[col_logs]) 
regr<-regr%>% rename(logmfi=mfi)

length(unique(regr$studyn)) #36 --> correct, is the number of positive participants for any of the ig-ag pairs and with pso
dim(regr) #108 (36 participants * 3 ig)
summary(regr)

#UNIVARIAT 
#Predictor or "x" variables:
l_varsx<-c("days_pso","participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy")
mfi<-c("logmfi") #log10(mfi)
ig_ag<-c("ig", "ag")
gr_adj<-c('ig', 'vars_x')
padj_method<-"BH" #verificar
arg_pval_arrange<-"raw_pval"

tab_uni_logmfi <- chmi.stat.lm_beta(
                          phen = regr, #change for db name
                          l_group = ig_ag, 
                          l_varsy = mfi, 
                          l_varsx = l_varsx, 
                          arg_multivar = 'univariate', #can be 'multivariate'
                          arg_info_model = c('aic', 'bic', 'adj_r_squared'),
                          arg_info_beta = c('std_error'),
                          gr_adj = gr_adj, #not change
                          format_pval = TRUE, #ajustat
                          arg_pval_arrange = arg_pval_arrange, #primer mirem el raw p-val
                          padj_method = padj_method) %>% #tipus de adjust de multiple comparisons
                          new_covidcat.beta_transform() %>% # es una funcio creada i carregada abans ( si fos log-log seria diferent)
                          mutate_at(vars(beta, ci_lower, ci_upper, raw_pval, aic, bic, adj_r_squared, std_error, adj_pval,
                              beta_transform_unit, beta_transform_percent, ci_lower_transform_unit, ci_lower_transform_percent, ci_upper_transform_unit, ci_upper_transform_percent), ~ round(.,4)) #round at 4 decimal

#Export to csv:
#write.csv (tab_uni_logmfi, "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/Tables/raw_univar_logmfi.csv")

    #IgA
tab_uni_logmfi_iga<-tab_uni_logmfi%>%
filter(ig=="IgA")
tab_uni_logmfi_iga<-arrange(tab_uni_logmfi_iga, ag)
#write.csv (tab_uni_logmfi_iga, "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/Tables/raw_univar_logmfi_iga.csv")

#p-value <0.05 for IgA - RBD
#loss_taste_smell
#p-value <0.05 for IgA - S
#allergy + fever 
#p-value <0.05 for IgA - RBD+S
#allergy + fever + loss_taste_smell


    #IgG
tab_uni_logmfi_igg<-tab_uni_logmfi%>%
filter(ig=="IgG")
tab_uni_logmfi_igg<-arrange(tab_uni_logmfi_igg, ag)
#write.csv (tab_uni_logmfi_igg, "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/Tables/raw_univar_logmfi_igg.csv")

#p-value <0.05 for IgG - RBD
#allergy + risk_factors_groupcold_yn +  upper_resp + participantsex + fever + age
#p-value <0.05 for IgG - S
#allergy + participantsex + loss_taste_smell
#p-value <0.05 for IgG - RBD+S
#allergy + risk_factors_groupcold_yn +  upper_resp + participantsex + fever + age

    #IgM
tab_uni_logmfi_igm<-tab_uni_logmfi%>%
filter(ig=="IgM")
tab_uni_logmfi_igg<-arrange(tab_uni_logmfi_igm, ag)
#write.csv (tab_uni_logmfi_igm, "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/Tables/raw_univar_logmfi_igm.csv")
#p-value <0.2 for IgM - RBD
#participantsex
#p-value <0.2 for IgM - S
#participantsex + bmi
#p-value <0.2 for IgM - RBD+S
#participantsex 

#MULTIVARIAT

#Predictor or "x" variables in the Multivariate for IgG-RBD:
l_varsx<-c("days_pso","participantsex", "age", "fever","upper_resp","risk_factors_groupcold_yn","allergy","loss_taste_smell", "bmi") #all isotype variables
mfi<-c("logmfi") #log10(mfi)
ig_ag<-c("ig", "ag")
gr_adj<-c('ig', 'vars_x', 'vars_y')
padj_method<-"BH" #verificar
arg_pval_arrange<-"raw_pval"


tab_multi_igg <- chmi.stat.lm_beta(
                          phen = regr, #change for db name
                          l_group = ig_ag, 
                          l_varsy = mfi, 
                          l_varsx = l_varsx, 
                          custom_formls= TRUE, #perquè volem X diferents per a cada isotip
                          formls = list(
                          'logmfi ~ days_pso +allergy + fever + loss_taste_smell', #IgA, make sure tu make it ordered (check levels)
                          'logmfi ~ days_pso +participantsex+age+fever+upper_resp+loss_taste_smell+risk_factors_groupcold_yn+allergy', #IgG
                          'logmfi ~ days_pso +participantsex + bmi'), #IgM
                          arg_multivar = 'multivariate', #can be 'multivariate'
                          arg_info_model = c('aic', 'bic', 'adj_r_squared'),
                          arg_info_beta = c('std_error'),
                          gr_adj = gr_adj, #do not change
                          format_pval = TRUE, #ajustat
                          arg_pval_arrange = arg_pval_arrange, #primer mirem el raw p-val
                          padj_method = padj_method) %>% #tipus de adjust de multiple comparisons
                          new_covidcat.beta_transform() %>% # es una funcio creada i carregada abans ( si fos log-log seria diferent)
                          mutate_at(vars(beta, ci_lower, ci_upper, raw_pval, aic, bic, adj_r_squared, std_error, adj_pval,
                              beta_transform_unit, beta_transform_percent, ci_lower_transform_unit, ci_lower_transform_percent, ci_upper_transform_unit, ci_upper_transform_percent), ~ round(.,4)) #round at 4 decimal
arrange(tab_multi_igg, raw_pval_text,ig, ag)
#write.csv (tab_multi_igg, "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/Tables/raw_multivar_logmfi.csv")

regr<-ddbb%>%
 filter(lab_sample_grouplab_sample_id!=952)%>% #we delete those ID's which are not complete
 filter(lab_sample_grouplab_sample_id!=999)%>%
 filter(studyn!=319)%>%
 filter(studyn!=1077)%>%
 filter(studyn!=325)%>%
 filter(studyn!=688)%>%
 filter(studyn!=981)%>%
 filter(Global_posneg=="pos")%>% # because we want to see what determines beig  positive for S once you're infected
 dplyr::select( "studyn","S_pos","participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy")%>% 
 unique() %>% droplevels()

#Get rid of: symptom_questions_groupsymptom_v,symptom_questions_groupsymptom_h because all answeres were "no"
length(unique(regr$studyn)) #56
summary(regr)

univar_outcome <- lapply(c("participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy"), function(var) {
                           
                           formula    <- as.formula(paste("S_pos~", var))
                           res.logist <- glm(formula, data = regr, family = binomial) #glm/lm depenent de l'outcome
                           #ci_univar_outcome <- confint(res.logist, level = .95)
                           results_univar_outcome <-summary(res.logist)$coefficients
                         })
univar_outcome
#Save the univariate in txt: change name, paste in a csv and separate by space, then re-name the columns:
#capture.output(univar_outcome, file = "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/univar_S_CI.txt")  

regr<-regr%>%filter(!is.na(lower_resp)) #have to delete 
length(unique(regr$studyn)) #55
model.outcome = glm(S_pos ~ participantsex+lower_resp+upper_resp+bmi, #glm/lm depenent de l'outcome
                    data=regr,
                    family = binomial)
summary(model.outcome)

# Now stepwise regression
step.model.outcome <- stepAIC(model.outcome)
coef(step.model.outcome)

model.outcome.final = glm(S_pos ~ participantsex+upper_resp+bmi,
                          data=regr,
                          family = binomial)
summary(model.outcome.final)

tab_model(model.outcome.final) #opens in an html with ORs

#I've created the RBD_pos variable at the begining of the code

#Now I get a dataframe from the big one with the variables I want to put in the model
regr<-ddbb%>%
 filter(lab_sample_grouplab_sample_id!=952)%>% #we delete those ID's which are not complete
 filter(lab_sample_grouplab_sample_id!=999)%>%
 filter(studyn!=319)%>%
 filter(studyn!=1077)%>%
 filter(studyn!=325)%>%
 filter(studyn!=688)%>%
 filter(studyn!=981)%>%
 filter(Global_posneg=="pos")%>% # because we want to see what determines beig  positive for S once you're infected
 dplyr::select( "studyn","RBD_pos","participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy")%>% 
 unique() %>% droplevels()
#Get rid of: symptom_questions_groupsymptom_v,symptom_questions_groupsymptom_h because all answeres were "no"

length(regr$RBD_pos) #56 --> correct, is the number of positive participants for any of the ig-ag pairs
summary(regr)

#Do once for the CIs and once for the results
univar_outcome <- lapply(c("participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy"), function(var) {                        
                           formula    <- as.formula(paste("RBD_pos~", var))
                           res.logist <- glm(formula, data = regr, family = binomial) #glm/lm depenent de l'outcome
                           results_univar_outcome <-summary(res.logist)$coefficients
                           #ci_univar_outcome <- confint(res.logist, level = .95)
                         })
univar_outcome

#Save the univariate in txt: change name, paste in a csv and separate by space, then re-name the columns:
#capture.output(univar_outcome, file = "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/univar_RBD.txt")  

# Variables with p-value <0.2 in the univariate tests enter the multivariate regression
model.outcome = glm(RBD_pos ~ age+participantsex+bmi+allergy, #glm/lm depenent de l'outcome
                    data=regr,
                    family = binomial)
summary(model.outcome)

# Now stepwise regression
step.model.outcome <- stepAIC(model.outcome)
coef(step.model.outcome)

model.outcome.final = glm(RBD_pos ~ participantsex+bmi+allergy, #glm/lm depenent de l'outcome
                    data=regr, family = binomial)
summary(model.outcome.final)
tab_model(model.outcome.final) #opens in an html with ORs

model<-as.data.frame(model.outcome.final)
model

model.outcome.final$Hochberg =
      p.adjust(model.outcome.final$ Raw.p,
               method = "hochberg")

#Now I get a dataframe from the big one with the variables I want to put in the model
regr<-ddbb%>%
 filter(lab_sample_grouplab_sample_id!=952)%>% #we delete those ID's which are not complete
 filter(lab_sample_grouplab_sample_id!=999)%>%
 filter(studyn!=319)%>%
 filter(studyn!=1077)%>%
 filter(studyn!=325)%>%
 filter(studyn!=688)%>%
 filter(studyn!=981)%>%
 filter(Global_posneg=="pos")%>% # because we want to see what determines being  positive for S once you're infected
 dplyr::select( "studyn","IgA_posneg","participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy")%>% 
 unique() %>% droplevels()
#Get rid of: symptom_questions_groupsymptom_v,symptom_questions_groupsymptom_h because all answeres were "no"

length(regr$IgA_posneg) #56 --> correct, is the number of positive participants for any of the ig-ag pairs
summary(regr)

univar_outcome <- lapply(c("participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy"), function(var) {
                           formula    <- as.formula(paste("IgA_posneg~", var))
                           res.logist <- glm(formula, data = regr, family = binomial) #glm/lm depenent de l'outcome
                           results_univar_outcome <-summary(res.logist)$coefficients
                           #ci_univar_outcome <- confint(res.logist, level = .95)
                         })
univar_outcome

#Save the univariate in txt: change name, paste in a csv and separate by space, then re-name the columns:
#capture.output(univar_outcome, file = "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/univar_IgA.txt")  

# Variables with p-value <0.2 in the univariate tests enter the multivariate regression
regr<-regr%>%filter(!is.na(risk_factors_groupsmoker_yn))
length(unique(regr$studyn))
model.outcome = glm(IgA_posneg ~ participantsex+age+fever+risk_factors_groupsmoker_yn+bmi, #glm/lm depenent de l'outcome
                    data=regr,
                    family = binomial)
summary(model.outcome)

# Now stepwise regression
step.model.outcome <- stepAIC(model.outcome)
coef(step.model.outcome)

model.outcome.final = glm(IgA_posneg ~ participantsex+age+risk_factors_groupsmoker_yn, #glm/lm depenent de l'outcome
                    data=regr, family = binomial)
summary(model.outcome.final)

#To transform the coefficient into % change the Beta coefficient in log-linear (categ or num):
#symptom_grouplost_tasteyes 
((10^0.15768)-1)*100

#To transform the coefficient into % change the Beta coefficient in log-log:
#((10^(-0.461256*log10(1.1)))-1)*100

tab_model(model.outcome.final) #opens in an html with ORs

#Now I get a dataframe from the big one with the variables I want to put in the model
regr<-ddbb%>%
 filter(lab_sample_grouplab_sample_id!=952)%>% #we delete those ID's which are not complete
 filter(lab_sample_grouplab_sample_id!=999)%>%
 filter(studyn!=319)%>%
 filter(studyn!=1077)%>%
 filter(studyn!=325)%>%
 filter(studyn!=688)%>%
 filter(studyn!=981)%>%
 filter(Global_posneg=="pos" )%>% # because we want to see what determines being  positive for S once you're infected 
 dplyr::select( "studyn","IgM_posneg","participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy")%>% 
 unique() %>% droplevels()
#Get rid of: symptom_questions_groupsymptom_v,symptom_questions_groupsymptom_h because all answeres were "no", autoinmune disease

length(regr$IgM_posneg) #36 --> correct, is the number of positive participants for any of the ig-ag pairs and with pso
summary(regr)

univar_outcome <- lapply(c("participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy"), function(var) {            
                           formula    <- as.formula(paste("IgM_posneg~", var))
                           res.logist <- glm(formula, data = regr, family = binomial) #glm/lm depenent de l'outcome
                           results_univar_outcome <-summary(res.logist)$coefficients
                           #ci_univar_outcome <- confint(res.logist, level = .95)
                         })
univar_outcome
#Save the univariate in txt: change name, paste in a csv and separate by space, then re-name the columns:
#capture.output(univar_outcome, file = "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/univar_IgM.txt")  

# Variables with p-value <0.2 in the univariate tests enter the multivariate regression:0

#Now I get a dataframe from the big one with the variables I want to put in the model
regr<-ddbb%>%
 filter(lab_sample_grouplab_sample_id!=952)%>% #we delete those ID's which are not complete
 filter(lab_sample_grouplab_sample_id!=999)%>%
 filter(studyn!=319)%>%
 filter(studyn!=1077)%>%
 filter(studyn!=325)%>%
 filter(studyn!=688)%>%
 filter(studyn!=981)%>%
 filter(Global_posneg=="pos" )%>% # because we want to see what determines being  positive for S once you're infected 
 dplyr::select( "studyn","IgG_posneg","participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy")%>% 
 unique() %>% droplevels()

regr$IgG_posneg<- as.character(regr$IgG_posneg) #we need 1 and 0, so I convert to character first
regr$IgG_posneg<-replace(regr$IgG_posneg, which(regr$IgG_posneg=="pos"),1)
regr$IgG_posneg[regr$IgG_posneg=="neg"]<-"0"
regr$IgG_posneg<- as.numeric(regr$IgG_posneg)
#Get rid of: symptom_questions_groupsymptom_v,symptom_questions_groupsymptom_h because all answeres were "no", autoinmune disease

length(regr$IgG_posneg) #36 --> correct, is the number of positive participants for any of the ig-ag pairs and with pso
summary(regr)

univar_outcome <- lapply(c("participantsex", "age", "fever","loss_taste_smell","lower_resp","upper_resp","dige","ache_fatigue",
               "risk_factors_groupcold_yn","risk_factors_groupsmoker_yn","bmi","allergy"), function(var) {            
                           formula    <- as.formula(paste("IgG_posneg~", var))
                           res.logist <- glm(formula, data = regr, family = binomial) #glm/lm depenent de l'outcome
                           results_univar_outcome <-summary(res.logist)$coefficients
                           # ci_univar_outcome <- confint(res.logist, level = .95)
                         })
univar_outcome
#Save the univariate in txt: change name, paste in a csv and separate by space, then re-name the columns:
#capture.output(univar_outcome, file = "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/univar_IgG.txt")  

# Variables with p-value <0.2 in the univariate tests enter the multivariate regression
regr<-regr%>%filter(!is.na(loss_taste_smell))
length(unique(regr$studyn))
model.outcome = glm(IgG_posneg ~ participantsex+age+fever+loss_taste_smell+dige+risk_factors_groupcold_yn+allergy, #glm/lm depenent de l'outcome
                    data=regr,
                    family = binomial)
summary(model.outcome)

# Now stepwise regression
step.model.outcome <- stepAIC(model.outcome)
coef(step.model.outcome)

model.outcome.final = glm(IgG_posneg ~ participantsex+age+risk_factors_groupcold_yn+allergy, #glm/lm depenent de l'outcome
                    data=regr, family = binomial)
summary(model.outcome.final)

tab_model(model.outcome.final) #opens in an html with ORs

#create a df with our vars and baseline visit only
occup<-ddbb%>%
dplyr::select('studyn','age','ageg','participantsex','socioecon_groupeducation','socioecon_groupoccupation','socioecon_groupemployed_before_y','socioecon_groupjob_loss_covid','socioecon_groupcovid_affect_fina','v137','socioecon_groupstruggle_expenses','v139','v140','v141')%>%
filter(age>=20)%>%
filter(socioecon_groupeducation!="other")%>%
mutate(affectation_yn=ifelse(is.na(socioecon_groupcovid_affect_fina),"No","Yes"))

#Filter out the retired: 
retired<-c(1, 42, 43, 54, 62, 74, 84, 85, 88, 134, 179, 205, 234, 301, 302, 311, 343, 352, 374, 390, 393, 398, 400, 401, 402, 405, 419, 468, 469, 487, 488, 541, 551, 572, 578, 582, 583, 599, 624, 626, 631, 658, 672, 681, 692, 702, 710, 716, 727, 766, 767, 803, 853, 861, 862, 874, 875, 882, 884, 889, 890, 905, 916,913, 917, 918, 929, 950, 1007, 1054, 1057, 1071, 1081, 1104, 1154, 1155,1557,1175, 1180, 1195, 1198, 1199, 1200, 1201, 1202, 1203,1205,1206) 
occup$socioecon_groupemployed_before_y<-as.character(occup$socioecon_groupemployed_before_y)
occup$socioecon_groupemployed_before_y[occup$studyn %in% c(retired)] <- "Retired"
occup$socioecon_groupemployed_before_y[occup$socioecon_groupoccupation==c("Ama de casa", "ama de casa")] <- "Ama de casa"
table(occup$socioecon_groupemployed_before_y)

#Those retired did not actually loose their job, but rather went into retirement, so we change the variable "do to COVID?" to NA, so as not to take it into account as lost job
occup$socioecon_groupjob_loss_covid<-replace(occup$socioecon_groupjob_loss_covid, which(occup$socioecon_groupemployed_before_y=="Retired"),NA)

occup$socioecon_groupeducation[occup$socioecon_groupeducation=="masters" | occup$socioecon_groupeducation=="PhD" ]<-"university_degree" #remove two categories
levels(occup$socioecon_groupeducation)<-list("Primary education"="school","Secondary education"='secondary_education',"Formación profesional"="technical_formacion_profesional", "University degree / Master / PhD"='university_degree')
table(occup$socioecon_groupeducation) # Other are all people underneath 20 and 1 that did not go to school

occup1<-occup%>% dplyr::select(-studyn, -socioecon_groupoccupation, -socioecon_groupcovid_affect_fina)

#We use the st function to directly do a table and tests:  @param group.test Set to \code{TRUE} to perform tests of whether each variable in the table varies over values of \code{group}.
#Only works with \code{group.long = FALSE}. Performs a joint F-test (using \code{anova(lm))}) for numeric variables, and a Chi-square test of independence (\code{chisq.test}) for categorical variables. 
#If you want to adjust things like which tests are used, significance star levels, etc., see the help file for \code{independence.test} and pass in a named list of options for that function.

occup2<-occup1
occup2[occup2=="unknown"]<-NA #we don't want unknown to be considered a category for the models
occup2<-droplevels(occup2) #to drop the unknown levels

#To see the p-value instead of the stars we have to create a named list with the opts from the format funcion in the independece.test:
my_list<-list(format='{name}:{stat}, Pr(>{name}):{pval}{stars}') #to get the pvalue and the stars coded


st(occup2, group = "socioecon_groupeducation", group.test=my_list, file="univariate_occup.csv", out="csv") #to save it in a csv and in the group.test we pass the list to get the numeric pvalues

#We use the st function to directly do a table and tests:  @param group.test Set to \code{TRUE} to perform tests of whether each variable in the table varies over values of \code{group}.
#Only works with \code{group.long = FALSE}. Performs a joint F-test (using \code{anova(lm))}) for numeric variables, and a Chi-square test of independence (\code{chisq.test}) for categorical variables. 
#If you want to adjust things like which tests are used, significance star levels, etc., see the help file for \code{independence.test} and pass in a named list of options for that function.

st(occup2, group = "participantsex", group.test=TRUE, file="univariate_occup.csv", out="csv") #to save it in a csv

exp<-ddbb%>%
dplyr::select('studyn','participantsex','age','socioecon_groupeducation','exposures_groupsanitary_centre','exposures_groupsanitary_centre_r','exposures_groupplaces_visited','exposures_groupplaces_visited_ou','exposures_groupsupermarket_pre','exposures_groupgrocery_pre','exposures_groupbank_pre','exposures_groupchurch_pre','exposures_grouphairdresser_pre','exposures_grouppamplona_pre','exposures_groupmadrid_pre','exposures_groupother_spain_pre','exposures_groupoutside_spain_pre','exposures_groupplaces_visited_su',
             'socioecon_groupwork_last_day','exposures_groupsanitary_work_dat','exposures_groupsanitary_other_da','exposures_groupsupermarket_date','exposures_groupgrocery_date','exposures_groupbank_date',
               'exposures_groupchurch_date','exposures_grouphairdresser_date','exposures_grouppamplona_date','exposures_groupmadrid_date','exposures_groupother_spain_date','exposures_groupoutside_spain_dat','exposures_groupcontact_covid_dat')%>%
filter(age>=18)

#Because the resposes are collected all together (because multiple answers were possible), I create variables with each of the responses separatedly and Y/No/
exp<-exp%>%
mutate(super=ifelse((grepl("supermarkets", exposures_groupplaces_visited) & exposures_groupsupermarket_date >= as.Date("2020-03-14")),"Yes", "No"))%>%
mutate(grocery=ifelse((grepl("small_grocery_store", exposures_groupplaces_visited)& exposures_groupgrocery_date >= as.Date("2020-03-14")),"Yes", "No"))%>%
mutate(bank=ifelse((grepl("bank", exposures_groupplaces_visited)& exposures_groupbank_date >= as.Date("2020-03-14")),"Yes", "No"))%>%
mutate(church=ifelse((grepl("church", exposures_groupplaces_visited)& exposures_groupchurch_date >= as.Date("2020-03-14")),"Yes", "No"))%>%
mutate(hairdresser=ifelse((grepl("hairdresser", exposures_groupplaces_visited)& exposures_grouphairdresser_date >= as.Date("2020-03-14")),"Yes", "No"))%>%
mutate(pamplona=ifelse((grepl("pamplona", exposures_groupplaces_visited)& exposures_grouppamplona_date >= as.Date("2020-03-14")),"Yes", "No"))%>%
mutate(madrid=ifelse((grepl("madrid", exposures_groupplaces_visited)& exposures_groupmadrid_date >= as.Date("2020-03-14")),"Yes", "No"))%>%
mutate(another_spain=ifelse((grepl("another_place_in_spain_specify", exposures_groupplaces_visited)& exposures_groupother_spain_date >= as.Date("2020-03-14")),"Yes", "No"))%>%
mutate(another_out_spain=ifelse((grepl("another_place_outside_spain_specify", exposures_groupplaces_visited)& exposures_groupoutside_spain_dat >= as.Date("2020-03-14")),"Yes", "No"))

#Reduce number of categories:
exp$exposures_groupsanitary_centre_r[exp$exposures_groupsanitary_centre_r=="accompanying_someone other"] <- "accompanying_someone"
exp$exposures_groupsanitary_centre_r[exp$exposures_groupsanitary_centre_r=="accompanying_someone symptoms_compatible_with_covid"] <- "accompanying_someone"
exp$exposures_groupsanitary_centre_r[exp$exposures_groupsanitary_centre_r=="accompanying_someone symptoms_incompatible_with_covid"] <- "accompanying_someone"
exp$exposures_groupsanitary_centre_r[exp$exposures_groupsanitary_centre_r=="symptoms_compatible_with_covid symptoms_incompatible_with_covid"] <- "symptoms_incompatible_with_covid"

#Keep the vars I want for the tables
cols<-c('exposures_groupsanitary_centre','exposures_groupsanitary_centre_r','super','grocery', 'bank', 'church', 'hairdresser', 'pamplona','madrid', 'another_spain', 'another_out_spain')
exp[cols]<-lapply(exp[cols],factor)

table<-exp%>%
dplyr::select('exposures_groupsanitary_centre','exposures_groupsanitary_centre_r','super','grocery', 'bank', 'church', 'hairdresser', 'pamplona','madrid', 'another_spain', 'another_out_spain')
#st(table, out="kable") #I paste the code in a new Marckdown cell and edit the state row
st(table, file="univariate_exposures.csv", out="csv") #to save it in a csv

a<-exp%>%
filter((super=="Yes" & exposures_groupsupermarket_date >= as.Date("2020-03-14")) | (grocery=="Yes"& exposures_groupgrocery_date >= as.Date("2020-03-14")))
length(unique(a$studyn)) #280
table(a$exposures_groupplaces_visited_su, useNA = "ifany")

table(exp$participantsex)

table_gender<-exp%>%
dplyr::select('participantsex','exposures_groupsanitary_centre','exposures_groupsanitary_centre_r','super','grocery', 'bank', 'church', 'hairdresser', 'pamplona','madrid', 'another_spain', 'another_out_spain')
#st(table_gender, group="participantsex", group.test=TRUE, out="kable") #I paste the code in a new Marckdown cell and edit the state row

#To see the p-value instead of the stars we have to create a named list with the opts from the format funcion in the independece.test:

my_list<-list(format='{name}:{stat}, Pr(>{name}):{pval}{stars}') #to get the pvalue and the stars coded
st(table_gender, group = "participantsex", group.test=my_list, file="univariate_exp.csv", out="csv") #to save it in a csv and in the group.test we pass the list to get the numeric pvalues

corr<-ddbb%>%
dplyr::select('Spike_IgM','RBD_IgM','Spike_IgG','RBD_IgG','Spike_IgA','RBD_IgA', "Global_posneg")
corr <- na.omit(corr)
#Log10 version because we see some outliers and the histogram is very skewed
vars<-c('Spike_IgM','RBD_IgM','Spike_IgG','RBD_IgG','Spike_IgA','RBD_IgA')
corr[vars]<-lapply(corr[vars],log10)
summary(corr)

corr.seropos<-corr %>% filter(Global_posneg=="pos") %>% dplyr::select(-Global_posneg)
plot <- rcorr(as.matrix(corr.seropos))
M <- plot$r
p_mat <- plot$P
corrplot(M, type = "upper", tl.col="black", tl.srt=45,
         p.mat = p_mat, sig.level = 0.05)# to order it:  order="hclust"

corr.seroneg<-corr %>% filter(Global_posneg=="neg") %>% dplyr::select(-Global_posneg)
plot <- rcorr(as.matrix(corr.seroneg))
M <- plot$r
p_mat <- plot$P
corrplot(M, type = "upper", tl.col="black", tl.srt=45,
         p.mat = p_mat, sig.level = 0.05)# to order it:  order="hclust"

#test normality after log
lapply(corr[vars], FUN=shapiro.test)

#Plot it with correlation coeff
corr.seropos<-corr %>% filter(Global_posneg=="pos") %>% dplyr::select(-Global_posneg)
plot <- rcorr(as.matrix(corr.seropos))
M <- plot$r
p_mat <- plot$P
corrplot.mixed(M, tl.col="black", upper = 'shade', p.mat = p_mat, sig.level = 0.05)# to order it:  order="hclust"


corr.seroneg<-corr %>% filter(Global_posneg=="neg") %>% dplyr::select(-Global_posneg)
plot <- rcorr(as.matrix(corr.seroneg))
M <- plot$r
p_mat <- plot$P
corrplot.mixed(M, tl.col="black", upper = 'shade', p.mat = p_mat, sig.level = 0.05)# to order it:  order="hclust"


#Another way to plot it, with Spearman correlation values and plots and histograms (we separate in colors pos and neg)
data<-corr

#We have to create the function for loess before (copied from https://www.r-bloggers.com/2016/02/multiple-regression-lines-in-ggpairs/):
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping, aes(color=Global_posneg, fill=Global_posneg, alpha=0.005)) + 
    geom_point() + 
    geom_smooth(method=loess) 
  p
}


plot<- GGally::ggpairs(corr, columns=1:6, aes(color=Global_posneg, fill=Global_posneg, alpha=0.005), lower = list(continuous = my_fn), upper = list(continuous = wrap("cor", method = 'spearman'))) +
#plot<- GGally::ggpairs(corr, columns=1:6, aes(color=Global_posneg, fill=Global_posneg, alpha=0.005), lower = list(continuous = wrap("smooth"),method=c('loess','spearman'),se= TRUE, alpha=0.01,color="black")) +
  scale_color_manual(values=c( "slateblue2", "mediumaquamarine"))+  
  scale_fill_manual(values=c( "slateblue2",  "mediumaquamarine"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #tilt the x titles
plot
#ggsave(plot, filename = "C:/Users/mribe/Documents/ESTUDIS/Master Global Health/Master Final Project/Seroprevalence of SARS-Cov-2 in Cizur/Paper Navarra/corr.png", width=8, height= 8)

symp_sex<-(ddbb%>% 
           filter(Global_posneg=="pos")%>%
           dplyr::select(participantsex,symptom_groupfever_37, symptom_groupfever_38, symptom_groupchills,symptom_grouplost_smell, symptom_grouplost_taste, symptom_groupcough,symptom_groupshort_breath, symptom_groupdyspnea, symptom_groupexpectorations, symptom_groupwheeze, symptom_groupchest_pain, symptom_groupsore_throat,
                                symptom_grouprunny_nose,symptom_groupnausea,symptom_groupab_pain,symptom_groupdairrhea,symptom_groupfatigue,symptom_groupmuscle_ache,symptom_groupheadache))
dim(symp_sex)
st(symp_sex, group = "participantsex", group.test=TRUE, file="univariate_symp_gender.csv", out="csv") #to save it in a csv

symp_sex<-(ddbb%>% 
           filter(Global_posneg=="pos")%>%
           dplyr::select(participantsex,symptom_groupfever_37, symptom_groupfever_38, symptom_groupchills,symptom_grouplost_smell, symptom_grouplost_taste, symptom_groupcough,symptom_groupshort_breath, symptom_groupdyspnea, symptom_groupexpectorations, symptom_groupwheeze, symptom_groupchest_pain, symptom_groupsore_throat,
                                symptom_grouprunny_nose,symptom_groupnausea,symptom_groupab_pain,symptom_groupdairrhea,symptom_groupfatigue,symptom_groupmuscle_ache,symptom_groupheadache))
dim(symp_sex)
st(symp_sex, group = "participantsex", group.test=TRUE, file="univariate_symp_gender_neg.csv", out="csv") #to save it in a csv

