# Load required packages
library(data.table)
library(readxl)
library(ggplot2)
library(devtools)

# Load custom SCM package
load_all('~/Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')


# Load dataset
dataset = read_xlsx('~/Dropbox/research/synth_control_paper/data/OIS_cali_replication_data/fillinpanel3.xlsx')
dataset <- data.table(dataset)

# Data preprocessing
dataset[is.na(ois), ois := 0]
dataset[is.na(unarmed), unarmed := 0]
dataset[is.na(mh), mh := 0]

# Create date variables and convert to year level
dataset[, mdate := as.Date(modate, format = "%Y-%m-%d")]
dataset[, year := year(mdate)]
dataset[, month := month(mdate)]

# Aggregate to year level
dataset <- dataset[, .(
  ois = sum(ois),
  unarmed = sum(unarmed),
  mh = sum(mh),
  officers = mean(officers),
  pop = mean(pop)
), by = c('stateid', 'year')]

# Calculate derived metrics
dataset[, oisp := (ois / pop) * 1000000]  # OIS per million population
dataset[, oiso := (ois / officers) * 1000] # OIS per thousand officers
dataset[, ofp := officers / pop]          # Officers per population

# Run specification curve analysis
sc <- spec_curve(
  dataset,
  outcomes = c('oiso', 'oisp'),
  col_name_unit_name = 'stateid',
  name_treated_unit = 'CA',

  # Different covariate specifications
  covagg = list(
    list(every_period = c('oisp', 'oiso')),
    list(every_period = c('oisp', 'oiso', 'ofp'))
  ),
  
  # Study parameters
  treated_period = 2020,
  min_period = 2015,
  end_period = 2022,
  col_name_period = 'year',
  num_pre_period_years = c(5),
  
  # Methodological choices
  feature_weights = c('uniform', 'optimize'),
  donor_sample = c('all', 'most_similar'),
  outcome_models = c('none', 'augsynth', 'lasso', 'ridge', 'ols'),
  constraints = list(list(name = 'simplex'), list(name = 'lasso'))
)

# Create visualization of specification curve
p <- plot_spec_curve(sc, name_treated_unit = 'CA')
print(p)

# Extract results into a dataframe for further analysis
sc_results_list <- list()

# Loop through all model specifications and compile results
for(outcome in names(sc)) {
  for(const in names(sc[[outcome]])) {
    for(fw in names(sc[[outcome]][[const]])) {
      for(feat in names(sc[[outcome]][[const]][[fw]])) {
        for(ds in names(sc[[outcome]][[const]][[fw]][[feat]])) {
          for(ny in names(sc[[outcome]][[const]][[fw]][[feat]][[ds]])) {
            # Get the current model estimate
            estee <- sc[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]$estimate
            
            # Compile inference results with metadata
            combined <- rbindlist(sc[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]$infer)
            combined$outcome <- outcome
            combined$const <- const
            combined$fw <- fw
            combined$data_sample <- ds
            combined$feat <- feat
            combined$num_pre_period_years <- ny
            
            # Calculate RMSE for pre-intervention fit
            combined$rmse <- sqrt(mean((estee$data$Y.pre - estee$est.results$Y.pre.fit)^2))
            
            # Add to results list
            sc_results_list[[length(sc_results_list) + 1]] <- combined
          }
        }
      }
    }
  }
}

# Combine all results into a single dataframe
sc_results_df <- rbindlist(sc_results_list)

# Analyze distribution of model fit
plot(density(sc_results_df$rmse))





















library(readstata13)
library(data.table)
library(readxl)
library(devtools)
library(ggplot2)

load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')


#### kaplan dataset
#### hogan dataset
#### german reunif
#### california












































dataset[, oisp := (ois / pop) * 1000000]
dataset[, oiso := (ois / officers) * 1000]
datset[, mdate := as.Date(modate, format = "%Y-%m-%d")] # Adjust format if needed
dataset[, mdate := as.Date(modate, format = "%Y-%m-%d")] # Adjust format if needed
dataset[, year := year(mdate)]
dataset[, month := month(mdate)]
dataset

# Set up panel data structure
dataset[, sid := .GRP, by = stateid]
setkey(dataset, sid, mdate)



# SCM Analysis
# Function to generate time periods for predictors
gen_predictors <- function(start_year, start_month, end_year, end_month) {
  dates <- seq(as.Date(paste(start_year, start_month, "01", sep = "-")),
               as.Date(paste(end_year, end_month, "01", sep = "-")), by = "month")
  return(paste0("ois(`=tm(", format(dates, "%Ym%d"), ")')"))
}
dataset = read_xlsx('~/Dropbox/research/synth_control_paper/data/OIS_cali_replication_data/fillinpanel3.xlsx')
dataset = data.table(dataset)



predictors <- gen_predictors(2015, 1, 2019, 12)  # Adjust years if needed

dataset[is.na(mh), mh:0]
dataset[is.na(ois), ois:=0]
dataset[is.na(unarmed), unarmed:=0]
dataset[, oisp := (ois / pop) * 1000000]
dataset[, oiso := (ois / officers) * 1000]
dataset[, mdate := as.Date(modate, format = "%Y-%m-%d")] # Adjust format if needed
dataset[, year := year(mdate)]
dataset[, month := month(mdate)]
# Set up panel data structure
dataset[, sid := .GRP, by = stateid]
setkey(dataset, sid, mdate)
dataset
as.Date('2020-01-01')
dataset[, .(ois=sum(ois), unarmed=sum(unarmed), mh=sum(mh),
            officers=mean(officers), pop=mean(pop)),
        by=c('stateid', 'year')]
dataset = dataset[, .(ois=sum(ois), unarmed=sum(unarmed), mh=sum(mh),
                      officers=mean(officers), pop=mean(pop)),
                  by=c('stateid', 'year')]
colMeans(is.na(dataset))
dataset[, ofp:=officers/pop]
sc = spec_curve(
  dataset,
  outcomes='oisp',
  col_name_unit_name='stateid',
  name_treated_unit='CA',
  covagg=list(
    list(every_period=c('oisp')),
    list(every_period=c('oisp', 'officers')),
    list(every_period=c('oisp', 'pop')),
    list(every_period=c('oisp', 'pop', 'officers'))),
  treated_period=2020,#as.Date('2020-01-01'),
  min_period=2015,
  end_period=2022,
  col_name_period='year',
  num_pre_period_years = c(5),
  feature_weights=c('uniform', 'optimize'),
  donor_sample = c('all', 'most_similar'),
  outcome_models=c('none', 'augsynth', 'lasso', 'ridge', 'ols'),
  constraints=list(
    list(name='simplex')#,
    # list(name='lasso')
  )
)
dataset[, oisp := (ois / pop) * 1000000]
dataset[, oiso := (ois / officers) * 1000]
dataset[, mdate := as.Date(modate, format = "%Y-%m-%d")] # Adjust format if needed
dataset[, year := year(mdate)]
dataset[, month := month(mdate)]
# Set up panel data structure
dataset[, sid := .GRP, by = stateid]
setkey(dataset, sid, mdate)
dataset = dataset[, .(ois=sum(ois), unarmed=sum(unarmed), mh=sum(mh),
                      officers=mean(officers), pop=mean(pop)),
                  by=c('stateid', 'year')]
dataset[, ofp:=officers/pop]
dataset[, ofp:=officers/pop]
dataset[, oisp := (ois / pop) * 1000000]
dataset[, oiso := (ois / officers) * 1000]
sc = spec_curve(
  dataset,
  outcomes='oisp',
  col_name_unit_name='stateid',
  name_treated_unit='CA',
  covagg=list(
    list(every_period=c('oisp')),
    list(every_period=c('oisp', 'officers')),
    list(every_period=c('oisp', 'pop')),
    list(every_period=c('oisp', 'pop', 'officers'))),
  treated_period=2020,#as.Date('2020-01-01'),
  min_period=2015,
  end_period=2022,
  col_name_period='year',
  num_pre_period_years = c(5),
  feature_weights=c('uniform', 'optimize'),
  donor_sample = c('all', 'most_similar'),
  outcome_models=c('none', 'augsynth', 'lasso', 'ridge', 'ols'),
  constraints=list(
    list(name='simplex')#,
    # list(name='lasso')
  )
)
warnings()
sc
p=plot_spec_curve(sc, name_treated_unit = 'CA')
p
sc = spec_curve(
  dataset,
  outcomes=c('oiso', 'oisp'),
  col_name_unit_name='stateid',
  name_treated_unit='CA',
  covagg=list(
    list(every_period=c('oisp', 'oiso')),
    list(every_period=c('oisp', 'oiso', 'ofp'))),
  treated_period=2020,#as.Date('2020-01-01'),
  min_period=2015,
  end_period=2022,
  col_name_period='year',
  num_pre_period_years = c(5),
  feature_weights=c('uniform', 'optimize'),
  donor_sample = c('all', 'most_similar'),
  outcome_models=c('none', 'augsynth', 'lasso', 'ridge', 'ols'),
  constraints=list(
    list(name='simplex')#,
    # list(name='lasso')
  )
)



p=plot_spec_curve(sc, name_treated_unit = 'CA')
p
sc$oiso$simplex$uniform[['every_period__c("oisp", "oiso")']]
sc$oiso$simplex$uniform[['every_period__c("oisp", "oiso")']]$most_similar$n_pp_years_5$infer
sc$oiso$simplex$uniform[['every_period__c("oisp", "oiso")']]$most_similar$n_pp_years_5$sc.est
sc$oiso$simplex$uniform[['every_period__c("oisp", "oiso")']]$most_similar$n_pp_years_5$sc_est
sc$oiso$simplex$uniform[['every_period__c("oisp", "oiso")']]$most_similar$n_pp_years_5$scm
sc$oiso$simplex$uniform[['every_period__c("oisp", "oiso")']]$most_similar$n_pp_years_5$estimate
a=sc$oiso$simplex$uniform[['every_period__c("oisp", "oiso")']]$most_similar$n_pp_years_5$estimate
a
a$est.results$A.hat
a$est.results$Y.pre.fit
a$est.results$Y.pre.fit - a$est.results$yhats$lasso$control
a$est.results$yhats$lasso$control
a$est.results$yhats$lasso$treated
a=sc$oiso$simplex$uniform$`every_period__c("oisp", "oiso")`
a$all$n_pp_years_5$estimate$est.results$yhats$lasso$control
a$all$n_pp_years_5$estimate$est.results$yhats$ridge$control
a$all$n_pp_years_5$estimate$est.results$yhats$ols$control
a
a$all$n_pp_years_5$estimate$est.results$Y.pre.fit
a$all$n_pp_years_5$estimate$data$Y.pre
sc_results_list = list()
for(outcome in names(sc)){
  for(const in names(sc[[outcome]])){
    for(fw in names(sc[[outcome]][[const]])){
      for(feat in names(sc[[outcome]][[const]][[fw]])){
        for(ds in names(sc[[outcome]][[const]][[fw]][[feat]])){
          for(ny in names(sc[[outcome]][[const]][[fw]][[feat]][[ds]])){
            combined = rbindlist(sc[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]$infer)
            combined$outcome = outcome
            combined$const = const
            combined$fw = fw
            combined$data_sample = ds
            combined$feat = feat
            combined$num_pre_period_years = ny
            sc_results_list[[length(sc_results_list) + 1]] = combined
          }
        }
      }
    }
  }
}
sc_results_df = rbindlist(sc_results_list)
sc_results_df






estee = sc[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]$estimate
sqrt(mean((estee$data$Y.pre-estee$est.results$Y.pre.fit)^2))
sc_results_list = list()
for(outcome in names(sc)){
  for(const in names(sc[[outcome]])){
    for(fw in names(sc[[outcome]][[const]])){
      for(feat in names(sc[[outcome]][[const]][[fw]])){
        for(ds in names(sc[[outcome]][[const]][[fw]][[feat]])){
          for(ny in names(sc[[outcome]][[const]][[fw]][[feat]][[ds]])){
            estee = sc[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]$estimate
            combined = rbindlist(sc[[outcome]][[const]][[fw]][[feat]][[ds]][[ny]]$infer)
            combined$outcome = outcome
            combined$const = const
            combined$fw = fw
            combined$data_sample = ds
            combined$feat = feat
            combined$num_pre_period_years = ny
            combined$rmse = sqrt(mean((estee$data$Y.pre-estee$est.results$Y.pre.fit)^2))
            sc_results_list[[length(sc_results_list) + 1]] = combined
          }
        }
      }
    }
  }
}
sc_results_df = rbindlist(sc_results_list)
plot(density(combined$rmse))
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
p=plot_spec_curve(sc, name_treated_unit = 'CA',.09)
p
p=plot_spec_curve(sc, name_treated_unit = 'CA',.04)
p



# load dataset
dataset = fread('Dropbox/research/synth_control_paper/repo/variation_scm/data/kaplan_dataset_processed.csv')
dataset[, hr_rate:=num_homicide/population]
# proprocess
dataset[, ori9:=gsub(',| |-|\\.', '_', ori9)]
dataset[, ori9:=gsub('_+', '_', ori9)]
dataset[, num_homicide:=as.numeric(num_homicide)]
dataset[, hr_rate:=as.numeric(hr_rate)]
gc()
# load dataset
dataset = fread('Dropbox/research/synth_control_paper/repo/variation_scm/data/kaplan_dataset_processed.csv')
dataset[, hr_rate:=num_homicide/population]
# proprocess
dataset[, ori9:=gsub(',| |-|\\.', '_', ori9)]
dataset[, ori9:=gsub('_+', '_', ori9)]
dataset[, num_homicide:=as.numeric(num_homicide)]
dataset[, hr_rate:=as.numeric(hr_rate)]
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
sc = spec_curve(
  dataset,
  outcomes=c('num_homicide'),
  col_name_unit_name='ori9',
  name_treated_unit='PAPEP0000',
  covagg=list(
    list(every_period=c('num_homicide', 'hr_rate')),
    list(every_period=c('num_homicide', 'cleared_cases', 'hr_rate'))),
  treated_period=2015,
  min_period=2010,
  end_period=2019,
  col_name_period='year',
  feature_weights=c('uniform', 'optimize'),
  donor_sample = c('all', 'most_similar'),
  outcome_models=c('none', 'augsynth', 'lasso', 'ridge', 'ols'),
  constraints=list(
    list(name='simplex'),
    list(name='lasso')
  )
)


warnings()
a=plot_spec_curve(sc, name_treated_unit='PAPEP0000', rmse_threshold = Inf)
a
bb=sc$num_homicide$simplex$uniform$`every_period__c("num_homicide", "hr_rate")`
bb$all$n_pp_years_5$estimate$data$Y.donors
bb$all$n_pp_years_5$estimate$treated_period
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a=plot_spec_curve(sc, name_treated_unit='PAPEP0000', rmse_threshold = Inf, normalize_outcomes=FALSE)
a
a=plot_spec_curve(sc, name_treated_unit='PAPEP0000', rmse_threshold = Inf, normalize_outcomes=TRUE)
bb$all$n_pp_years_5$estimate$data
bb$all$n_pp_years_5$estimate$data$Y.donors
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a=plot_spec_curve(sc, name_treated_unit='PAPEP0000', rmse_threshold = Inf, normalize_outcomes=TRUE)
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a=plot_spec_curve(sc, name_treated_unit='PAPEP0000', rmse_threshold = Inf, normalize_outcomes=FALSE)
a=plot_spec_curve(sc, name_treated_unit='PAPEP0000', rmse_threshold = Inf, normalize_outcomes=TRUE)
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a=plot_spec_curve(sc, name_treated_unit='PAPEP0000', rmse_threshold = Inf, normalize_outcomes=TRUE)
bb$all$n_pp_years_5$estimate$data$Y.donors
sd(bb$all$n_pp_years_5$estimate$data$Y.donors)
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a=plot_spec_curve(sc, name_treated_unit='PAPEP0000', rmse_threshold = Inf, normalize_outcomes=TRUE)
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a=plot_spec_curve(sc, name_treated_unit='PAPEP0000', rmse_threshold = Inf, normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = .05,
                    normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = .02,
                    normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = .01,
                    normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = .0001,
                    normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 1e-12,
                    normalize_outcomes=TRUE)
a
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a = plot_spec_curve(sc,
                    outcomes = c('num_homicide'),
                    name_treated_unit='PAPEP0000',
                    # rmse_threshold = 1,
                    normalize_outcomes=TRUE)
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a = plot_spec_curve(sc,
                    outcomes = c('num_homicide'),
                    name_treated_unit='PAPEP0000',
                    # rmse_threshold = 1,
                    normalize_outcomes=TRUE)
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a = plot_spec_curve(sc,
                    outcomes = c('num_homicide'),
                    name_treated_unit='PAPEP0000',
                    # rmse_threshold = Inf,
                    normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    outcomes = c('num_homicide'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 10,
                    normalize_outcomes=TRUE)
a
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a = plot_spec_curve(sc,
                    outcomes = c('num_homicide'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 10,
                    normalize_outcomes=TRUE)
a = plot_spec_curve(sc,
                    outcomes = c('num_homicide'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 5,
                    normalize_outcomes=TRUE)
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a = plot_spec_curve(sc,
                    outcomes = c('num_homicide'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 5,
                    normalize_outcomes=TRUE)
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
a = plot_spec_curve(sc,
                    outcomes = c('num_homicide'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 5,
                    normalize_outcomes=TRUE)
a = plot_spec_curve(sc,
                    outcomes = c('num_homicide'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 20,
                    normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    outcomes = c('num_homicide'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 100,
                    normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    outcomes = c('hr_rate'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = .01,
                    normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    outcomes = c('hr_rate'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 20,
                    normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    outcomes = c('hr_rate'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 1e-10,
                    normalize_outcomes=TRUE)
a
a = plot_spec_curve(sc,
                    outcomes = c('hr_rate'),
                    name_treated_unit='PAPEP0000',
                    rmse_threshold = 1e-12,
                    normalize_outcomes=TRUE)
a
load_all('Dropbox/research/synth_control_paper/repo/SCMs/SCMs/')
sc = spec_curve(
  dataset,
  outcomes=c('num_homicide', 'hr_rate'),
  col_name_unit_name='ori9',
  name_treated_unit='PAPEP0000',
  covagg=list(
    list(every_period=c('num_homicide', 'hr_rate')),
    list(every_period=c('num_homicide', 'cleared_cases', 'hr_rate'))),
  treated_period=2015,
  min_period=2010,
  end_period=2019,
  col_name_period='year',
  feature_weights=c('uniform', 'optimize'),
  donor_sample = c('all', 'most_similar'),
  outcome_models=c('none', 'augsynth', 'lasso', 'ridge', 'ols'),
  constraints=list(
    list(name='simplex'),
    list(name='lasso')
  )
)





