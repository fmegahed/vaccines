
# This code is written to fit an example of a logistic growth curve on vaccine data
# We capitalize on the example from 
# https://bscheng.com/2014/05/07/modeling-logistic-growth-data-in-r/

if(require(pacman)==FALSE) install.packages("pacman")

pacman::p_load(tidyverse, magrittr, janitor, lubridate, # data analysis packages
               car, minpack.lm,
               ggpubr)



# * Getting the Vaccine Data ----------------------------------------------

csvLink = "https://data.cdc.gov/api/views/8xkx-amqh/rows.csv?accessType=DOWNLOAD"
rawVaccinesData = read_csv(csvLink)

vacCSVtime = Sys.time()

vaccines = rawVaccinesData %>% 
  # Removing unknown (UNK) and territories/states outside of continental US
  filter(!Recip_State %in% c('AK', 'AS', 'FM', 'GU', 'HI', 
                             'MH', 'MP', 'PR', 'PW', 'UNK', 'VI')) %>% 
  # Converting Date from Char to Date format
  mutate(Date = mdy(Date)) %>% 
  # Ascending order by FIPS and Date
  arrange(FIPS, Date) %>% 
  # Reducing the columns to those potentially relevant
  select(Date, FIPS, MMWR_week, Recip_County, Recip_State,
         Series_Complete_Pop_Pct, Completeness_pct) %>% 
  # Converting Char to Factor cols  
  mutate_if(is.character, as_factor)


# * Extracting Butler County for our Example ------------------------------

butler_county = vaccines %>% 
  filter(Recip_State == 'OH' & Recip_County == 'Butler County')

write_rds(butler_county, file = 'butler_county_vaccines.rds')


# * Creating the Data for Model Fitting -----------------------------------

butler_county = butler_county %>% 
  janitor::clean_names() %>% 
  mutate(day = day(date),
         month = month(date))

butler_modeling_tbl = 
  butler_county %>% 
  dplyr::filter(day == 1) %>% 
  mutate(months_from_start = row_number()) 


# * Butler County -----------------------------------------------------

nls_model_butler = nls(
  series_complete_pop_pct~phi1/( 1+exp(- (phi2+phi3*months_from_start) ) ),
  start=list(phi1=50, phi2=-3, phi3=0.8),
  data=butler_modeling_tbl,trace=TRUE)

phi1 = coef(nls_model_butler)[1]
phi2 = coef(nls_model_butler)[2]
phi3 = coef(nls_model_butler)[3]

prediction_df = tibble(
  months_from_start= butler_modeling_tbl$months_from_start,
  predicted_pop_vaccine = phi1/(1+exp(-(phi2+phi3*months_from_start))) 
)

combined_tbl = 
  butler_modeling_tbl %>% 
  left_join(y = prediction_df, by = 'months_from_start')

colors_brewer = RColorBrewer::brewer.pal(n = 10, name = 'Dark2')

cols <- c("Actual"= colors_brewer[1], 'Logistic Growth Curve Fitted' = colors_brewer[2])

combined_tbl %>% 
  ggplot( aes(x = months_from_start ))+
  geom_point(aes(y = series_complete_pop_pct, color = 'Actual', shape = 'Actual'), size = 4)+
  geom_point(aes(y = predicted_pop_vaccine, color = 'Logistic Growth Curve Fitted',
                 shape = 'Logistic Growth Curve Fitted'),  stroke = 1.25, size = 4)+
  geom_line(aes(y = predicted_pop_vaccine, color = 'Logistic Growth Curve Fitted'), 
            size=1, linetype = 5) +
  scale_color_manual(name = 'Data', values = cols) +
  scale_shape_manual(name = 'Data', values = c("Actual"= 19, 'Logistic Growth Curve Fitted' = 1)) +
  labs(x = 'Months from Start', y = '% Vaccinated',  title = 'Butler County (Ohio) with complete (i.e., no missing) data') +
  scale_x_continuous(breaks = scales::pretty_breaks(n= 20)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n= 8), limits = c(0,65)) +
  theme_bw() +
  theme(legend.position =c(.85,.2)) -> p_complete

p_complete


# * Harris County ---------------------------------------------------------

harris_county = vaccines %>% 
  filter(Recip_State == 'TX' & Recip_County == 'Harris County')

harris_county = harris_county %>% 
  janitor::clean_names() %>% 
  mutate(day = day(date),
         month = month(date))

harris_modeling_tbl = 
  harris_county %>% 
  dplyr::filter(day == 1) %>% 
  mutate(months_from_start = row_number(),
         series_complete_pop_pct = ifelse(series_complete_pop_pct == 0,
                                          NA, series_complete_pop_pct)) 

nls_model_harris = nls(
  series_complete_pop_pct~phi1/( 1+exp(- (phi2+phi3*months_from_start) ) ),
  start=list(phi1=65, phi2=0, phi3=0.1),
  data=harris_modeling_tbl,trace=TRUE)

phi1 = coef(nls_model_harris)[1]
phi2 = coef(nls_model_harris)[2]
phi3 = coef(nls_model_harris)[3]

prediction_df = tibble(
  months_from_start= harris_modeling_tbl$months_from_start,
  predicted_pop_vaccine = phi1/(1+exp(-(phi2+phi3*months_from_start))) 
)

cols <- c("Actual"= colors_brewer[1], 'Bayesian Imputed Curve' = colors_brewer[3])

combined_tbl = 
  harris_modeling_tbl %>% 
  left_join(y = prediction_df, by = 'months_from_start')

combined_tbl %>% 
  ggplot( aes(x = months_from_start ))+
  geom_point(aes(y = series_complete_pop_pct, color = 'Actual', shape = 'Actual'), size = 4)+
  geom_point(aes(y = predicted_pop_vaccine, color = 'Bayesian Imputed Curve',
                 shape = 'Bayesian Imputed Curve'),  stroke = 1.25, size = 4)+
  geom_line(aes(y = predicted_pop_vaccine, color = 'Bayesian Imputed Curve'), 
            size=1, linetype = 5) +
  scale_color_manual(name = 'Data', values = cols) +
  scale_shape_manual(name = 'Data', values = c("Actual"= 19, 'Bayesian Imputed Curve' = 1)) +
  labs(x = 'Months from Start', y = '% Vaccinated',  title = 'Harris County (Texas) with incomplete data') +
  scale_x_continuous(breaks = scales::pretty_breaks(n= 20)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n= 8), limits = c(0,65)) +
  theme_bw() +
  theme(legend.position =c(.85,.2)) -> p_incomplete

p_incomplete



# * County in California with Completely Missing Data ---------------------

mariposa_county = vaccines %>% 
  filter(Recip_State == 'CA' & Recip_County == 'Mariposa County')

mariposa_county = mariposa_county %>% 
  janitor::clean_names() %>% 
  mutate(day = day(date),
         month = month(date))

mariposa_modeling_tbl = 
  mariposa_county %>% 
  dplyr::filter(day == 1) %>% 
  mutate(months_from_start = row_number(),
         series_complete_pop_pct = ifelse(series_complete_pop_pct == 0,
                                          NA, series_complete_pop_pct)) 

# we cannot fit a model to mariposa since it has no data
# so we will be using an example from another county to talk about the expected
# curve
prediction_df = tibble(
  months_from_start= harris_modeling_tbl$months_from_start,
  predicted_pop_vaccine = phi1/(1+exp(-(phi2+phi3*months_from_start))) 
)

combined_tbl = 
  mariposa_modeling_tbl %>% 
  left_join(y = prediction_df, by = 'months_from_start')

cols <- c( 'Bayesian Imputed Curve' = colors_brewer[3])

combined_tbl %>% 
  ggplot( aes(x = months_from_start ))+
  # geom_point(aes(y = series_complete_pop_pct, color = 'Actual', shape = 'Actual'), size = 4)+
  geom_point(aes(y = predicted_pop_vaccine, color = 'Bayesian Imputed Curve',
                 shape = 'Bayesian Imputed Curve'),  stroke = 1.25, size = 4)+
  geom_line(aes(y = predicted_pop_vaccine, color = 'Bayesian Imputed Curve'), 
            size=1, linetype = 5) +
  scale_color_manual(name = 'Data', values = cols) +
  scale_shape_manual(name = 'Data', values = c('Bayesian Imputed Curve' = 1)) +
  labs(x = 'Months from Start', y = '% Vaccinated',  title = 'Mariposa County (California) with censored data') +
  scale_x_continuous(breaks = scales::pretty_breaks(n= 20)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n= 8), limits = c(0,65)) +
  theme_bw() +
  theme(legend.position =c(.85,.2)) -> p_total_missing

p_total_missing


ggarrange(p_complete, p_incomplete, p_total_missing, ncol = 1, common.legend = F)

ggsave(filename = 'hypothetical_scenarios.png')
