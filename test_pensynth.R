library(devtools)
load_all('.')

# Load basque data
basque_data <- read.csv('inst/extdata/basque.csv')

# Check data availability
basque_subset <- basque_data[basque_data$regionname == 'Basque Country (Pais Vasco)', ]
print("Basque Country data:")
print(basque_subset[1:10, c('year', 'gdpcap')])

# Check years with non-missing gdpcap
non_missing_years <- basque_subset[!is.na(basque_subset$gdpcap), 'year']
print(paste("Years with data:", min(non_missing_years), "to", max(non_missing_years)))

# Create scdata object with available years
basque_scdata <- scdata(
  df = basque_data,
  id.var = 'regionname',
  outcome.var = 'gdpcap', 
  time.var = 'year',
  period.pre = 1955:1969,
  period.post = 1970:1975,  # Shorter period to avoid missing data
  unit.tr = 'Basque Country (Pais Vasco)',
  unit.co = c('Cataluna', 'Madrid (Comunidad De)', 'Andalucia')  # Just a few units
)

print('Testing pensynth with fixed lambda...')
result <- scest(data = basque_scdata, w.constr = list(name = 'pensynth', lambda = 0.1))
print('Success!')