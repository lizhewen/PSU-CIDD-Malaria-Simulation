import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from constant import HEADLINE, ENCODINGDB, REPORTDAYS

def xaxis_label_ticker(scale_x=365, burnin_year=10):
  # Ticker function that labels x-axis in years
  #   and reflect burn-in phase of the simulation
  # Input:
  #   `scale_x` - how many days in a year, default=365
  #   `burnin_year` - how long is the burn-in phase in years, default=10
  # Return:
  #   A ticker function that labels the first day after burn-in as 0 on x-axis

  return (ticker.FuncFormatter(lambda x, pos: '{0:g}'.format((x-burnin_year*365)/scale_x)))

def IQR_compute(rawfilepattern, parsedatapattern, settherapyinfo):
  # This function takes in 100 output files from 
  #   simulation with the same setting, parsing them 
  #   by adding appropriate headers, and then compute
  #   the IQR for the genotype frequency and returns 
  #   three dataframes (LQ, Median, and UQ).
  allgenodf = []
  for i in range(131):
    # 361 rows of data are in output file
    allgenodf.append(pd.DataFrame(columns = range(361)))
  from constant import HEADLINE, ENCODINGDB, REPORTDAYS
  for j in range(1,101):
    file_name = rawfilepattern % j
    dummy_file = parsedatapattern % j
    # open original file in read mode and dummy file in write mode
    with open(file_name, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
      # Write given line to the dummy file
      write_obj.write(HEADLINE + '\n')
      # Read lines from original file one by one and append them to the dummy file
      for line in read_obj:
        write_obj.write(line)
    df = pd.read_csv(dummy_file, sep='\t')
    # For each file - append genotype freq to 128 sep dfs
    for (genotype, dfcounter) in zip(ENCODINGDB, range(128)):
      allgenodf[dfcounter] = allgenodf[dfcounter].append(df[genotype], ignore_index=True)
    # For each file - append bsp, population, and ntf to the last df
    allgenodf[128] = allgenodf[128].append(df['population'], ignore_index=True)
    allgenodf[129] = allgenodf[129].append(df['blood_slide_prevalence'], ignore_index=True)
    allgenodf[130] = allgenodf[130].append(df['monthly_ntf_raw'], ignore_index=True)
  
  L_IQRdf = pd.DataFrame(columns = ['time_elapsed', 'population', 'blood_slide_prevalence', 'monthly_ntf_raw', 'ntf_percent'] + ENCODINGDB)
  M_IQRdf = pd.DataFrame(columns = ['time_elapsed', 'population', 'blood_slide_prevalence', 'monthly_ntf_raw', 'ntf_percent'] + ENCODINGDB)
  U_IQRdf = pd.DataFrame(columns = ['time_elapsed', 'population', 'blood_slide_prevalence', 'monthly_ntf_raw', 'ntf_percent'] + ENCODINGDB)

  # Recombine population IQR results to three dfs
  population_IQR_result = allgenodf[128].astype('float').quantile([.25, .5, .75]).values.tolist()
  L_IQRdf['population'] = population_IQR_result[0]
  M_IQRdf['population'] = population_IQR_result[1]
  U_IQRdf['population'] = population_IQR_result[2]

  # Recombine bsp IQR results to three dfs
  bsp_IQR_result = allgenodf[129].astype('float').quantile([.25, .5, .75]).values.tolist()
  L_IQRdf['blood_slide_prevalence'] = bsp_IQR_result[0]
  M_IQRdf['blood_slide_prevalence'] = bsp_IQR_result[1]
  U_IQRdf['blood_slide_prevalence'] = bsp_IQR_result[2]

  # Recombine NTF IQR results to three dfs
  ntf_IQR_result = allgenodf[130].astype('float').quantile([.25, .5, .75]).values.tolist()
  L_IQRdf['monthly_ntf_raw'] = ntf_IQR_result[0]
  M_IQRdf['monthly_ntf_raw'] = ntf_IQR_result[1]
  U_IQRdf['monthly_ntf_raw'] = ntf_IQR_result[2]

  # Recombine genotype freq IQR results to three dfs
  for k in range(128):
    temp_IQR_result = allgenodf[k].astype('float').quantile([.25, .5, .75]).values.tolist()
    L_IQRdf[ENCODINGDB[k]] = temp_IQR_result[0]
    M_IQRdf[ENCODINGDB[k]] = temp_IQR_result[1]
    U_IQRdf[ENCODINGDB[k]] = temp_IQR_result[2]

  L_IQRdf['time_elapsed'] = REPORTDAYS
  M_IQRdf['time_elapsed'] = REPORTDAYS
  U_IQRdf['time_elapsed'] = REPORTDAYS

  L_IQRdf['ntf_percent'] = L_IQRdf['monthly_ntf_raw'] / L_IQRdf['population'] * 100
  M_IQRdf['ntf_percent'] = M_IQRdf['monthly_ntf_raw'] / M_IQRdf['population'] * 100
  U_IQRdf['ntf_percent'] = U_IQRdf['monthly_ntf_raw'] / U_IQRdf['population'] * 100

  L_IQRdf.to_csv('IQR_data/%s_L.csv' % settherapyinfo, index=False)
  M_IQRdf.to_csv('IQR_data/%s_M.csv' % settherapyinfo, index=False)
  U_IQRdf.to_csv('IQR_data/%s_U.csv' % settherapyinfo, index=False)
