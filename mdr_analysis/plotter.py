import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from plot_helper import xaxis_label_ticker, coloring_legend, df_col_replace

FIRST_ROW_AFTER_BURNIN = 120
ANNOTATION_X_LOCATION = 3833

def fig1_plot_IQR(ax, df_l, df_m, df_u, drug, annoty=None):
  df_m = df_m.iloc[FIRST_ROW_AFTER_BURNIN:]
  df = df_col_replace(df_m, drug, option=1)
  df_nl = df_col_replace(df_l, drug, option=1)
  df_nl = df_nl.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_nu = df_col_replace(df_u, drug, option=1)
  df_nu = df_nu.iloc[FIRST_ROW_AFTER_BURNIN:]
  col_names = list(df.columns)
  mdr_cases = col_names[:-5] # Last five columns are not about MDR
  mdr_cases.reverse() # reverse to plot most-dangerous type latest
  # Loop thru each MDR case
  for mdr_case in mdr_cases:
    color = coloring_legend(mdr_case, 1)
    ax.plot(df['time_elapsed'], df[mdr_case], 
            color=color)
    ax.fill_between(df_m['time_elapsed'], df_nl[mdr_case], df_nu[mdr_case], 
                    color=color, alpha=0.25)

def fig1_plot_vars(ax, dflist, drug):
  # Highest is 2-2 for DHA-PPQ
  if drug == 'DHA-PPQ':
    for df in dflist:
      df = df_col_replace(df, drug, option=1)
      df = df.iloc[FIRST_ROW_AFTER_BURNIN:]
      ax.plot(df['time_elapsed'], df['2-2'], color='#F88379', alpha=0.1)
  # Highest is 2-4 for ASAQ & AL
  else:
    for df in dflist:
      df = df_col_replace(df, drug, option=1)
      df = df.iloc[FIRST_ROW_AFTER_BURNIN:]
      ax.plot(df['time_elapsed'], df['2-4'], color='k', alpha=0.1)

def fig2_dangerous_triple(ax, df_l, df_m, df_u, df_ll, df_uu, pattern, annoty=None, ntf=None):
  # iloc gives data after burn-in only
  df_ll = df_ll.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_l = df_l.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_m = df_m.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_u = df_u.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_uu = df_uu.iloc[FIRST_ROW_AFTER_BURNIN:]
  ax.plot(df_m['time_elapsed'], df_m.filter(regex=pattern, axis=1).\
            sum(axis=1), color='#800080')
  ax.fill_between(df_m['time_elapsed'], 
                            df_l.filter(regex=pattern, axis=1).sum(axis=1), 
                            df_u.filter(regex=pattern, axis=1).sum(axis=1), 
                            color='#800080', alpha=0.25)
  ax.fill_between(df_m['time_elapsed'], 
                            df_ll.filter(regex=pattern, axis=1).sum(axis=1), 
                            df_uu.filter(regex=pattern, axis=1).sum(axis=1), 
                            color='#800080', alpha=0.1)
  
  # Mark where it exceeds the 10%
  try:      
    frn_tenp = df_m[df_m.filter(regex=pattern, axis=1).sum(axis=1).\
                gt(0.1)].index[0]
    ax.plot(df_m.loc[frn_tenp, 'time_elapsed'], 
            df_m.loc[frn_tenp].filter(regex=pattern).sum(), 'o', color='k')
  except:
    pass # Do not annotate if not reached

  try:
    # First Row Number that exceeds 1%
    frn = df_m[df_m.filter(regex=pattern, axis=1).sum(axis=1).\
                gt(0.01)].index[0]
    t_01 = df_m.loc[frn, 'time_elapsed']
    # Mark where it reaches 1% with dot
    ax.plot(t_01, df_m.loc[frn].filter(regex=pattern).sum(), 'o', color='k')
    t_01 = round(t_01/365-10, 1)
  except IndexError:
    t_01 = 'N/A'

  if annoty is not None:
    # for most_dange_trip type 1 
    # calculate genotype freq at last day
    x_20 = df_m.filter(regex=pattern, axis=1).sum(axis=1).tail(1).values[0]
    # calculate area under median curve, with IQR
    yaxis = df_m.filter(regex=pattern, axis=1).sum(axis=1).values
    xaxis = df_m['time_elapsed'].values
    auc = np.trapz(yaxis, x=xaxis)
    # lower quartile
    yaxis_l = df_l.filter(regex=pattern, axis=1).sum(axis=1).values
    xaxis_l = df_l['time_elapsed'].values
    auc_l = np.trapz(yaxis_l, x=xaxis_l)
    auc_l = round(auc_l, 2)
    # upper quartile
    yaxis_u = df_u.filter(regex=pattern, axis=1).sum(axis=1).values
    xaxis_u = df_u['time_elapsed'].values
    auc_u = np.trapz(yaxis_u, x=xaxis_u)
    auc_u = round(auc_u, 2)

    # annotate
    x_20 = round(x_20, 3)
    auc = round(auc, 2)
    annotation_string = r"$x_{20}$ = %s" % x_20
    annotation_string += "\n"
    annotation_string += r"$T_{.01}$ = %s" % t_01
    annotation_string += "\n"
    annotation_string += "AUC = %s (%s-%s)" % (auc, auc_l, auc_u)
    if ntf is not None:
      annotation_string += "\n"
      annotation_string += "NTF = %s" % ntf
    ax.text(ANNOTATION_X_LOCATION, annoty*0.65, annotation_string)

def fig2_double_and_higher(ax, df_l, df_m, df_u, df_ll, df_uu, drug, annoty=None):
  option = 1
  
  df_ll = df_col_replace(df_ll, drug, option)
  df_ll = df_ll.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_l = df_col_replace(df_l, drug, option)
  df_l = df_l.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_m = df_col_replace(df_m, drug, option)
  df_m = df_m.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_u = df_col_replace(df_u, drug, option)
  df_u = df_u.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_uu = df_col_replace(df_uu, drug, option)
  df_uu = df_uu.iloc[FIRST_ROW_AFTER_BURNIN:]

  col_names = list(df_m.columns)
  mdr_cases = col_names[:-5] # Last five columns are not about MDR
  for mdr_case in mdr_cases:
    # only plot double-or-higher resistant
    if mdr_case[0:1] == '2':
      color = coloring_legend(mdr_case,option)
      ax.plot(df_m['time_elapsed'], df_m[mdr_case], label=mdr_case, 
                      color=color)
      ax.fill_between(df_m['time_elapsed'], df_l[mdr_case], df_u[mdr_case], 
                      color=color, alpha=0.25)
      ax.fill_between(df_m['time_elapsed'], df_ll[mdr_case], df_uu[mdr_case], 
                      color=color, alpha=0.1)

  # get the most dangerous double type for annotation
  most_dang_type = mdr_cases[-1]
  # Mark where it exceeds the 10%
  try:      
    frn_tenp = df_m[df_m[most_dang_type].gt(0.1)].index[0]
    ax.plot(df_m.loc[frn_tenp, 'time_elapsed'], 
            df_m.loc[frn_tenp, most_dang_type], 'o', color='k')
  except:
    pass # Do not annotate if not reached

  try:
    # First Row Number that exceeds 1%
    frn = df_m[df_m[most_dang_type].gt(0.01)].index[0]
    t_01 = df_m.loc[frn, 'time_elapsed']
    # Mark where it reaches 1% with dot
    ax.plot(t_01, df_m.loc[frn, most_dang_type], 'o', color='k')
    t_01 = round(t_01/365-10, 1)
  except IndexError:
    t_01 = 'N/A'

  if annoty is not None:    
    # for most-dangerous double
    # calculate genotype freq at last day
    x_20 = df_m[most_dang_type].tail(1).values[0]
    # calculate area under median curve
    yaxis = df_m[most_dang_type].values
    xaxis = df_m['time_elapsed'].values
    auc = np.trapz(yaxis, x=xaxis)
    # lower quartile
    yaxis_l = df_l[most_dang_type].values
    xaxis_l = df_l['time_elapsed'].values
    auc_l = np.trapz(yaxis_l, x=xaxis_l)
    auc_l = round(auc_l, 2)
    # upper quartile
    yaxis_u = df_u[most_dang_type].values
    xaxis_u = df_u['time_elapsed'].values
    auc_u = np.trapz(yaxis_u, x=xaxis_u)
    auc_u = round(auc_u, 2)

    # annotate
    x_20 = round(x_20, 3)
    auc = round(auc, 2)
    annotation_string = r"$x_{20}$ = %s" % x_20
    annotation_string += "\n"
    annotation_string += r"$T_{.01}$ = %s" % t_01
    annotation_string += "\n"
    annotation_string += "AUC = %s (%s-%s)" % (auc, auc_l, auc_u)
    ax.text(ANNOTATION_X_LOCATION, annoty*0.65, annotation_string)