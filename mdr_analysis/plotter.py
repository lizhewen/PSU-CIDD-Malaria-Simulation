import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from plot_helper import xaxis_label_ticker, coloring_legend, df_col_replace

FIRST_ROW_AFTER_BURNIN = 120
ANNOTATION_X_LOCATION = 3833

def fig1_plot_IQR(ax, df_l, df_m, df_u, drug, annoty=None):
  df = df_col_replace(df_m, drug, option=1)
  df_nl = df_col_replace(df_l, drug, option=1)
  df_nu = df_col_replace(df_u, drug, option=1)
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

def fig1_plot_ten_vars(ax, dflist, drug):
  # Highest is 2-2 for DHA-PPQ
  if drug == 'DHA-PPQ':
    for df in dflist:
      df = df_col_replace(df, drug, option=1)
      ax.plot(df['time_elapsed'], df['2-2'], color='#F88379')
  # Highest is 2-4 for ASAQ & AL
  else:
    for df in dflist:
      df = df_col_replace(df, drug, option=1)
      ax.plot(df['time_elapsed'], df['2-4'], color='k')

def fig2_dangerous_triple(ax, df_l, df_m, df_u, pattern, annoty=None):
  # iloc gives data after burn-in only
  df_l = df_l.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_m = df_m.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_u = df_u.iloc[FIRST_ROW_AFTER_BURNIN:]
  ax.plot(df_m['time_elapsed'], df_m.filter(regex=pattern, axis=1).\
            sum(axis=1), color='#800080')
  ax.fill_between(df_m['time_elapsed'], 
                            df_l.filter(regex=pattern, axis=1).sum(axis=1), 
                            df_u.filter(regex=pattern, axis=1).sum(axis=1), 
                            color='#800080', alpha=0.25)
  if annoty is not None:
    # for most_dange_trip type 1 
    # calculate genotype freq at last day
    x_20 = df_m.filter(regex=pattern, axis=1).sum(axis=1).tail(1).values[0]
    # calculate time until it's 1% of total genotype freq
    # Get first row # that's bigger than threshold
    threshold = 0.01
    try:
      # First Row Number that exceeds the threshold
      frn = df_m[df_m.filter(regex=pattern, axis=1).sum(axis=1).\
                  gt(threshold)].index[0]
      t_01 = df_m.loc[frn, 'time_elapsed']
      t_01 = round(t_01/365-10, 1)
    except IndexError:
      t_01 = 'N/A'
    # calculate area under median curve
    yaxis = df_m.filter(regex=pattern, axis=1).sum(axis=1).values
    xaxis = df_m['time_elapsed'].values
    auc = np.trapz(yaxis, x=xaxis)

    # annotate
    x_20 = round(x_20, 3)
    auc = round(auc, 2)
    annotation_string = r"$x_{20}$ = %s" % x_20
    annotation_string += "\n"
    annotation_string += r"$T_{.01}$ = %s" % t_01
    annotation_string += "\n"
    annotation_string += "AUC = %s" % auc
    ax.text(ANNOTATION_X_LOCATION, annoty*0.65, annotation_string)

def fig2_double_and_higher(ax, df_l, df_m, df_u, drug, annoty=None):
  option = 1
  
  df_l = df_col_replace(df_l, drug, option)
  df_l = df_l.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_m = df_col_replace(df_m, drug, option)
  df_m = df_m.iloc[FIRST_ROW_AFTER_BURNIN:]
  df_u = df_col_replace(df_u, drug, option)
  df_u = df_u.iloc[FIRST_ROW_AFTER_BURNIN:]

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

  if annoty is not None:
    # get the most dangerous double type for annotation
    most_dang_type = mdr_cases[-1]
    # for most-dangerous double
    # calculate genotype freq at last day
    x_20 = df_m[most_dang_type].tail(1).values[0]
    # calculate time until it's 1% of total genotype freq
    # Get first row # that's bigger than threshold
    threshold = 0.01
    try:
      frn = df_m[df_m[most_dang_type].gt(threshold)].index[0]
      t_01 = df_m.loc[frn, 'time_elapsed']
      t_01 = round(t_01/365-10, 1)
    except IndexError:
      t_01 = 'N/A'
    # calculate area under median curve
    yaxis = df_m[most_dang_type].values
    xaxis = df_m['time_elapsed'].values
    auc = np.trapz(yaxis, x=xaxis)

    # annotate
    x_20 = round(x_20, 3)
    auc = round(auc, 2)
    annotation_string = r"$x_{20}$ = %s" % x_20
    annotation_string += "\n"
    annotation_string += r"$T_{.01}$ = %s" % t_01
    annotation_string += "\n"
    annotation_string += "AUC = %s" % auc
    ax.text(ANNOTATION_X_LOCATION, annoty*0.65, annotation_string)