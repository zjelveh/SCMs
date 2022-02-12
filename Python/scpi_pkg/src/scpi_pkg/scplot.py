# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 07:55:26 2021

@author: Filippo Palomba
"""
# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas
import numpy
from plotnine import ggplot, aes, geom_point, geom_errorbar, geom_vline, geom_line, theme, theme_bw
from plotnine import element_blank, labs, guide_legend, scale_color_manual, ggtitle


def scplot(result,
           col_dots_t=None,
           col_line_t=None,
           col_dots_s=None,
           col_line_s=None,
           x_lab=None,
           y_lab=None,
           e_out=True,
           e_method=None,
           save_data=None):

    '''

    Parameters
    ----------

    result
    a class `scpi_est' object, obtained by calling scest, or a class
   `scpi_pi' object, obtained by calling scpi

    col_dots_t
    string indicating the color of the time series marker for treated unit

    col_line_t
    string indicating the color of the time series line for treated unit

    col_dots_s
    string indicating the color of the time series marker for synthetic control unit

    col_line_s
    string indicating the color of the time series line for synthetic control unit

    x_lab
    string indicating x axis title

    y_lab
    string indicating y axis title

    e_out
    a logical specifying whether out-of-sample uncertainty should be included in the plot(s).

    e_method
    a string specifying the type of uncertainty estimation used for out-of-sample uncertainty quantification.

    save_data
    a character specifying the name (and the folder) of the saved dataframe containing the processed data used to
    produce the plot. The data is saved in .csv format and the folder specified.

    Returns
    ----------
    plot plotnine object that can be further modified.

    References
    ----------
    Cattaneo, M. D., Feng, Y., and Titiunik, R. (2021), “Prediction Intervals for Synthetic Control
    Methods,” Journal of the American Statistical Association, 116, 1865-1880.

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2022), “Uncertainty Quantification in
    Synthetic Controls with Staggered Treatment Adoption,” working paper.

    See Also
    --------
    scdata, scest, scpi

    '''
    class_input = result.__class__.__name__

    if class_input not in ['scest_output', 'scpi_output']:
        raise Exception("The object 'result' should be the output of scest or scpi!")

    if x_lab is None:
        x_lab = 'Time'
    else:
        if not isinstance(x_lab, str):
            raise Exception("The option 'x_lab' should be a string!")

    if y_lab is None:
        y_lab = 'Outcome Variable'
    else:
        if not isinstance(y_lab, str):
            raise Exception("The option 'y_lab' should be a string!")

    period_pre = result.period_pre
    period_post = result.period_post
    T0 = period_pre[len(period_pre) - 1]

    if col_dots_t is None:
        col_dots_t = 'gray'
    if col_line_t is None:
        col_line_t = 'gray'
    if col_dots_s is None:
        col_dots_s = 'blue'
    if col_line_s is None:
        col_line_s = 'blue'

    # Store data for treated unit
    time = numpy.concatenate([period_pre, period_post])
    y_act = pandas.concat([result.Y_pre, result.Y_post]).to_numpy().flatten()
    data_points_act = pandas.DataFrame({'time': time,
                                        'y_act': y_act
                                        })

    # Store data for synthetic control unit
    y_sc_df = pandas.concat([result.Y_pre_fit, result.Y_post_fit])

    if class_input == 'scest_output':

        y_sc_na = pandas.DataFrame(numpy.array([numpy.nan] * len(time)))

        # Check if some periods have missing point estimate/missing CI
        not_miss_plot = [t in y_sc_df.index.tolist() for t in time]
        y_sc_na.loc[not_miss_plot, ] = y_sc_df.iloc[:, [0]].to_numpy()

        data_points_act = pandas.DataFrame({'time': time,
                                            'y_act': y_act
                                            })

        data_points = pandas.concat([data_points_act, y_sc_na], axis=1)
        data_points.columns = ['time', 'y_act', 'y_sc']
        data_points['tr'] = ['Treated'] * len(y_sc_na)
        data_points['sc'] = ['Synthetic Control'] * len(y_sc_na)

        plot_struc = (ggplot(data_points) +
                      theme_bw() +
                      theme(panel_grid=element_blank(),
                            legend_position=(.5, .05),
                            legend_direction='horizontal',
                            legend_box='horizontal',
                            legend_title=element_blank(),
                            subplots_adjust={'bottom': 0.2}) +
                      labs(x=x_lab, y=y_lab)
                      )

        plot = (plot_struc +
                geom_point(mapping=aes(x='time', y='y_act', color='tr'), shape='o', fill='white', na_rm=False)
                + geom_point(mapping=aes(x='time', y='y_sc', color='sc'), shape='o', na_rm=False) +
                geom_line(mapping=aes(x='time', y='y_act', color='tr'), na_rm=False) +
                geom_line(mapping=aes(x='time', y='y_sc', color='sc'), linetype='dashed', na_rm=False) +
                geom_vline(xintercept=T0, linetype='dotted') +
                scale_color_manual(name="", values=[col_dots_t, col_dots_s],
                                   labels=["Treated", "Synthetic Control"],
                                   guide=guide_legend(override_aes={'linetype': ['solid', 'dashed'],
                                                                    'shape': ['o', 'o']})))
        print(plot)

        return plot

    elif class_input == 'scpi_output':
        if e_method is None:
            e_method = result.e_method

        sc_l_0 = result.CI_in_sample.iloc[:, [0]].to_numpy()
        sc_r_0 = result.CI_in_sample.iloc[:, [1]].to_numpy()

        if e_method == 'gaussian':
            sc_l_1 = result.CI_all_gaussian.iloc[:, [0]].to_numpy()
            sc_r_1 = result.CI_all_gaussian.iloc[:, [1]].to_numpy()

        if e_method == 'ls':
            sc_l_1 = result.CI_all_ls.iloc[:, [0]].to_numpy()
            sc_r_1 = result.CI_all_ls.iloc[:, [1]].to_numpy()

        if e_method == 'qreg':
            sc_l_1 = result.CI_all_qreg.iloc[:, [0]].to_numpy(dtype='float64')
            sc_r_1 = result.CI_all_qreg.iloc[:, [1]].to_numpy(dtype='float64')

        y_sc_na = pandas.DataFrame(numpy.array([numpy.nan] * len(time)))
        sc_l_0_na = pandas.DataFrame(numpy.array([numpy.nan] * len(time)))
        sc_r_0_na = pandas.DataFrame(numpy.array([numpy.nan] * len(time)))
        sc_l_1_na = pandas.DataFrame(numpy.array([numpy.nan] * len(time)))
        sc_r_1_na = pandas.DataFrame(numpy.array([numpy.nan] * len(time)))

        # Check if some periods have missing point estimate/missing CI
        not_miss_plot = [t in y_sc_df.index.tolist() for t in time]
        not_miss_ci = [t in result.CI_in_sample.index.tolist() for t in time]

        y_sc_na.loc[not_miss_plot, ] = y_sc_df.iloc[:, [0]].to_numpy()
        sc_l_0_na.loc[not_miss_ci, ] = sc_l_0
        sc_r_0_na.loc[not_miss_ci, ] = sc_r_0
        sc_l_1_na.loc[not_miss_ci, ] = sc_l_1
        sc_r_1_na.loc[not_miss_ci, ] = sc_r_1

        data_points_act = pandas.DataFrame({'time': time,
                                            'y_act': y_act
                                            })

        data_points = pandas.concat([data_points_act, y_sc_na,
                                    sc_l_0_na, sc_r_0_na,
                                    sc_l_1_na, sc_r_1_na], axis=1)
        data_points.columns = ['time', 'y_act', 'y_sc', 'lb0', 'ub0', 'lb1', 'ub1']

        data_points['tr'] = ['Treated'] * len(y_sc_na)
        data_points['sc'] = ['Synthetic Control'] * len(y_sc_na)

        plot_struc = (ggplot(data_points) +
                      theme_bw() +
                      theme(panel_grid=element_blank(),
                            legend_position=(.5, .05),
                            legend_direction='horizontal',
                            legend_box='horizontal',
                            legend_title=element_blank(),
                            subplots_adjust={'bottom': 0.2}) +
                      labs(x=x_lab, y=y_lab)
                      )

        plot_lines = (plot_struc +
                      geom_point(mapping=aes(x='time', y='y_act', color='tr'), shape='o', fill='white', na_rm=False) +
                      geom_point(mapping=aes(x='time', y='y_sc', color='sc'), shape='o', na_rm=False) +
                      geom_line(mapping=aes(x='time', y='y_act', color='tr'), na_rm=False) +
                      geom_line(mapping=aes(x='time', y='y_sc', color='sc'), linetype='dashed', na_rm=False) +
                      geom_vline(xintercept=T0, linetype='dotted') +
                      scale_color_manual(name="", values=[col_dots_t, col_dots_s],
                                         labels=["Treated", "Synthetic Control"],
                                         guide=guide_legend(override_aes={'linetype': ['solid', 'dashed'],
                                                                          'shape': ['o', 'o']})))

        if e_out is False:
            plot = plot_lines + geom_errorbar(mapping=aes(x='time', ymin='lb0', ymax='ub0', color='sc'),
                                              size=0.5, linetype='solid') + ggtitle('In-sample Uncertainty')
            print(plot)

        if e_out is True:
            if e_method == 'gaussian':
                title_str = 'In and Out of Sample Uncertainty - Subgaussian Bounds'
                plot = plot_lines + geom_errorbar(mapping=aes(x='time', ymin='lb1', ymax='ub1', color='sc'),
                                                  size=0.5, linetype='solid') + ggtitle(title_str)
                print(plot)

            if e_method == 'ls':
                title_str = 'In and Out of Sample Uncertainty - Location-scale Model'
                plot = plot_lines + geom_errorbar(mapping=aes(x='time', ymin='lb1', ymax='ub1', color='sc'),
                                                  size=0.5, linetype='solid') + ggtitle(title_str)
                print(plot)

            if e_method == 'qreg':
                title_str = 'In and Out of Sample Uncertainty - Quantile Regression'
                plot = plot_lines + geom_errorbar(mapping=aes(x='time', ymin='lb1', ymax='ub1', color='sc'),
                                                  size=0.5, linetype='solid') + ggtitle(title_str)
                print(plot)

        # Save data to reproduce plot
        if save_data is not None:
            data_name = save_data + ".csv"
            data_points.to_csv(data_name)

        return plot
