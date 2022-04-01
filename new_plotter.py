import argparse
import utils
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pysam import VariantFile
import pandas as pd
import numpy as np
import pysam
import statistics
from random import randrange
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import typing
import matplotlib

"""
Global Settings
"""
STD_POINT_SIZE = 2
STD_ARROW_SIZE = 5
DEFAULT_BLUE = '#5683d6'
DEFAULT_GREY = '#a5a5a5'
DEFAULT_RED = '#ff5025'
SUBPLOT_TITLE_SIZE = 10
SMALL_TEXT_SIZE = 7
STD_COVERAGE_Y_LIM = (-6, 6)
CALL_BOX_ALPHA = .5
DUP_COLOR = '#1f78b4'
DEL_COLOR = '#b2df8a'
ULTRA_RARE_COLOR = '#1b9e77'
RARE_COLOR = '#d95f02'
COMMON_COLOR = '#7570b3'
REPEAT_COLORS = ['#fc8d62', '#8da0cb', '#e78ac3', '#a6d854']
DUP_ALPHA = 0.7
DEL_ALPHA = 0.3


def get_args():
    """
    Gather the command line parameters, returns a namespace
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--calls',
                        '-i',
                        dest='calls',
                        help='BGZIPED and TABIXED bed file of calls',
                        action='append')

    parser.add_argument('--region',
                        '-r',
                        dest='region',
                        required=True,
                        help='genomic location to plot in format 1:1000-2000')

    parser.add_argument('--window',
                        '-w',
                        dest='window',
                        type=int,
                        default=0,
                        help='window size to be plotted around the specified region (--region, -r)')

    parser.add_argument('--sample',
                        '-s',
                        dest='sample',
                        required=True,
                        help='sample identifier used to select the specific calls to plot')

    parser.add_argument('--coverage_scores',
                        '-c',
                        dest='coverage_scores',
                        help='Name of output file')

    parser.add_argument('--sites',
                        dest='sites',
                        help='BED file that hase been bgzipped and tabixed')

    parser.add_argument('--gnomad',
                        dest='gnomad',
                        help='gnomad SV BED file that hase been bgzipped and tabixed')

    parser.add_argument('--vardb',
                        dest='vardb',
                        help='vardb common variants BED file that hase been bgzipped and tabixed')

    parser.add_argument('--depth',
                        dest='depth',
                        help='bed file site name, mean and standard deviation depth coverages')

    parser.add_argument('--title',
                        dest='title',
                        default=None,
                        help='what to title the plot, "# Probes XX" will be appended to it.'
                             'Default will be SAMPLE_ID Chr:Start-End # Probes XX')

    parser.add_argument('--repeatmasker',
                        dest='repeatmasker',
                        help='repeat masker bed file that has been bgzipped and tabixed')

    parser.add_argument('--output',
                        '-o',
                        dest='output',
                        help='Name of output file')

    return parser.parse_args()


def ss(_line: str) -> typing.List[str]:
    """
    Strip new line characters and split the give line on tabs
    """
    return _line.strip().split('\t')


def parse_region(_region: utils.Interval, _window: int = 0) -> utils.Interval:
    """
    param region: string of genomic location example "1:1000-2000" or "X:1000-2000"
    param window: int window to include on either side, default 0
    """
    _chrom = _region.split(':')[0]
    _target_region = utils.Interval(
        chrom=_chrom,
        start=max(0,
                  int(_region.split(':')[1].split('-')[0]) - _window),
        end=int(_region.split(':')[1].split('-')[1]) + _window,
        data=None)
    return _target_region


def get_calls(_files: typing.List[str], _region: utils.Interval, _with_sample_id: str = None,
              _without_sample_id: str = None) -> typing.List[utils.Interval]:
    """
    params _files: list of gz files with an accompanying tbi file of calls
    param _region: Interval object
    param _with_sample_id: sample id to specifically include
    param _without_sample_id: sample id to specifically exclude
    return list of Interval objects, one for each call that intersects the region of interest
    """
    _calls = []
    for _f in _files:
        print(_f)
        # get all calls that overlap the region
        try:
            _tmp = utils.get_intervals_in_region(_region, _f)
        except ValueError:
            print('Warning: no overlapping calls in ' + _f)
            continue
        _f_name = _f.split('/')[-1].split('.')[0]
        [x.data.append(_f_name) for x in _tmp]
        _keepers = []
        # get calls with a specific sample
        if _with_sample_id is not None:
            _keepers = [x for x in _tmp if sum(_with_sample_id in y for y in x.data) > 0]
        # get calls without a specific sample
        if _without_sample_id is not None:
            _keepers = _keepers + [x for x in _tmp if sum(_without_sample_id not in y for y in x.data) > 0]
        _calls = _calls + _keepers
    return _calls


def add_interval_box(_ax: plt.Axes, _interval: utils.Interval, _color: str, _a: float = 0.25, _ymin: float = 0.0,
                     _ymax: float = 1.0, _hatch: str = None) -> None:
    """
    Mark an interval with a box on an axis
    param ax: matplotlib ax object the boxes should be added to
    param _interval: interval object to mark with a box
    param _color: string color to be used at the color parameter in matplotlib
    param _a: float 0.0-1.0 for the alpha value, default 0.25
    param _ymin: float 0.0-1.0 bottom start point for the box default 0
    param _ymax: float 0.0-1.0 top end point for the box default 1
    param _hatch: string, the hatch (patterned fill) to be used in the box
    return None
    """
    _ax.axvspan(_interval.start,
                _interval.end,
                _ymin, _ymax,
                alpha=_a,
                color=_color,
                hatch=_hatch
                )


def mark_calls(_ax: plt.Axes, _calls: typing.List[utils.Interval], _hatch_map: typing.Dict) -> None:
    """
    Mark all the _calls on the given _ax, color based on _color_map
    param _ax: matplotlib Axes object to plot on
    param _calls: list of intervals to be marked
    param _color_map: dictionary mapping file name to color
    param _hatch_map: dictionary mapping file names to hatch (patterned fill)
    return: None
    """
    # mark the calls on the top plot
    _bottom = 0
    _increment = 1.0 / len(_calls)
    for _c in _calls:
        _this_calls_color = None
        if 'Dup' in _c.data[0] or 'DUP' in _c.data[0] or 'dup' in _c.data[0]:
            _this_calls_color = DUP_COLOR
        elif 'Del' in _c.data[0] or 'DEL' in _c.data[0] or 'del' in _c.data[0]:
            _this_calls_color = DEL_COLOR
        else:
            _this_calls_color = DEFAULT_GREY
        _this_hatch = _hatch_map[_c.data[-1]]
        add_interval_box(_ax, _c, _this_calls_color, _a=.25, _ymin=_bottom, _ymax=_bottom + _increment,
                         _hatch=_this_hatch)
        _bottom = _bottom + _increment


def plot_sample_coverage_stats(_ax: plt.Axes, _region: utils.Interval, _scores_file: str,
                               _sample_id: str,
                               _sites: str = None) -> None:
    """
    param _ax: matplotlib axes object to plot onto
    param _ax_legend: matplotlib axes object put the legend in
    param _region: Interval object
    param _scores_file: bed file that has been bgzipped and tabixed and has normalized coverage stats
    param _sample_id: sample id
    param _sites: path to site/probe/exon BED file that have been bgzipped and tabixed
    """
    _sites_tabix = None
    if _sites is not None:
        _sites_tabix = pysam.TabixFile(_sites)

    _scores = utils.get_intervals_in_region(_region, _scores_file)

    _site_pos = []
    _means = []
    _stdevs = []
    _ys = []
    _non_samp_xs = []
    _non_samp_pos = []
    _samples = utils.get_header(_scores_file)[0].split('\t')[4:]
    _sample_column_index = _samples.index(_sample_id)
    _colors = []
    _color_legend = {}
    for _site in _scores:
        # the the color representing quality from the probe file
        _added = False
        if _sites_tabix is not None:
            res = _sites_tabix.fetch(_site.chrom, _site.start, _site.end)
            for _line in res:
                _colors.append(ss(_line)[3])
                if _colors[-1] not in _color_legend:
                    _color_legend[_colors[-1]] = ss(_line)[4]
                _added = True
                break
        # if there is not color information provided or there is no probe matching this region use the default color
        if _sites_tabix is None or not _added:
            _colors.append(DEFAULT_BLUE)
            _color_legend[DEFAULT_BLUE] = 'No info'
        _tmp_scores = [float(x) for x in _site.data[1:]]
        _site_pos.append(_site.end)
        _means.append(np.mean(_tmp_scores))
        _stdevs.append(np.std(_tmp_scores))
        _ys.append(_tmp_scores[_sample_column_index])
        _non_samp_xs += _tmp_scores
        _non_samp_pos += [_site.end] * len(_tmp_scores)

    # the non-sample points plotted in the background
    _non_sample_plot = _ax.plot(_non_samp_pos,
                                _non_samp_xs,
                                '-o',
                                lw=0,
                                markersize=STD_POINT_SIZE,
                                label='Cohort Sample Coverage',
                                c=DEFAULT_GREY, alpha=.3)

    # the mean of the cohort (not database)
    _mean_plot = None
    for i in range(len(_colors)):
        _tmp = _ax.plot(_site_pos[i],
                        _means[i],
                        marker='_',
                        markersize=5,
                        lw=0,
                        label='Mean Pop. Coverage',
                        c=_colors[i])
        if _mean_plot is None:
            _mean_plot = _tmp

    # plot the standard deviations of the cohort normalized coverage z-scores
    _std_plots = None
    for i in range(len(_means)):
        _tmp = _ax.plot([_site_pos[i], _site_pos[i]],
                        [_means[i] - _stdevs[i], _means[i] + _stdevs[i]],
                        lw=1,
                        color=_colors[i],
                        alpha=CALL_BOX_ALPHA)
        if _std_plots is None:
            _std_plots = _tmp

    # plot this samples normalized coverage z-scores as points
    _sample_cov_plot = None
    for i in range(len(_ys)):
        _tmp = _ax.plot(_site_pos,
                        _ys,
                        '-o',
                        lw=0,
                        markersize=STD_POINT_SIZE,
                        label='Sample Coverage',
                        c=_colors[i])
        if _sample_cov_plot is None:
            _sample_cov_plot = _tmp
    _ax.set_title('Normalized Coverage', loc='left', fontsize=SUBPLOT_TITLE_SIZE)

    # standardize the y-axis limits to match all other plots, only if it falls within that range
    if _ax.get_ylim()[0] >= STD_COVERAGE_Y_LIM[0] and _ax.get_ylim()[1] <= STD_COVERAGE_Y_LIM[1]:
        _ax.set_ylim(STD_COVERAGE_Y_LIM)
    remove_borders(_ax)
    remove_ticks(_ax, remove_y=False)
    _ax.tick_params(axis='both', which='both', labelsize=SMALL_TEXT_SIZE)
    _ax.set_ylabel('Z-score', size=SMALL_TEXT_SIZE)


def remove_borders(_ax: plt.Axes) -> None:
    """
    Remove all borders from the given matplotlib axes object
    param _ax: axes object
    """
    _ax.spines['top'].set_visible(False)
    _ax.spines['right'].set_visible(False)
    _ax.spines['bottom'].set_visible(False)
    _ax.spines['left'].set_visible(False)


def remove_ticks(_ax: plt.Axes, remove_x: bool = True, remove_y: bool = True) -> None:
    """
    Remove tick marks from a matplotlib axes object
    param _ax: axes object
    param remove_x=True: If true remove the x ticks
    param remove_y=True: If true remove the y ticks
    """
    if remove_x:
        _ax.get_xaxis().set_ticks([])
    if remove_y:
        _ax.get_yaxis().set_ticks([])


def plot_other_samples_calls(_ax: plt.Axes, _region: utils.Interval, _files: str, _sample_id: str,
                             _x_range: typing.Tuple[int, int], _color_map: typing.Dict,
                             _hatch_map: typing.Dict) -> None:
    """
    Plot boxes of the non-sample calls
    param _ax: matplotlib axes object to plot onto
    param _region: Interval object
    params _files: list of gz files with an accompanying tbi file of calls
    param _sample_id: sample id
    param _x_range: tuple contain the range of the x axis
    param _color_map: dictionary of calls files and colors they should be marked as
    param _hatch_map: dictionary mapping file names to hatch (patterned fill)
    """

    # find all calls that overlap the region that are not the current sample
    _calls = get_calls(_files, _region, _without_sample_id=_sample_id)
    _range_off_set = (_x_range[1] - _x_range[0]) * 0.01
    # plot each call
    _calls = list(_calls)
    _hatches = []
    _colors = []
    for _i, _c in enumerate(_calls):
        _this_calls_alpha = None
        _this_calls_color = DEFAULT_GREY
        _colors.append(DEFAULT_GREY)
        _hatches.append(None)
        if 'Dup' in _c.data[0] or 'DUP' in _c.data[0] or 'dup' in _c.data[0]:
            _this_calls_alpha = DUP_ALPHA
        elif 'Del' in _c.data[0] or 'DEL' in _c.data[0] or 'del' in _c.data[0]:
            _this_calls_alpha = DEL_ALPHA
        else:
            _this_calls_color = DEFAULT_GREY
            _colors.append(DEFAULT_GREY)
            _hatches.append(None)
    _Y, _Y_info = stack_plots(_ax, _calls, _hatches, _colors)
    _hatch_idx = 0
    _color_idx = 1
    for _i, _row in enumerate(_Y):
        for _j, _call in enumerate(_row):
            _ax.hlines(xmin=_call.start, xmax=_call.end, y=_i, color=_Y_info[_i][_j][_color_idx], linewidth=3, alpha=.5)
        # if the call extends past the x limits, add arrow
            if _call.start < _x_range[0]:
                # plot arrows if call goes out of range
                _ax.plot(_x_range[0] + _range_off_set, _i, '<', markersize=STD_ARROW_SIZE, c='black')
                _ax.plot(_x_range[0] + (_range_off_set * 2.5), _i, '>', markersize=STD_ARROW_SIZE, c='black')
            if _call.end > _x_range[1]:
                # plot arrows if call goes out of range
                _ax.plot(_x_range[1] - _range_off_set, _i, '>', markersize=STD_ARROW_SIZE, c='black')
                _ax.plot(_x_range[1] - (_range_off_set * 2.5), _i, '>', markersize=STD_ARROW_SIZE, c='black')

    # make the x_lim match the top plot
    _ax.set_xlim(_x_range)
    # add some extra space to the top of the plot for arrows are not cut off
    _lim = _ax.get_ylim()
    _ax.set_ylim(_lim[0] - 1, _lim[1] + 1)
    remove_borders(_ax)
    remove_ticks(_ax)
    _ax.set_title('Cohort Calls', loc='left', fontsize=SUBPLOT_TITLE_SIZE)


def plot_gnomad(_ax: plt.Axes, _ax_legend: plt.Axes, _region: utils.Interval, _file: str,
                _x_range: typing.Tuple[int, int]) -> None:
    """
    param _ax: matplotlib axes object to plot onto
    param _ax_legend: matplotlib axes object put the legend in
    param _region: Interval object
    param _file: gnomad sv file bgzipped and tabixed
    param x_range: tuple holding x axis limit
    """
    _gnomad_tbx = pysam.TabixFile(_file)
    _res = _gnomad_tbx.fetch(_region.chrom, _region.start, _region.end)

    _allele_freqs = []
    _regions = []
    _colors = []
    _hatches = []
    for _line in _res:
        _row = ss(_line)
        _regions.append(utils.Interval(chrom=_row[0], start=int(_row[1]), end=int(_row[2]), data=_row[3:]))
        # _colors.append(DEFAULT_GREY)
        _hatches.append(None)
        _afs = [float(_x) for _x in _row[37].split(',')]
        _af = None
        if len(_afs) > 1:
            _af = sum(_afs) / len(_afs)
        else:
            _af = _afs[0]

        if _af <= .01:
            # ultra rare
            _colors.append(ULTRA_RARE_COLOR)
        elif _af > .01 and _af < 0.05:
            # rare
            _colors.append(RARE_COLOR)
        else:
            # polymorphic
            _colors.append(COMMON_COLOR)

    _Y, _Y_info = stack_plots(_ax, _regions, _hatches, _colors)
    _hatch_idx = 0
    _color_idx = 1
    for _i, _row in enumerate(_Y):
        for _j, _call in enumerate(_row):
            _ax.hlines(xmin=_call.start, xmax=_call.end, y=_i, color=_Y_info[_i][_j][_color_idx], linewidth=3)

    _ax.set_xlim(_x_range)
    _ax.set_ylim(-1, len(_Y)-1)
    remove_ticks(_ax_legend)
    remove_borders(_ax_legend)
    remove_borders(_ax)
    remove_ticks(_ax)
    _ax.set_title('gnomAD SVs Allele Freq', loc='left', fontsize=SUBPLOT_TITLE_SIZE)


def plot_dbvar(_ax: plt.Axes, _ax_legend: plt.Axes, _region: utils.Interval, _file: str,
               _x_range: typing.Tuple[float, float]) -> None:
    """
    Plot VarDB common variants
    param _ax: matplotlib axes object to plot onto
    param _ax_legend: matplotlib axes object put the legend in
    param _region: Interval object
    param _file: VarDB Common Variants file bgzipped and tabixed
    param x_range: tuple holding x axis limit
    """

    _tbx = pysam.TabixFile(_file)
    _res = _tbx.fetch(_region.chrom, _region.start, _region.end)

    _regions = []
    _hatches = []
    _colors = []
    for _line in _res:
        _row = ss(_line)
        _regions.append(utils.Interval(chrom=_row[0], start=int(_row[1]), end=int(_row[2]), data=_row[3:]))
        _colors.append(DEFAULT_GREY)
        _hatches.append(None)

    # _ax.hlines(xmin=_x_starts, xmax=_x_ends, y=list(range(len(_x_starts))), color=DEFAULT_GREY, linewidth=3)
    _Y, _Y_info = stack_plots(_ax, _regions, _hatches, _colors)
    _hatch_idx = 0
    _color_idx = 1

    for _i, _row in enumerate(_Y):
        for _j, _call in enumerate(_row):
            _ax.hlines(xmin=_call.start, xmax=_call.end, y=_i, color=_Y_info[_i][_j][_color_idx], linewidth=3)

    _ax.set_xlim(_x_range)
    _ax.set_ylim(-1, len(_Y) - 1)
    remove_ticks(_ax_legend)
    remove_ticks(_ax)
    remove_borders(_ax_legend)
    remove_borders(_ax)
    _ax.set_title('dbVar Common Variants', loc='left', fontsize=SUBPLOT_TITLE_SIZE)


def plot_chr_mean_std(_ax_mean: plt.Axes, _ax_std: plt.Axes, _region: utils.Interval, _file: str) -> None:
    """
    param _ax: matplotlib axes object to plot onto
    param _ax_legend: matplotlib axes object put the legend in
    param _region: Interval object
    param _file: adjusted scores file bgzipped and tabixed
    """
    _xs = []
    _means = []
    _stds = []
    for _line in open(_file, 'r'):
        _row = ss(_line)
        if _row[0] != _region.chrom:
            continue
        _xs.append(int(_row[1]))
        _means.append(float(_row[4]))
        _stds.append(float(_row[5]))

    _ax_mean.scatter(_xs, _means, s=.1, color=DEFAULT_GREY)
    _ax_std.scatter(_xs, _stds, s=.1, color=DEFAULT_GREY)
    _ax_mean.set_title('Coverage', loc='left', fontsize=SUBPLOT_TITLE_SIZE)
    _ax_std.set_title('Coverage', loc='left', fontsize=SUBPLOT_TITLE_SIZE)
    _ax_std.set_ylabel('Std',size=SMALL_TEXT_SIZE)
    _ax_mean.set_ylabel('Mean',size=SMALL_TEXT_SIZE)
    _ax_std.tick_params(axis='both', which='both', labelsize=SMALL_TEXT_SIZE)
    _ax_mean.tick_params(axis='both', which='both', labelsize=SMALL_TEXT_SIZE)
    remove_borders(_ax_std)
    remove_borders(_ax_mean)
    remove_ticks(_ax_mean, remove_y=False)
    remove_ticks(_ax_std, remove_y=False)
    _ax_mean.axvspan(_region.start, _region.end, alpha=.5, color=DEFAULT_RED)
    _ax_std.axvspan(_region.start, _region.end, alpha=.5, color=DEFAULT_RED)
    _ax_std.set_xlabel('Chromosome ' + str(_region.chrom))
    # use only integer in y lim


def convert_size_to_string(size: int):
    """
    param size: number of base pairs
    return the size converted from an integer to a string e.g. 2000 -> 2Kb or 123,000,000 -> 123Mb
    """
    if size < 1000:
        return str(size) + ' bases'
    elif size >= 1000 and size < 1000000:
        return str(round(size / 1000,1)) + 'Kb'
    elif size >= 1000000 and size < 1000000000:
        return str(round(size / 1000000, 1)) + 'Mb'
    else:
        return str(round(size / 1000000000, 1)) + 'Gb'


def configure_title(_fig: plt.Figure, _region: utils.Interval, _scores_file: str, _sample_id: str, _title: str = None):
    """
    param _fig: matplotlib figure object for the whole entire plot
    param _region: Interval object
    param _scores_file: bed file that has been bgzipped and tabixed and has normalized coverage stats
    param _sample_id: sample id
    param _title: string what the plot should be titled, defaults to sample_id, region.
        # probes is appended to the end of all titles
    """
    _scores = utils.get_intervals_in_region(_region, _scores_file)
    if _title is None:
        _title = '{} {}:{}-{} # Probes {}'.format(_sample_id, str(_region.chrom), str(_region.start), str(_region.end),
                                                  str(len(_scores)))
        _title = _title + '\nTarget size: {}'.format(convert_size_to_string( _region.end - _region.start))
        _fig.suptitle(_title)
    else:
        _title = _title + ' # Probes {}'.format(str(len(_scores)))
        _title = _title + '\nTarget size: {}'.format(convert_size_to_string( _region.end - _region.start))
        _fig.suptitle(_title)

    # control the amount fo space below the title
    _fig.subplots_adjust(top=0.92)


def plot_probe_qualities(_ax: plt.Axes, _region: utils.Interval, _scores_file: str,
                         _sample_id: str,
                         _sites: str = None):
    """
    param _ax: matplotlib axes object to plot onto
    param _ax_legend: matplotlib axes object put the legend in
    param _region: Interval object
    param _scores_file: bed file that has been bgzipped and tabixed and has normalized coverage stats
    param _sample_id: sample id
    param _sites: path to site/probe/exon BED file that have been bgzipped and tabixed
    """

    _scores = utils.get_intervals_in_region(_region, _scores_file)
    _sites_tabix = pysam.TabixFile(_sites)
    _colors = []
    for _site in _scores:
        # the the color representing quality from the probe file
        _added = False
        if _sites_tabix is not None:
            res = _sites_tabix.fetch(_site.chrom, _site.start, _site.end)
            for _line in res:
                _colors.append(ss(_line)[3])

    _ax.set_title('Probe quality', loc='left', fontsize=SUBPLOT_TITLE_SIZE)
    for _i in range(len(_colors)):
        _ax.axvspan(_i, _i + .7, alpha=1.0, color=_colors[_i])
    _ax.set_ylim(0, 1)
    remove_ticks(_ax)
    remove_borders(_ax)


def a_single_legend_to_rule_them_all(_ax: plt.Axes, _region: utils.Interval, _scores_file: str, _hatch_map: typing.Dict,
                                     _sites: str = None) -> None:
    """
    param _ax: matplotlib Axes object to plot onto
    param _region: Interval object
    param _scores_file: bed file that has been bgzipped and tabixed and has normalized coverage stats
    param _sites: path to site/probe/exon BED file that have been bgzipped and tabixed
    """
    _sites_tabix = None
    if _sites is not None:
        _sites_tabix = pysam.TabixFile(_sites)
    _scores = utils.get_intervals_in_region(_region, _scores_file)

    _colors = []
    _color_legend = {}
    _color_legend['#ffffff'] = 'Probe Quality'
    for _site in _scores:
        # the the color representing quality from the probe file
        _added = False
        if _sites_tabix is not None:
            res = _sites_tabix.fetch(_site.chrom, _site.start, _site.end)
            for _line in res:
                _colors.append(ss(_line)[3])
                if _colors[-1] not in _color_legend:
                    _color_legend[_colors[-1]] = ss(_line)[4]
                _added = True
                break
    # Normalized coverage:
    # create the legend
    _legend_elements = []
    _legend_labels = []
    _color_legend[DEFAULT_GREY] = 'Cohort'
    for _color in _color_legend.keys():
        _legend_elements.append(Line2D([0], [0],
                                       marker='o',
                                       color=_color,
                                       # label=_color_legend[_color],
                                       markerfacecolor=_color,
                                       markersize=5,
                                       linewidth=0))
        _legend_labels.append(_color_legend[_color])

    # _legend_elements.append(Line2D([0], [0],
    #                                marker='_',
    #                                color=DEFAULT_GREY,
    #                                # label='Mean Pop.\nCoverage',
    #                                markerfacecolor=DEFAULT_GREY,
    #                                markersize=10,
    #                                linewidth=0))
    # _legend_labels.append('Mean Pop.\nCoverage')
    # _legend_elements.append(Line2D([0], [0],
    #                                marker='|',
    #                                color=DEFAULT_GREY,
    #                                # label='Std Pop.\nCoverage',
    #                                markerfacecolor=DEFAULT_GREY,
    #                                markersize=10,
    #                                linewidth=0))
    # _legend_labels.append('Std Pop.\nCoverage')

    """
    Call Groups
    """
    p = matplotlib.patches.Patch(facecolor='#ffffff', hatch=None, alpha=.5)
    _legend_labels.append('Call Groups')
    _legend_elements.append(p)
    for _key in _hatch_map.keys():
        p = matplotlib.patches.Patch(facecolor=DEFAULT_GREY, hatch=_hatch_map[_key], alpha=.5)
        _legend_labels.append(_key[:10])
        _legend_elements.append(p)

    _legend_elements.append(matplotlib.patches.Patch(facecolor='#ffffff', alpha=.5))
    _legend_labels.append('SV Type')
    _legend_elements.append(matplotlib.patches.Patch(facecolor=DEL_COLOR, alpha=.5))
    _legend_labels.append('Del')
    _legend_elements.append(matplotlib.patches.Patch(facecolor=DUP_COLOR, alpha=.5))
    _legend_labels.append('Dup')
    _legend_elements.append(matplotlib.patches.Patch(facecolor='#ffffff', alpha=.5))
    _legend_labels.append('Pop. Freq.')
    _legend_elements.append(matplotlib.patches.Patch(facecolor=ULTRA_RARE_COLOR))
    _legend_labels.append('Very Rare')
    _legend_elements.append(matplotlib.patches.Patch(facecolor=RARE_COLOR))
    _legend_labels.append('Rare')
    _legend_elements.append(matplotlib.patches.Patch(facecolor=COMMON_COLOR))
    _legend_labels.append('Polymorphic')
    # plot the legend
    _leg = _ax.legend(_legend_elements, _legend_labels,
                      frameon=False,
                      prop={'size': SMALL_TEXT_SIZE},
                      borderaxespad=0)

    remove_borders(_ax)
    remove_ticks(_ax)


def get_overlap(a, b):
    """
    param a: interval object
    param b: interval object
    return the amount of overlap two Interval objects have assuming they are on the same chromosome
    """
    return max(0, min(a.end, b.end) - max(a.start, b.start))


def stack_plots(_ax, _regions, _hatches, _colors):
    # list of lists, each internal list represent a row in the plot
    _Y = []
    _Y_info = []
    for _i, _call in enumerate(_regions):
        _added_to_a_row = False
        for _j, _row in enumerate(_Y):
            _fits_in_row = True
            for _r_call in _row:
                # check if the two intervals overlap
                if get_overlap(_r_call, _call) != 0:
                    _fits_in_row = False
                    break
            if _fits_in_row:
                _Y[_j].append(_call)
                _Y_info[_j].append((_hatches[_i], _colors[_i]))
                _added_to_a_row = True
                break
        # if it reaches here and has not been added to a row, start a new row
        if not _added_to_a_row:
            _Y.append([_call])
            _Y_info.append([(_hatches[_i], _colors[_i])])
    return _Y, _Y_info


def plot_repeat_masker(_ax, _ax_legend, _file, _region):
    _tbx = pysam.TabixFile(_file)
    _res = _tbx.fetch(_region.chrom, _region.start, _region.end)
    _repeat_element = []
    _segdup = []
    _selfchain = []
    _simple_repeat = []
    """
    repeatElement
    segdup
    selfchain
    simpleRepeat
    """
    for _line in _res:
        _row = ss(_line)
        _pair = (int(_row[1]), int(_row[2]))
        try:
            if 'repeatElement' in _row[4]:
                _repeat_element.append(_pair)
            elif 'segdup' in _row[4]:
                _segdup.append(_pair)
            elif 'selfchain' in _row[4]:
                _selfchain.append(_pair)
            elif 'simpleRepeat' in _row[4]:
                _simple_repeat.append(_pair)
            else:
                print('Unrecognized {}'.format(_row[4]))
        except IndexError:
            print(_row)
    _colors_index = 0
    for _pair in _repeat_element:
        _ax.hlines(xmin=_pair[0], xmax=_pair[1], y=_colors_index, color=REPEAT_COLORS[_colors_index], linewidth=3, alpha=.5)
    _colors_index = 1
    for _pair in _segdup:
        _ax.hlines(xmin=_pair[0], xmax=_pair[1], y=_colors_index, color=REPEAT_COLORS[_colors_index], linewidth=3,
                   alpha=.5)
    _colors_index = 2
    for _pair in _selfchain:
        _ax.hlines(xmin=_pair[0], xmax=_pair[1], y=_colors_index, color=REPEAT_COLORS[_colors_index], linewidth=3,
                   alpha=.5)
    _colors_index = 3
    for _pair in _simple_repeat:
        _ax.hlines(xmin=_pair[0], xmax=_pair[1], y=_colors_index, color=REPEAT_COLORS[_colors_index], linewidth=3,
                   alpha=.5)
    _ax.set_title('Repeat Masker', loc='left', fontsize=SUBPLOT_TITLE_SIZE)
    remove_ticks(_ax)
    remove_ticks(_ax_legend)
    remove_borders(_ax)
    remove_borders(_ax_legend)
    _ax.set_yticks([0,1,2,3])
    _ax.set_yticklabels(['Repeat Elements', 'SegDup', 'Self Chain', 'Simple Repeat'])
    _ax.tick_params(axis='both', which='both', labelsize=SMALL_TEXT_SIZE)


def main():
    """
    Main function, here just to ensure proper variable scopes
    """
    args = get_args()

    fig = plt.figure()
    fig.set_size_inches(8, 10, forward=True)
    # fig.patch.set_facecolor('#f6f6f6')

    """
    Format the axes
    """
    gs1 = GridSpec(8, 2, width_ratios=[4, 1], height_ratios=[12, 4, 4, 4, 4, 2, 4, 4], left=0.12, right=0.98,
                   wspace=0.10,
                   hspace=0.5)

    all_axes = []

    # top panel
    coverage_ax = fig.add_subplot(gs1[0, 0])
    # coverage_ax_legend = fig.add_subplot(gs1[0, 1])
    all_axes.append(coverage_ax)
    # all_axes.append(coverage_ax_legend)

    # middle panel
    ax2 = fig.add_subplot(gs1[1, 0])
    # ax2_legend = fig.add_subplot(gs1[1, 1])
    all_axes.append(ax2)
    # all_axes.append(ax2_legend)

    # legend
    first_legend = fig.add_subplot(gs1[:2, 1])
    all_axes.append(first_legend)

    # gnomad panel
    gnomad_ax = fig.add_subplot(gs1[2, 0])
    gnomad_ax_legend = fig.add_subplot(gs1[2, 1])
    all_axes.append(gnomad_ax)
    all_axes.append(gnomad_ax_legend)

    # vardb panel
    vardb_ax = fig.add_subplot(gs1[3, 0])
    vardb_ax_legend = fig.add_subplot(gs1[3, 1])
    all_axes.append(vardb_ax)
    all_axes.append(vardb_ax_legend)

    # vardb panel
    repeat_ax = fig.add_subplot(gs1[4, 0])
    repeat_ax_legend = fig.add_subplot(gs1[4, 1])
    all_axes.append(repeat_ax)
    all_axes.append(repeat_ax_legend)

    # probe quality report
    # site_qc_ax = fig.add_subplot(gs1[5, 0])
    # all_axes.append(site_qc_ax)

    # bottom panel
    mean_ax = fig.add_subplot(gs1[5, 0])
    std_ax = fig.add_subplot(gs1[6, 0])
    all_axes.append(mean_ax)
    all_axes.append(std_ax)

    """
    Plot the calls and a window size
    """
    target_window = parse_region(args.region, args.window)
    target = parse_region(args.region, 0)
    in_sample_calls = get_calls(args.calls, target_window, _with_sample_id=args.sample)
    # input_file_colors_options = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c']
    input_file_hatch_options = ['////', '....', '||||', '----', '++++', 'xxxx', 'oooo', '****', 'OOOO', '\\\\\\\\']
    # index x.data[-1] is the file the call came from
    f_names = list(set([x.data[-1] for x in in_sample_calls]))
    f_names.sort()
    # input_file_color_map = {y: input_file_colors_options[i] for i, y in enumerate(f_names)}
    input_file_hatch_map = {y: input_file_hatch_options[i] for i, y in enumerate(f_names)}


    """
    Mark this sample's coverage stats
    """
    plot_sample_coverage_stats(coverage_ax, target_window, args.coverage_scores, args.sample,
                               args.sites)
    mark_calls(coverage_ax, in_sample_calls, input_file_hatch_map)
    """
    Plot the calls made on other cohort members in this region
    """
    plot_other_samples_calls(_ax=ax2, _region=target_window, _files=args.calls,
                             _sample_id=args.sample, _x_range=coverage_ax.get_xlim(), _color_map=None,
                             _hatch_map=input_file_hatch_map)

    """
    Plot gnomAD
    """
    plot_gnomad(gnomad_ax, gnomad_ax_legend, _region=target_window, _file=args.gnomad, _x_range=coverage_ax.get_xlim())

    """
    Plot VarDB
    """
    plot_dbvar(vardb_ax, vardb_ax_legend, _region=target_window, _file=args.vardb, _x_range=coverage_ax.get_xlim())

    """
    Plot Mean and Standard Deviation Coverages across whole chromosome 
    """
    plot_chr_mean_std(mean_ax, std_ax, target_window, args.depth)

    """
    Make the title
    """
    configure_title(fig, target, args.coverage_scores, args.sample, args.title)

    """
    Add a non-to scale probe quality report
    """
    # plot_probe_qualities(site_qc_ax, target_window, args.coverage_scores, args.sample,
    #                      args.sites)

    """
    """
    a_single_legend_to_rule_them_all(first_legend, target_window, args.coverage_scores, input_file_hatch_map,
                                     args.sites)

    """
    Add the repeat masker lines
    """
    plot_repeat_masker(repeat_ax, repeat_ax_legend, _file=args.repeatmasker, _region=target_window)

    """
    Final Plotting Options
    """
    plt.savefig(args.output, bbox_inches='tight', dpi=600)
    plt.show()


if __name__ == '__main__':
    main()


"""
Huge Example
-i
ExampleData/calls6000.sorted.bed.gz
-i
ExampleData/calls_savvy1000.sorted.bed.gz
-o
test_figure_2.png
-r
15:2036300-22014250
-s
WES100_S13
-w
50000
--coverage_scores
ExampleData/adj_scores.bed.gz
--sites
ExampleData/wes_probes.sorted.bed.gz
--gnomad
ExampleData/gnomad_v2.1_sv.sites.bed.gz
--vardb
ExampleData/vardb.sorted.bed.gz
--depth
ExampleData/WES100_S13_probe.cover.mean.stdev.bed
--title
"WES100_S13 15:2036300-22014250 Duplication"
--repeatmasker
ExampleData/genomicRepeats.sorted.bed.gz
"""

"""
Small Example

-i
ExampleData/calls6000.sorted.bed.gz
-i
ExampleData/calls_savvy1000.sorted.bed.gz
-o
test_figure_17:18324000-18458000.png
-r
17:18324000-18458000
-s
WES100_S13
-w
50000
--coverage_scores
ExampleData/adj_scores.bed.gz
--sites
ExampleData/wes_probes.sorted.bed.gz
--gnomad
ExampleData/gnomad_v2.1_sv.sites.bed.gz
--vardb
ExampleData/vardb.sorted.bed.gz
--depth
ExampleData/WES100_S13_probe.cover.mean.stdev.bed
--title
"WES100_S13 17:18324000-18458000 Duplication"
--repeatmasker
ExampleData/genomicRepeats.sorted.bed.gz
"""