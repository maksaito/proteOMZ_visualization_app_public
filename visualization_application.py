# coding: utf-8

import pandas as pd
import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool, PanTool, WheelZoomTool, BoxSelectTool, TapTool, OpenURL
from bokeh.models import GMapPlot, GMapOptions, Circle, DataRange1d, Range1d
from bokeh.io import curdoc
from bokeh.layouts import row, column, widgetbox, gridplot
from bokeh.models.widgets import Select, Slider, TextInput, DataTable, TableColumn, Div, Select
import itertools
import os

import logging
logger = logging.getLogger(__name__)

# data file locations
MAP_DATA_FILE = 'proteomz_stns.csv'
PROTEOMZ_DATA_FILE = 'ExampleDataset.csv'
TAXA_KEYS_DATA_FILE = 'Taxa_keys.csv'

# plot tools
TOOLS = "box_zoom, pan, xwheel_zoom, reset" 

# global visual parameters
TITLE_TEXT_SIZE = '20px'

# introduction text
INTRODUCTION_HTML = """<b>Ocean Proteomics data from <i>Falkor ProteOMZ Expedition</i> in January-February of 2016.
The prototype interactive display can explore millions of protein values from over a hundred 300 liter samples
collected in the Central Pacific Oxygen Minimum Zone to depths greater than 1 kilometer. Use the sliders and 
menus to enable scientific discovery within this novel dataset. <u>*NOTE: This is an example dataset containing <i>shuffled protein
annotations.</i> Public release of this dataset coming soon. </u> """

INTRODUCTION_WIDTH = 380
INTRODUCTION_HEIGHT = 130

# map visual parameters
MAP_WIDTH = 400
MAP_HEIGHT = 750
MAP_TITLE = 'ProteOMZ EXPEDITION FALKOR 2015'

MAP_LAT = 7.29
MAP_LON = -145.73
MAP_ZOOM = 4
MAP_API_KEY = 'AIzaSyDKbDUfDV3MhmMbUHVVClWdypnU4UrHVCA'
MAP_TYPE = 'hybrid'

DESELECTED_STATION_COLOR = 'white'
SELECTED_STATION_COLOR = 'red'

# profile visual parameters
PROFILE_TITLE = 'The Vertical Distribution of Microbial Proteins'
PROFILE_X_LABEL = 'Relative Abundance (Spectral Counts)'
PROFILE_Y_LABEL = 'Depth in the Ocean (meters)'
PROFILE_LINE_COLOR = 'red'
MAX_PROFILES = 1200
PROFILE_WIDTH = 600 
PROFILE_HEIGHT = 1100

# histogram visual parameters
HISTOGRAM_TITLE = 'All Spectra/IDs'
HISTOGRAM_X_LABEL = 'Sum of Proteins/Spectra'
HISTOGRAM_WIDTH = 400
HISTOGRAM_HEIGHT = 1100

# bar chart visual parameters
TAXA_BAR_TITLE = 'The Diversity of Microbial Proteins'
TAXA_BAR_WIDTH = 600
TAXA_BAR_HEIGHT = 350
TAXA_BAR_COLORS = ["#e6ab02", "#1f78b4", "#b2182b", "#7570b3", "#e7298a", "#66a61e", 
     "#d95f02", "#666666"] #, "#1b9e77"]

#table settings
TAXON_TABLE_WIDTH=600
TAXON_TABLE_HEIGHT=750

# initial selections
ALL = 'ALL'
INIT_TAXA_GROUP = ALL
INIT_EC_GROUP = ALL
INIT_PCTILE = 95
INIT_NUT = 'N+N'
INIT_PROT = 'P1'

ST_SELECT_TITLE = 'Station'
NUT_SELECT_TITLE = 'Select Hydrographic Parameter for Correlation'
TN_SELECT_TITLE = 'Select Microbial Taxon'
EC_SELECT_TITLE = 'Major Enzyme Classes'
PERCENTILE_SLIDER_TITLE = 'Percentile (Note: be patient below 90%)'

EC_GROUPS = ['Oxidoreductases','Transferases', 'Hydrolases', 'Lyases', 'Isomerases', 'Ligases']

# computing axis ranges

def compute_profile_axis_ranges(z, station_counts):
    # compute plot axis ranges for profile plot
    max_z, min_z = z.max(), z.min()
    min_c, max_c = 0, station_counts.max().max()
    return (max_z, min_z), (min_c, max_c)

def compute_histogram_axis_ranges(histogram_datasource):
    # compute plot axis ranges for histogram
    min_h = 0 
    max_h = max(histogram_datasource.data['prot_cts']) * 1.5
    return (min_h, max_h)

# main container
class Visualization(object):
    def __init__(self):
        """read data and construct plot elements in their initial state"""
        self.read_data()
        z, station_counts, hydrography_counts, all_counts, selected_nut = self.select_initial_data(self.stations[0])
        self.construct_datasources(z, station_counts, hydrography_counts, all_counts, selected_nut)
        # create plots and widgets
        self.make_plots(z, station_counts, hydrography_counts, selected_nut)
        self.make_widgets()

    def read_data(self):
        """read data and transform into dataframes"""
        self._read_map_data()
        self._read_proteomz_with_metadata()

    def _read_map_data(self):
        # second column data source for map stn/lat/long only, single point per stn
        self.stn_coor = pd.read_csv(MAP_DATA_FILE, index_col=None)

    def _read_proteomz_with_metadata(self):
        """read the large spreadsheet CSV, extract sections, and reorganize into
        meaningful dataframes"""
        df = pd.read_csv(PROTEOMZ_DATA_FILE, low_memory=False)

        # extract metadata section of spreadsheet containing station and depth information
        self.cruise_metadata = df[df.columns[:11]][:103]
        # stations are in that column of cruise_metadata
        self.stations = self.cruise_metadata.Station.unique().astype(int)

        # extract counts section of spreadsheet
        all_counts = df[df.columns[21:]][:103].transpose()
        self.all_counts = all_counts.dropna().astype(float)

        #extract hydrographic data
        hydrography = df[df.columns[4:17]][:103].transpose()
        self.hydrography = hydrography.dropna().astype(float)

        # extract metadata section of spreadsheet containing prot id information
        data = df[103:]
        data.index=data.pop('ID')
        for col in data.columns[:10]:
            data.pop(col)
        prot_metadata = data.transpose()

        ### For taxonomy information we read a different file
        taxa_df = pd.read_csv(TAXA_KEYS_DATA_FILE) 
        prot_metadata.best_hit_taxon_id = pd.to_numeric(prot_metadata['best_hit_taxon_id'], errors='coerce')
        taxa_df.best_hit_taxon_id = pd.to_numeric(taxa_df['best_hit_taxon_id'], errors='coerce')
        prot_metadata_taxa = pd.merge(prot_metadata, taxa_df, how='left')
        prot_metadata_taxa.index = prot_metadata.index
        prot_metadata = prot_metadata_taxa
        # get taxon groups from SelectGroup column
        prot_metadata.loc[prot_metadata['SelectGroup'].isnull(), 'SelectGroup'] = 'Unknown'

        # get ec group from ec column using substring
        self.ec_group = prot_metadata['EC'].dropna().str.slice(0,1).astype(int)
        self.unique_ec_groups = sorted(self.ec_group.unique())

        self.taxa_group = prot_metadata['SelectGroup']
        self.unique_taxa_groups = sorted(prot_metadata['SelectGroup'].unique())
        self.prot_metadata = prot_metadata

        self.selected_nut = hydrography

    def get_map_data(self, station=None):
        """get the dict needed to construct the ColumnDataSource for the map"""
        if station is not None: # optionally for a selected station
            self.stn_coor.iloc[:,-1] = DESELECTED_STATION_COLOR
            stnindex = (self.stn_coor['Stn'][self.stn_coor['Stn'] == station].index.values[0])
            self.stn_coor.iloc[stnindex, -1] = SELECTED_STATION_COLOR 
        return dict(stn=self.stn_coor['Stn'],
            long_single=self.stn_coor['Long'],
            lat_single=self.stn_coor['Lat'],
            color = self.stn_coor['Color'],
            z=self.stn_coor['Depth'])

    def get_profile_data(self, z, station_counts, hydrography_counts, all_counts):
        import itertools
        """given depths and counts, produce the dict necessary
        to construct the ColumnDataSource for the profile plot"""

        # construct the multiline
        zs, cs = [], []

        #add variables for scatter. hs is all hydrography data
        hs,ps = [],[]

        # add annotation and taxon for hover
        ann, tax = [], []

        # limit to maximum number of profiles to display
        sc_trunc = station_counts[:MAX_PROFILES]

        for prot in sc_trunc.index:
            cs.append(list(sc_trunc.loc[prot]))
            hs.append(list(hydrography_counts))
            zs.append(list(z))
            ann.append(self.prot_metadata.best_hit_annotation[prot])
            tax.append(self.prot_metadata.best_hit_species[prot])  

        return dict(zs=zs, cs=cs, ann=ann, tax=tax, hs=hs)

    def get_taxa_data(self, station_counts):
        """given counts, return the dict necessary to construct the
        ColumnDataSource for the taxa bar plot"""

        station_sums = station_counts.sum(axis=1)
        station_sums_df = station_sums.to_frame()
       
        newdf = pd.merge(station_sums_df, self.prot_metadata, left_index = True, right_index = True)
        newdf.columns.values[0] = 'spec_sum'

        station_sums = newdf.spec_sum
        select_group = newdf.SelectGroup 
        sum_group = pd.DataFrame({
            'spec_sum': station_sums,
            'select_group': select_group.fillna('Unknown')
        })
        
        fillTaxaZero = [0] * len(self.unique_taxa_groups)
        taxaFilldf = pd.DataFrame({'spec_sum' : fillTaxaZero, 'select_group' : self.unique_taxa_groups})

        sum_group = pd.merge(sum_group, taxaFilldf, how='outer')
        taxadf = sum_group.groupby(["select_group"], as_index=False).sum()
        
        select_group = taxadf.select_group.tolist()
        spec_sum = taxadf.spec_sum.tolist()

        return dict(select_group=select_group, spec_sum=spec_sum, barColors=TAXA_BAR_COLORS)

    def get_histogram_data(self, z, station_counts):
        """given depths and counts, produce the dict necessary to construct
        the ColumnDataSource for the histogram plot"""

        spec_cts= list(station_counts.sum(axis=0))
        prot_cts = list(station_counts.astype(bool).sum(axis=0))
        depth = z.tolist()
        sizes = list(map(lambda x: x * 0.002, spec_cts))
        return dict(spec_cts = spec_cts, prot_cts = prot_cts, depth = depth, sizes = sizes)

    def select_data(self, station, percentile, selected_taxa_group, selected_ec_group, selected_nut):
        """given user selections, return the real depths of this station along with
        the counts for this station above the percentile threshold"""

        # depth are the real depths of this station
        z = self.cruise_metadata[self.cruise_metadata.Station==station]['Real Depth']

        # counts for this station are the slice of columns corresponding
        # to the associated columns in the cruise_metadata df
        station_counts = self.all_counts[list(z.index)]
        all_counts = self.all_counts
        #return  hydrographic data for the station
        hydrography_counts = self.hydrography[list(z.index)]

        if selected_taxa_group != ALL:
            station_taxa_group = self.taxa_group[station_counts.index]
            selected_prots = station_taxa_group[station_taxa_group == selected_taxa_group]
            station_counts = station_counts.loc[selected_prots.index]

        if selected_ec_group != ALL:
            station_ec_group = self.ec_group[station_counts.index].dropna().astype(int)
            selected_prots = station_ec_group[station_ec_group == int(selected_ec_group)]
            station_counts = station_counts.loc[selected_prots.index]

        if selected_nut != ALL:
            station_nut_group = hydrography_counts.loc[selected_nut]
            hydrography_counts = station_nut_group

        # filter all proteins above certain count sum threshold
        s = station_counts.sum(axis=1)
        args = s.sort_values(ascending=False).index
        station_counts = station_counts.loc[args]
        station_counts = station_counts[s > s.quantile(percentile/100.)]
        return z, station_counts, hydrography_counts, all_counts, selected_nut
 
    def select_initial_data(self, init_station, init_pctile=INIT_PCTILE, init_taxa_group=INIT_TAXA_GROUP, init_ec_group=INIT_EC_GROUP, init_nut=INIT_NUT):
        """Get the depths and counts for initial user selections"""

        z, station_counts, hydrography_counts, all_counts, selected_nut = self.select_data(init_station, init_pctile, init_taxa_group, init_ec_group, init_nut)
        return z, station_counts, hydrography_counts, all_counts, selected_nut

    def construct_datasources(self, z, station_counts, hydrography_counts, all_counts, selected_nut):
        """construct plot ColumnDataSources based on user-selected depths/counts"""

        self.map_datasource = ColumnDataSource(self.get_map_data())
        self.profile_datasource = ColumnDataSource(self.get_profile_data(z, station_counts, hydrography_counts, all_counts))
        self.selected_nut = selected_nut
        self.histogram_datasource = ColumnDataSource(self.get_histogram_data(z, station_counts))
        self.taxa_datasource = ColumnDataSource(self.get_taxa_data(station_counts))

    def make_map(self):
        """construct the google map plot"""
        map_options = GMapOptions(lat=MAP_LAT, lng=MAP_LON, map_type=MAP_TYPE, zoom=MAP_ZOOM)
        plot = GMapPlot(
            x_range=Range1d(), y_range=Range1d(), map_options=map_options,
            plot_width = MAP_WIDTH, plot_height = MAP_HEIGHT
        )

        plot.title.text = MAP_TITLE
        plot.title.text_font_size = TITLE_TEXT_SIZE
        plot.api_key = MAP_API_KEY

        #show unselected stations in grey
        # second column data source for map stn/lat/long only, single point per stn
        circle2 = Circle(x="long_single", y="lat_single", fill_alpha=0.8, size=15, fill_color="color", line_color="black")
        plot.add_glyph(self.map_datasource, circle2)
        # had to convert to tap tool for simplicity, so we can only select one station at a time. V2 may be able to figure out how to select multiple stations using box select

        plot.add_tools(PanTool(), WheelZoomTool())

        return plot
 
    def make_profile(self, z, station_counts):
        """Make a figure with a multiline showing protein counts as profiles,
        given user selected depths and counts"""

        # figure with multiline
        TOOLS = "box_zoom, pan, xwheel_zoom, reset, tap" 
        (max_z, min_z), (min_c, max_c) = compute_profile_axis_ranges(z, station_counts)
        p = figure(plot_width=PROFILE_WIDTH, plot_height=PROFILE_HEIGHT, x_range=(min_c,max_c), y_range=(max_z,min_z),
            x_axis_label=PROFILE_X_LABEL, y_axis_label=PROFILE_Y_LABEL, title=PROFILE_TITLE, tools=TOOLS)
        p.title.text_font_size = TITLE_TEXT_SIZE
        p.multi_line(source=self.profile_datasource, xs='cs', ys='zs', line_width=2, line_alpha=0.5,
            line_color=PROFILE_LINE_COLOR, hover_line_alpha=1, hover_line_color=PROFILE_LINE_COLOR)

        hover = HoverTool(tooltips=[
                ('depth', '$y m'),
                ('count', '$x'),
                ('species', '@tax[$index]'),
                ('annotation', '@ann[$index]'),
            ])

        p.add_tools(hover)

        return p

    def make_correlation(self, z, station_counts):
        import itertools
        import unicodedata

        p = figure(plot_height=600, plot_width=400, title='Protein vs. Hydrographic Data')
        #plot as a multiline because list of lists but make lines invisible and markers visible
        #p.multi_line(source=self.profile_datasource, xs='hs', ys='cs', line_width=0)
        p.title.text_font_size = TITLE_TEXT_SIZE
        p.multi_line(xs='hs',ys='cs', source=self.profile_datasource, line_color=PROFILE_LINE_COLOR, hover_line_alpha=1, hover_line_color=PROFILE_LINE_COLOR)

        hover = HoverTool(tooltips=[
                ('depth', '$y m'),
                ('count', '$x'),
                ('species', '@tax[$index]'),
                ('annotation', '@ann[$index]'),
            ])
        p.xaxis.axis_label = 'Selected Hydrographic Data'
        p.yaxis.axis_label = 'Protein Abundance'
        p.add_tools(hover)

        return p

    def make_histogram(self, max_z):
        """Make a "histogram" showing number of proteins and number of spectra
        by depth"""

        h = figure(plot_width=400, plot_height=500, y_range=(max_z,0), title="Unique Proteins and PSMs",
            x_axis_label=HISTOGRAM_X_LABEL, tools=TOOLS)
        h.title.text_font_size = TITLE_TEXT_SIZE
        h.line("prot_cts", "depth", color = "black", legend = "# Proteins Identified", source=self.histogram_datasource)   
        h.circle("prot_cts", "depth", size="sizes", color="navy", alpha=0.5, legend = "# Spectra Identified", source=self.histogram_datasource)
        h.legend.location = "bottom_left"

        return h

    def make_taxa_bar(self):
        """Make a barchart of bins of taxonomic classes"""

        b = figure(x_range=self.unique_taxa_groups, plot_height=TAXA_BAR_HEIGHT, plot_width=TAXA_BAR_WIDTH, 
            title=TAXA_BAR_TITLE, toolbar_location="right", tools=TOOLS)
        b.title.text_font_size = TITLE_TEXT_SIZE
        b.vbar(x='select_group', top='spec_sum', width=0.9, fill_color='barColors',
               line_width=1, line_color='barColors', source=self.taxa_datasource)
        b.xaxis.major_label_text_font_size='12pt'
        b.xaxis.major_label_text_align='left'
        b.xaxis.major_label_orientation=-45
        b.left[0].formatter.use_scientific = False

        hover = HoverTool(tooltips=[
            ('taxa', '@select_group'),
            ('spectral counts', '@spec_sum'),
        ])

        b.add_tools(hover)

        return b

    def make_taxon_table(self):
        """Make the table showing proteins and taxa information"""

        # protein table
        columns = [
                TableColumn(field="ann", title="Protein Name",),
                TableColumn(field="tax", title="Taxon"),
                ]
        return DataTable(source=self.profile_datasource, columns=columns,
            width=TAXON_TABLE_WIDTH, height=TAXON_TABLE_HEIGHT)

    def make_plots(self, z, station_counts, hydrography_counts, selected_nut):
        """Construct the plots"""

        (max_z, _), (_, _) = compute_profile_axis_ranges(z, station_counts)
        self.map = self.make_map()
        self.profile = self.make_profile(z, station_counts)
        self.taxa_bar = self.make_taxa_bar()
        self.histogram = self.make_histogram(max_z)
        self.taxon_table = self.make_taxon_table()
        self.correlation = self.make_correlation(z, station_counts)

    def update_selected(self):
        """This is called when the user changes a selection"""

        # read selected values from widgets
        station = int(self.st_select.value)
        pct = self.percentile_slider.value
        selected_ec_group = self.ec_groups_by_name[self.ec_select.value]
        selected_taxa_group = self.tn_select.value
        selected_nut = self.nut_select.value

        # select data and return depths / counts
        z, station_counts, hydrography_counts, all_counts, selected_nut= self.select_data(station, pct, selected_taxa_group, selected_ec_group, selected_nut)

        # produce the data driving the plots
        profile_data = self.get_profile_data(z, station_counts, hydrography_counts, all_counts)
        histogram_data = self.get_histogram_data(z, station_counts)
        taxa_data = self.get_taxa_data(station_counts)
        map_data = self.get_map_data(station)

        # inject that data into the backing datasources, which will trigger the UI to update the plots
        self.profile_datasource.data = profile_data
        self.taxa_datasource.data = taxa_data
        self.map_datasource.data = map_data
        self.histogram_datasource.data = histogram_data
        self.nutrientname = selected_nut
        self.selected_nut= selected_nut

        # now compute plot axis ranges
        (max_z, min_z), (min_c, max_c) = compute_profile_axis_ranges(z, station_counts)
        (min_h, max_h) = compute_histogram_axis_ranges(self.histogram_datasource)

        # change them in the UI
        self.profile.y_range.start = max_z+100
        self.profile.y_range.end = 0
        self.profile.x_range.start = min_c
        self.profile.x_range.end = max_c

        self.histogram.y_range.start = max_z+100
        self.histogram.y_range.end = 0
        self.histogram.x_range.start = min_h
        self.histogram.x_range.end = max_h

    def make_widgets(self):
        """Construct the widgets and attach the update callback to them"""

        # station selector
        options = list(self.stations.astype(str))
        self.st_select = Select(title=ST_SELECT_TITLE, options=options, value=options[0])

        # taxon selector
        optionsTaxa = [ALL] + self.unique_taxa_groups
        self.tn_select = Select(title=TN_SELECT_TITLE, options=optionsTaxa, value=INIT_TAXA_GROUP)

        self.ec_groups_by_name = dict(zip(EC_GROUPS,range(1,len(EC_GROUPS)+1)))
        self.ec_groups_by_name.update({ALL:ALL})

        options = [ALL] + EC_GROUPS
        self.ec_select = Select(title=EC_SELECT_TITLE, options=options, value=INIT_EC_GROUP)

        # percentile slider
        self.percentile_slider = Slider(title=PERCENTILE_SLIDER_TITLE, value=INIT_PCTILE, start=1, end=100, step=1)

        def update_selected_callback(attr, old, new):
            self.update_selected()

        # select nutrient to correlate to
        nutrients = [ALL] + ['N+N','NH4','NO2','NO3','O2_2','PO4']
        self.nut_select = Select(title=NUT_SELECT_TITLE, options=nutrients, value=nutrients[0])

        # attach callback to widgets
        self.st_select.on_change('value', update_selected_callback)
        self.percentile_slider.on_change('value', update_selected_callback)
        self.ec_select.on_change('value', update_selected_callback)
        self.tn_select.on_change('value', update_selected_callback)
        self.nut_select.on_change('value', update_selected_callback)

    def make_layout(self):
        """construct the  layout"""
        # introductory Text
        introduction_div = Div(text=INTRODUCTION_HTML, width=INTRODUCTION_WIDTH, height=INTRODUCTION_HEIGHT)

        # create layout
        select_widgetbox = widgetbox(self.st_select, self.tn_select, self.ec_select, self.percentile_slider, self.nut_select)
        col1 = column(introduction_div, select_widgetbox, self.map)
        col2 = self.profile
        col3 = column(self.correlation, self.histogram)
        col4 = column(self.taxa_bar, widgetbox(self.taxon_table))
        return row(col1, col2, col3, col4)

# beginning of bokeh app

viz = Visualization()
curdoc().add_root(viz.make_layout())